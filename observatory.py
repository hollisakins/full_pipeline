from astropy.io import fits # fits module for opening and writing to fits files
from astropy import wcs # world coordinate system module for converting to RA/Dec
from astroquery.vizier import Vizier # for looking up stars in catalogs listed in Vizier
import astropy.coordinates as coord # for inputting coordinates into Vizier
import astropy.units as u # for units
import numpy as np 
import math 
from datetime import datetime
from time import strftime, gmtime, strptime, sleep
import os 
import warnings 
import sep 
import csv
import sys
import shutil
from collections import OrderedDict
from photutils import CircularAperture, CircularAnnulus, aperture_photometry, Background2D, MedianBackground
from astropy.stats import SigmaClip

# defines the width of the console to print output more clearly 
rows, columns = os.popen('stty size', 'r').read().split()
termsize = int(columns)

# astropy gives warning for a depricated date format in TheSkyX fits header,
# we dont need to see that so these two lines supress all warnings
# comment them out when testing
warnings.catch_warnings() 
warnings.simplefilter('ignore')

slow = False # if slow=True it will pause between printing each line just to make it easier to read

# function for printing the output consistently 
def prnt(indent,strng,filename=False):
    # if filename=True it will print the name of the file at the left
    if slow:
        if not filename:
            print(' '*len(indent+': ')+strng)
            sleep(0.3)
        else: 
            print(indent+': '+strng)
            sleep(0.3)
    else:
        if not filename:
            print(' '*len(indent+': ')+strng)
        else: 
            print(indent+': '+strng)

# simple function to print a line at the top of the screen that shows what process is going on 
def header(i):
    print('-'*int((termsize-len(i)-2)/2)+' '+i+' '+'-'*int((termsize-len(i)-2)/2))
    print('')


def makeMasters(writeOver=False):
    '''Index calibration files and generate masters. 

    Searches for calibration files in the most recent date directory under 'ArchCal/'. 
    For example, makeMasters() will search 'ArchCal/20180701/' if that is the most recent date directory available.

    Opt. Argument writeOver=False can be changed to True to allow older bias & dark frames to be
    over written by more recent ones.
    '''
    path_to_cal = 'ArchCal/'
    dates = [f for f in os.listdir(path_to_cal) if not f.startswith('.')] # index date folders in ArchCal
    path_to_cal += max(dates)+'/' # specify path as most recent date
    filenames = [f for f in os.listdir(path_to_cal) if os.path.isfile(os.path.join(path_to_cal,f)) if not f.startswith('.')] # list of filenames to process
    
    print('\033c')
    header('Making Masters')
    print('\tSearching %s for calibraton files...' % path_to_cal)
    print('\tIndexed %s files' % len(filenames))
    
    bias1,dark1,Red1,Green1,Blue1,R1,V1,B1,Halpha1,Lum1,filters1 = [],[],[],[],[],[],[],[],[],[],[] # initialize lists
    bias2,dark2,Red2,Green2,Blue2,R2,V2,B2,Halpha2,Lum2,filters2 = [],[],[],[],[],[],[],[],[],[],[] # initialize lists
    bias3,dark3,Red3,Green3,Blue3,R3,V3,B3,Halpha3,Lum3,filters3 = [],[],[],[],[],[],[],[],[],[],[] # initialize lists
    bias4,dark4,Red4,Green4,Blue4,R4,V4,B4,Halpha4,Lum4,filters4 = [],[],[],[],[],[],[],[],[],[],[] # initialize lists
    # lists are used to store the data for each calibration file and then combine into a master
    # each element in these lists will be a NxN numpy array of the fits data

    print('\tSorting files...')

    # sort the calibration images by type and store them in arrays
    for filename in filenames:
        with fits.open(path_to_cal+filename) as hdulist: 
            img = hdulist[0].data # split into data and header
            hdr = hdulist[0].header
            typ = hdr['IMAGETYP'] # save image type as variable
            binn = hdr['XBINNING'] # save binning as variable
            if typ=='Bias Frame':
                exec('bias'+str(binn)+'_header=hdr') # save the header to write back into the master
                exec('bias'+str(binn)+'.append(img)') # add the data to the list with respective type/binning
            if typ=='Dark Frame':
                exec('dark'+str(binn)+'_header=hdr')
                exec('dark'+str(binn)+'.append(img)')
            if typ=='Flat Field':
                exec(hdr['FILTER']+str(binn)+'_header=hdr')
                exec('filters'+str(binn)+".append(hdr['FILTER'])") # store the filters found in this directory in a list
                # so that we don't attempt to create new master flats with filters we did not have raw flats for
                exec(hdr['FILTER']+str(binn)+'.append(img)') 
    
    print('')
    print('\tIndexed files:        Binning1x1  Binning2x2  Binning3x3  Binning4x4')
    print('\t\tBias:             %s          %s          %s          %s' % (len(bias1),len(bias2),len(bias3),len(bias4)))
    print('\t\tDark:             %s          %s          %s          %s' % (len(dark1),len(dark2),len(dark3),len(dark4)))
    print('\t\tRed Flat:         %s          %s          %s          %s' % (len(Red1),len(Red2),len(Red3),len(Red4)))
    print('\t\tGreen Flat:       %s          %s          %s          %s' % (len(Green1),len(Green2),len(Green3),len(Green4)))
    print('\t\tBlue Flat:        %s          %s          %s          %s' % (len(Blue1),len(Blue2),len(Blue3),len(Blue4)))
    print('\t\tR Flat:           %s          %s          %s          %s' % (len(R1),len(R2),len(R3),len(R4)))
    print('\t\tV Flat:           %s          %s          %s          %s' % (len(V1),len(V2),len(V3),len(V4)))
    print('\t\tB Flat:           %s          %s          %s          %s' % (len(B1),len(B2),len(B3),len(B4)))
    print('\t\tHalpha Flat:      %s          %s          %s          %s' % (len(Halpha1),len(Halpha2),len(Halpha3),len(Halpha4)))
    print('\t\tLum Flat:         %s          %s          %s          %s' % (len(Lum1),len(Lum2),len(Lum3),len(Lum4)))
    print('')

    ## make the masters
    for i in ['1','2','3','4']: # for each binning factor 
        exec('s=np.size(bias'+i+')') # define var s as the size of the list
        if not s==0: # if the size is nonzero, there is data for that type & binning
            exec('bias'+i+'_master=np.median(np.array(bias'+i+'),axis=0)') # define bias master as the 
            print('\tConstructed a master bias with binning %sx%s' % (i,i))

    for i in ['1','2','3','4']:
        exec('s=np.size(dark'+i+')')
        if not s==0:
            try: # try to  make dark master
                exec('dark'+i+'_master=np.median(np.array(dark'+i+')-bias'+i+'_master,axis=0)') # make dark master by removing the bias first
                print('\tConstructed a scalable master dark with binning %sx%s' % (i,i))
            except NameError: # if you get a NameError:
                print('\tNo bias master for binning %sx%s, failed to create scalable dark. Wrote to DR_errorlog.txt' % (i,i))
                with open('DR_errorlog.txt','a') as erlog:
                    erlog.write('Failed to create scalable dark with binning %sx%s, no bias master present at'+strftime("%Y%m%d %H:%M GMT", gmtime()))

    for j in ['1','2','3','4']: 
        exec('f=np.unique(filters'+j+')') # establish unique filters 
        for i in f: # for each UNIQUE filter
            exec('s=np.size('+i+j+')')
            if not s==0: 
                exec(i+j+"_master = np.median("+i+j+",axis=0)/np.max(np.median("+i+j+",axis=0))")  # normalize flat field and make master
                print('\tConstructed master %s flat with binning %sx%s' % (i,j,j))
    

    # write the masters to fits files
    for i in ['1','2','3','4']:
        for j in ['bias','dark']: # for now: do not overwrite old bias / dark masters
            if j+i+'_master' in locals():
                try:
                    code = "fits.writeto('MasterCal/binning"+i+'/'+j+"_master.fit',"+j+i+'_master, header='+j+i+'_header,overwrite='+str(writeOver)+')'
                    exec(code)
                    print('\tWrote master %s to file MasterCal/binning%s/%s_master.fit' % (j,i,j))   
                except:
                    print('\tBias or dark master already exists, no new file written')

    for i in ['1','2','3','4']:
        exec('f=np.unique(filters'+i+')')
        for j in f: # only overwrite flats for the unique filters that we chose to update that night
            code = "fits.writeto('MasterCal/binning"+i+'/'+"flat_master_"+j+".fit',"+j+i+"_master,header="+j+i+"_header,overwrite=True)"
            exec(code)   
            print('\tWrote master %s flat to file MasterCal/binning%s/flat_master_%s.fit' % (j,i,j))
    
    print('\n\tComplete')
    sleep(3)
    print('\033c')



class Field:
    def __init__(self):
        # when a Field object is created, define some variables
        self.calibrated_path = 'Calibrated Images/' 
        self.uncalibrated_path = 'ArchSky/'
        self.path_to_cal = 'MasterCal/'
        self.calibrated = False
        self.aperture_size = 80.0
        self.max_temp = -3.0

    def openFits(self,filename,calibrated=False):
        self.filename = filename
        if not calibrated: # if it hasnt been calibrated we need the uncalibrated path 
            with fits.open(self.uncalibrated_path+self.filename) as hdulist:
                self.hdr = hdulist[0].header
                img = hdulist[0].data
                self.img = np.array(img,dtype='<f4') 
        else: # otherwise we need the calibrated path
            with fits.open(self.calibrated_path+self.filename.replace('.fits','_calibrated.fits')) as hdulist:
                self.hdr = hdulist[0].header
                img = hdulist[0].data
                self.img = np.array(img,dtype='<f4')

    def calibrate(self,h,image): # check to see whether or not we need to calibrate the file
        if np.size(image)==8487264 or np.size(image)==942748: # need to add number for 2x2 binning
            if h['CCD-TEMP']<=self.max_temp: 
                if h.get('CALSTAT',default=0)==0: 
                    return True # True means we calibrate
                if h.get('CALSTAT',default=0)=='BDF' or h.get('CALSTAT',default=0)=='DF':
                    return 'Redundant' # Redundant calibration
                if h.get('CALSTAT',default=0)=='D':
                    return 'OnlyDark' # already had an auto dark
            else:
                return 'Temp'    
        else:
            return 'Size'

    @staticmethod
    def write_to_header(head):
        if head.get('CALSTAT',default=0)==0: # if there is no calstat field in the header
            head.append(('CALSTAT','BDF','Status of Calibration')) # add one
        else:
            head['CALSTAT']='BDF' # otherwise set the value of calstat to BDF

    def save_file(self,head,data,filename):
        if not os.path.exists(self.calibrated_path): 
            os.makedirs(self.calibrated_path) # make a directory if there isnt one
        
        fits.writeto(self.calibrated_path+filename.replace(".fit","_calibrated.fit"),data,head,overwrite=True)
        prnt(self.filename,'Wrote file to '+self.calibrated_path)
        print(' ')
        self.calibrated = True # now its calibrated so we change this variable to True

    def initialize(self):
        '''Index the files we need to calibrate'''
        print("\033c") # clear the screen
        header('Initialization') # print the header
        self.columnsWritten = True # if we need to write the columns into the sources.csv file
        ## specify source files
        self.dates = [f for f in os.listdir(self.uncalibrated_path) if not f.startswith('.')] # index date folders in ArchSky
        self.uncalibrated_path += max(self.dates)+'/' # specify both paths as most recent date
        self.calibrated_path += max(self.dates)+'/' 
        if not os.path.exists(self.calibrated_path):
            os.makedirs(self.calibrated_path)
        all_files = [f for f in os.listdir(self.uncalibrated_path) if os.path.isfile(os.path.join(self.uncalibrated_path,f)) and not f.startswith('.')]
        self.list_of_files  = [f.rstrip() for f in all_files if f.endswith('.fits') or f.endswith('.fit')]

        src_files = [f for f in all_files if f.endswith('.SRC')]
        print('\tSearching %s for sky images...' % self.uncalibrated_path)
        sleep(0.3)
        print('\tSearching %s for calibration files...' % self.path_to_cal)
        sleep(0.3)
        print('\033c')


        for filename in src_files:
            shutil.copy(self.uncalibrated_path+filename, self.calibrated_path) # copy over SRC files since we need them later but obviously dont need to calibrate them


    def Reduce(self):
        print('\033c')
        header('Calibration & Source Extraction')
        
        light_h,light = self.hdr,self.img # bring up the hdr and image
        prnt(self.filename,'Successfully opened '+light_h['FILTER']+' image in '+self.uncalibrated_path,filename=True)
        self.path_to_cal = 'MasterCal/binning'+str(light_h['XBINNING'])+'/' # search for calibration files in binning-specific folder
     
        # open bias frame
        try: 
            bias_fits = fits.open(self.path_to_cal+'bias_master.fit') 
            prnt(self.filename,'Successfully opened bias master %s' % self.path_to_cal+'bias_master.fit')
        except: # if you encounter error
            prnt(self.filename,'Failed to open bias master %s' % self.path_to_cal+'bias_master.fit. Wrote to DR_errorlog.txt')
            with open('DR_errorlog.txt','a') as erlog: # open error log and write to it
                erlog.write('Missing bias master at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'. Auto DR halted.\n')
            sys.exit() # exit the program since you can't calibrate files without a bias frame

        bias_h = bias_fits[0].header # split into header and data
        bias = bias_fits[0].data

        # open dark frame
        try:
            dark_fits = fits.open(self.path_to_cal+'dark_master.fit') 
            prnt(self.filename,'Successfully opened dark master %s' % self.path_to_cal+'dark_master.fit')
        except:
            prnt(self.filename,'Failed to open dark master %s' % self.path_to_cal+'dark_master.fit. Wrote to DR_errorlog.txt')
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Missing dark master at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'. Auto DR halted.\n')
            sys.exit()

        dark_h = dark_fits[0].header
        dark = dark_fits[0].data
        
        dxptime = dark_h['EXPTIME'] # store the exposure time for the dark master for scaling purposes
        exptime = light_h['EXPTIME'] # store light image exposure time

        # open filter-specific flat field
        try: 
            flat_fits = fits.open(self.path_to_cal+'flat_master_'+light_h['FILTER']+'.fit') 
            prnt(self.filename,'Successfully opened '+self.path_to_cal+'flat_master_'+light_h['FILTER']+'.fit')
        except:
            prnt(self.filename,'Failed to open flat master %s' % self.path_to_cal+'flat_master_'+light_h['FILTER']+'.fit. Wrote to DR_errorlog.txt')
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Missing '+light_h['FILTER']+'flat master at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'. Auto DR halted.\n')
            sys.exit()
        
        flat_h = flat_fits[0].header
        flat = flat_fits[0].data

        # perform the actual data reduction
        if self.calibrate(light_h,light)==True: # if we need to calibrated
            prnt(self.filename,'Calibrating image...' )

            bias_corrected_image = light - bias # subtract the bias
            dark_corrected_image = bias_corrected_image - (exptime/dxptime)*dark # scale the dark linearly w/ exptime and subtract
            final_image = dark_corrected_image / flat # divide by the flat field (already normalized)
            
            self.write_to_header(light_h)
            self.save_file(light_h, final_image,self.filename)


        elif self.calibrate(light_h,light)=='OnlyDark': # if we only had an auto dark
            prnt(self.filename,'Calibrating image...' )
            
            final_image = light / flat # divide by the flat field

            self.write_to_header(light_h)
            self.save_file(light_h, final_image,self.filename)


        elif self.calibrate(light_h,light)=='Redundant': # if it was already calibrated
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Attempted redundant calibration on '+self.filename+' at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'\n')
            prnt(self.filename,'Image already calibrated')
            self.save_file(light_h, light,self.filename) # still save the file because we can still use it

        elif self.calibrate(light_h,light)=='Temp': # if its temp is wrong
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Image '+self.filename+' temp '+light_h['CCD-TEMP']+' degrees C, rejected calibration at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'.')
            prnt(self.filename,'Image taken at > '+str(self.max_temp)+' degrees C')

        elif self.calibrate(light_h,light)=='Size':
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Image '+self.filename+' not full size, rejected calibration at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'.')
            prnt(self.filename,'Image not full size')

        del flat_fits,bias_fits,dark_fits

    def Source(self): # gathers source extraction data from .SRC file
        src = np.loadtxt(self.calibrated_path+self.filename.replace('.fits','.SRC'))
        objects = src[:,0:2] # pixel X,Y coordinates of the objects in question
        X_pos = src[:,0]
        Y_pos = src[:,1]
        A = src[:,8]
        B = src[:,9]
        theta = src[:,10]
        prnt(self.filename,'Gathered source data')
        return {'obj':objects,'X':X_pos,'Y':Y_pos,'A':A,'B':B,'theta':theta} # return source data as dict

    def Convert(self): # converts obj list in pixel coordinate to RA-dec coordinates
        hdr = self.hdr
        w = wcs.WCS(hdr) # gets WCS matrix from the header
        objects = self.source['obj']
        world = w.wcs_pix2world(objects, 1) # World Coorindate System function converts matrix in fits header to RA/Dec
        prnt(self.filename,'Converted coordinates to RA/Dec')
        return world

    def Photometry(self):
        ### perform aperture photometry
        hdr,img = self.hdr,self.img
        egain = float(hdr['EGAIN'])

        bkg = sep.Background(img) # get background noise from image (maybe need to looki into issues with this?)
        img_sub = img - bkg # subtract background
        flux, fluxerr, flag = sep.sum_circle(img_sub, self.source['X'], self.source['Y'], self.aperture_size, err=bkg.globalrms)
        # get flux values from source extraction package
        instmag = -2.5*np.log10(flux) # convert flux to instrumental magnitude
        snr = np.sqrt(flux*egain)
        instmag_err = 1/snr

        for j in flux:
            if j<0:
                print('Negative')

        # snr = np.sqrt(final_sum*egain) 

        # instmag_err = 1/snr
        # instmag = -2.5*np.log10(final_sum)


        ### retrieve magnitudes from catalog
        time = hdr['DATE-OBS'] # time image was taken
        time = datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%f') # convert string to datetime object
        filt = hdr['FILTER']
        objects = self.world

        # lookup data in the UCAC4 catalog by querying Vizier
        v = Vizier(columns=['UCAC4','+_r','RAJ2000','DEJ2000','Bmag','Vmag','rmag'])
        output = OrderedDict([('id',[]),('RA_C',[]),('DEC_C',[]),('RA_M',[]),('DEC_M',[]),('DIF',[]),('MAG_R',[]),('MAG_V',[]),('MAG_B',[]),('MAG_err',[]),('CMAG_R',[]),('CMAG_V',[]),('CMAG_B',[]),('DATETIME',[]),('IMGNAME',[])])
        output['MAG_err'] = instmag_err
        cmags = [] # catalog magnitudes list
        misfires = 0 # number of errors
        
        objects_indices_matched = []

        for n in range(len(objects)):
            catalog = 'UCAC4' # specify catalog 
            #(if needed, we can implement a method to change this based on the object, which is why it is defined *inside* the loop)
            
            result = v.query_region(coord.SkyCoord(ra=objects[n,0], dec=objects[n,1],
            unit=(u.degree, u.degree), frame='fk5'),radius='2s',catalog=catalog) # submit query at object coordinates with a radius of 2 arcseconds
            
            try:
                # query_region returns result which is a TableList and we only need the first Table in the List
                result = result[0] # try to get the first result from the list of results (which is usually just 1 element)
                # but if there are NO tables in the list...
            except: # !! important, if we do not find a match we still save the data as this may be an anomoly or an object like an asteroid
                prnt(self.filename,'No star match within 2 arcseconds')
                misfires += 1 
                output['id'].append('nan')
                output['RA_C'].append('nan')
                output['DEC_C'].append('nan')
                output['RA_M'].append(objects[n,0])
                output['DEC_M'].append(objects[n,1])
                output['DIF'].append('nan')
                output['DATETIME'].append(time)
                output['IMGNAME'].append(self.filename)
                cmags.append('nan') 
                continue # skip to the next object by moving to the next iteration of the loop

            ids = np.array(result['UCAC4'],str) # get array of all the stars identified
            ra = np.array(result['RAJ2000'],float)
            dec = np.array(result['DEJ2000'],float) # catalog RA and Dec
            dif = np.array(result['_r'],float)
            
            fluxtype = filt+'mag' # get a variable for fluxtype to match to catalog magnitude types
            if filt=='R':
                fluxtype = 'rmag'
            
            flux = np.array(result[fluxtype],float)

            for i in range(len(ids)): # for all the stars matched, 
                if dif[i] <= 2 and i==np.argmin(dif) and ids[i] not in output['id']: # pick the star w the min residual value and less than 2 arcsec off
                    prnt(self.filename,'Star match in %s, mag %s, residual %s arcsec' % (catalog,flux[i],dif[i]))
                    output['id'].append(ids[i]) # add this data to the output dictionary 
                    output['RA_C'].append(ra[i])
                    output['DEC_C'].append(dec[i])
                    output['RA_M'].append(objects[n,0])
                    output['DEC_M'].append(objects[n,1])
                    output['DIF'].append(dif[i])
                    output['DATETIME'].append(time)
                    output['IMGNAME'].append(self.filename)
                    cmags.append(flux[i])
                    objects_indices_matched.append(n)

                else:
                    prnt(self.filename,'No star match within 2 arcseconds')
                    misfires += 1
                    output['id'].append('nan')
                    output['RA_C'].append('nan')
                    output['DEC_C'].append('nan')
                    output['RA_M'].append(objects[n,0])
                    output['DEC_M'].append(objects[n,1])
                    output['DIF'].append('nan')
                    output['DATETIME'].append(time)
                    output['IMGNAME'].append(self.filename)
                    cmags.append('nan') 
                    continue


        prnt(self.filename,'Output %s stars' % len(output['id']))
        prnt(self.filename,'Output %s unique stars' % len(set(output['id'])))
        prnt(self.filename,'Missed %s objects' % misfires)
        

        instmags_to_median = [instmag[m] for m in objects_indices_matched]
        cmags_nonan = [k for k in cmags if not math.isnan(float(k))] 
        
        if not len(instmags_to_median)==len(cmags_nonan):
            raise Exception('Catalog comparison list not same length as instrumental magnitude list')

        # with open(self.filename+'.csv', 'a') as outfile:
        #     writer = csv.writer(outfile)
        #     writer.writerows(magnitudes)

        # median_i = np.median(np.array(instmags_to_median))
        # median_c = np.median(np.array(cmags_nonan)) # median of catalog magnitude values
        d = np.array(cmags_nonan) - np.array(instmags_to_median)
        d = float(np.median(d))

        # d = float(median_c) - float(median_i) # difference is the difference between the medians


        for i in ['R','V','B']:
            magtype = 'MAG_'+i
            if i==filt:
                output[magtype] = instmag
            else:
                output[magtype] = np.full(np.shape(instmag),'---',dtype="S3")

        for i in ['R','V','B']:
            magtype = 'CMAG_'+i
            if i==filt:
                output[magtype] = cmags
            else:
                output[magtype] = np.full(np.shape(cmags),'---',dtype="S3")

        output['MAG_'+filt] += d # offset magnitudes by that difference 
        
        prnt(self.filename,'Wrote magnitude data to sources.csv')
        sleep(3)
        print(' ')
        sleep(1)
        print("\033c")
        header('Calibration & Source Extraction')
        
        return output

    def Plot(self):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse
        import matplotlib
        
        hdr, img = self.hdr, self.img
        proj = wcs.WCS(hdr)
        fig = plt.figure(figsize=(13,10)) 
        ax = fig.add_subplot(111,projection=proj)
        m, s = np.mean(img), np.std(img)
        im = ax.imshow(img, interpolation='nearest', cmap='gray',
                    vmin=m-s, vmax=m+s, origin='lower')

        overlay = ax.get_coords_overlay('fk5')
        overlay.grid(color='white', ls='dotted')
        overlay[0].set_axislabel('Right Ascension (J2000)')
        overlay[1].set_axislabel('Declination (J2000)')

        # plot an ellipse for each object
        for i in range(len(self.source['X'])):
            e = Ellipse(xy=(self.source['X'][i], self.source['Y'][i]),
                        width=6*self.source['A'][i],
                        height=6*self.source['B'][i],
                        angle=self.source['theta'][i])
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        name = self.filename.replace('Calibrated Images/','plots/')
        name = name.replace('.fits','.jpg')
        plt.savefig(name)
        print('    Created plot of %s' % self.filename)



    def writeData(self):
        with open('sources.csv', 'a') as outfile:
            writer = csv.writer(outfile)
            if self.columnsWritten==False:
                writer.writerow(self.output.keys())
                self.columnsWritten = True
            writer.writerows(zip(*self.output.values()))
        print('\033c')
   
    def Extract(self):
        self.source = self.Source()
        self.world = self.Convert()
        self.output = self.Photometry()
        self.writeData()

