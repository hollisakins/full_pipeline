from astropy.io import fits # fits module for opening and writing to fits files
from astropy import wcs # world coordinate system module for converting to RA/Dec
from astroquery.vizier import Vizier # for looking up stars in catalogs listed in Vizier
import astropy.coordinates as coord # for inputting coordinates into Vizier
import astropy.units as u # for units
import numpy as np 
import math 
from datetime import datetime
from time import strftime, gmtime, strptime, sleep
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os 
import warnings 
import sep 
import csv
import sys
import shutil

rows, columns = os.popen('stty size', 'r').read().split()
termsize = int(columns)

warnings.catch_warnings() 
warnings.simplefilter('ignore')

slow = False
def prnt(indent,strng,filename=False):
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


############ implement binning into this
def makeMasters(path_to_cal,writeOver=False):
    dates = [f for f in os.listdir(path_to_cal) if not f.startswith('.')] # index date folders in ArchCal
    path_to_cal += max(dates)+'/' # specify path as most recent date
    filenames = [f for f in os.listdir(path_to_cal) if not f.startswith('.')] # list of filenames to process
    print('Searching %s for calibraton files...' % path_to_cal)

    bias1,dark1,Red1,Green1,Blue1,R1,V1,B1,Halpha1,Lum1,filters1 = [],[],[],[],[],[],[],[],[],[],[] # initialize lists
    bias2,dark2,Red2,Green2,Blue2,R2,V2,B2,Halpha2,Lum2,filters2 = [],[],[],[],[],[],[],[],[],[],[] # initialize lists
    bias3,dark3,Red3,Green3,Blue3,R3,V3,B3,Halpha3,Lum3,filters3 = [],[],[],[],[],[],[],[],[],[],[] # initialize lists

    # lists are used to store the data for each calibration file and then combine into a master
    ## sort the calibration images by type and store them in arrays
    for filename in filenames:
        with fits.open(path_to_cal+filename) as hdulist:
            img = hdulist[0].data # split into data and header
            hdr = hdulist[0].header
            typ = hdr['IMAGETYP']
            binn = hdr['XBINNING']
            if typ=='Bias Frame':
                exec('bias'+str(binn)+'_header=hdr')
                exec('bias'+str(binn)+'.append(img)')
            if typ=='Dark Frame':
                exec('dark'+str(binn)+'_header=hdr')
                exec('dark'+str(binn)+'.append(img)')
            if typ=='Flat Field':
                exec(hdr['FILTER']+str(binn)+'_header=hdr')
                exec('filters'+str(binn)+".append(hdr['FILTER'])") # store the filters found in this directory in a list
                # so that we don't attempt to create new master flats with filters we did not have raw flats for
                exec(hdr['FILTER']+str(binn)+'.append(img)') # string operations to add data to filter-specific list
    print('')
    print('Indexed files:        Binning1x1  Binning2x2  Binning3x3')
    print('\tBias:             %s          %s          %s' % (len(bias1),len(bias2),len(bias3)))
    print('\tDark:             %s          %s          %s' % (len(dark1),len(dark2),len(dark3)))
    print('\tRed Flat:         %s          %s          %s' % (len(Red1),len(Red2),len(Red3)))
    print('\tGreen Flat:       %s          %s          %s' % (len(Green1),len(Green2),len(Green3)))
    print('\tBlue Flat:        %s          %s          %s' % (len(Blue1),len(Blue2),len(Blue3)))
    print('\tR Flat:           %s          %s          %s' % (len(R1),len(R2),len(R3)))
    print('\tV Flat:           %s          %s          %s' % (len(V1),len(V2),len(V3)))
    print('\tB Flat:           %s          %s          %s' % (len(B1),len(B2),len(B3)))
    print('\tHalpha Flat:      %s          %s          %s' % (len(Halpha1),len(Halpha2),len(Halpha3)))
    print('\tLum Flat:         %s          %s          %s' % (len(Lum1),len(Lum2),len(Lum3)))
    print('')

    ## make the masters
    for i in ['1','2','3']:
        exec('s=np.size(bias'+i+')')
        if not s==0:
            exec('bias'+i+'_master=np.median(np.array(bias'+i+'),axis=0)')
            print('Constructed a master bias with binning %sx%s' % (i,i))

    for i in ['1','2','3']:
        exec('s=np.size(dark'+i+')')
        if not s==0:
            try:
                exec('dark'+i+'_master=np.median(np.array(dark'+i+')-bias'+i+'_master,axis=0)')
                print('Constructed a scalable master dark with binning %sx%s' % (i,i))
            except NameError:
                print('* No bias master for binning %sx%s, failed to create scalable dark. Wrote to DR_errorlog.txt' % (i,i))
                with open('DR_errorlog.txt','a') as erlog:
                    erlog.write('Failed to create scalable dark with binning %sx%s, no bias master present at'+strftime("%Y%m%d %H:%M GMT", gmtime()))


    for j in ['1','2','3']:
        exec('f=np.unique(filters'+j+')')
        for i in f: # for each UNIQUE filter
            exec('s=np.size('+i+j+')')
            if not s==0:
                exec(i+j+"_master = np.median("+i+j+",axis=0)/np.max(np.median("+i+j+",axis=0))")  # more string operations
                # normalize flat field
                print('Constructed master %s flat with binning %sx%s' % (i,j,j))
    ## write the masters to fits files
    for i in ['1','2','3']:
        for j in ['bias','dark']: # for now: do not overwrite old bias / dark masters
            if j+i+'_master' in locals():
                try:
                    code = "fits.writeto('MasterCal/binning"+i+'/'+j+"_master.fit',"+j+i+'_master, header='+j+i+'_header,overwrite='+str(writeOver)+')'
                    exec(code)
                    print('Wrote master %s to file MasterCal/binning%s/%s_master.fit' % (j,i,j))   
                except:
                    print('Bias or dark master already exists, no new file written')

    for i in ['1','2','3']:
        exec('f=np.unique(filters'+i+')')
        for j in f: # only overwrite flats for the unique filters that we chose to update that night
            code = "fits.writeto('MasterCal/binning"+i+'/'+"flat_master_"+j+".fit',"+j+i+"_master,header="+j+i+"_header,overwrite=True)"
            exec(code)   
            print('Wrote master %s flat to file MasterCal/binning%s/flat_master_%s.fit' % (j,i,j))



class Field:
    def __init__(self):
        self.calibrated_path = 'Calibrated Images/'
        self.uncalibrated_path = 'ArchSky/'
        self.path_to_cal = 'MasterCal/'
        self.c = False
        self.aperture_size = 10.0

    def openFits(self,filename,calibrated=False):
        self.filename = filename
        if not calibrated: # if it hasnt been calibrated we need the uncalibrated path 
            with fits.open(self.uncalibrated_path+self.filename) as hdulist:
                self.hdr = hdulist[0].header
                img = hdulist[0].data
                self.img = np.array(img,dtype='<f4')
        else:
            with fits.open(self.calibrated_path+self.filename.replace('.fits','_calibrated.fits')) as hdulist:
                self.hdr = hdulist[0].header
                img = hdulist[0].data
                self.img = np.array(img,dtype='<f4')

    @staticmethod
    def calibrate(h,image):
        if np.size(image)==8487264 or np.size(image)==942748: # need to add number for 2x2 binning
            if h['CCD-TEMP']<=-3.0:
                if h.get('CALSTAT',default=0)==0: 
                    return True
                if h.get('CALSTAT',default=0)=='BDF' or h.get('CALSTAT',default=0)=='DF':
                    return 'Redundant' 
                if h.get('CALSTAT',default=0)=='D':
                    return 'OnlyDark'
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

    def save_file(self,head,data,day,filename):
        if not os.path.exists('Calibrated Images/'+day): 
            os.makedirs('Calibrated Images/'+day) 
        
        fits.writeto('Calibrated Images/'+day+'/'+filename.replace(".fit","_calibrated.fit"),data,head,overwrite=True)
        prnt(self.filename,'Wrote file to Calibrated Images/'+day)
        print(' ')
        self.c = True

    def initialize(self):
        '''Index the files we need to calibrate'''
        print("\033c")
        print('-'*termsize)
        self.columnsWritten = True
        ## specify source files
        self.dates = [f for f in os.listdir(self.uncalibrated_path) if not f.startswith('.')] # index date folders in ArchSky
        self.uncalibrated_path += max(self.dates)+'/' # specify path as most recent date
        self.calibrated_path += max(self.dates)+'/'
        if not os.path.exists(self.calibrated_path):
            os.makedirs(self.calibrated_path)
        all_files = [f for f in os.listdir(self.uncalibrated_path) if os.path.isfile(os.path.join(self.uncalibrated_path,f)) and not f.startswith('.')]
        self.list_of_files  = [f for f in all_files if not f.endswith('.SRC')]
        src_files = [f for f in all_files if f.endswith('.SRC')]

        print('Searching %s for sky images...' % self.uncalibrated_path)
        print('Searching %s for calibration files...' % self.path_to_cal)

        for filename in src_files:
            shutil.copy(self.uncalibrated_path+filename, self.calibrated_path)


    def Reduce(self):                
        light_h,light = self.hdr,self.img
        prnt(self.filename,'Successfully opened '+light_h['FILTER']+' image in '+self.uncalibrated_path,filename=True)
        self.path_to_cal = 'MasterCal/'
        self.path_to_cal += 'binning'+str(light_h['XBINNING'])+'/'
        try: 
            self.bias_fits = fits.open(self.path_to_cal+'bias_master.fit') 
            prnt(self.filename,'Successfully opened bias master %s' % self.path_to_cal+'bias_master.fit')
        except: # if you encounter error
            prnt(self.filename,'Failed to open bias master %s' % self.path_to_cal+'bias_master.fit. Wrote to DR_errorlog.txt')
            with open('DR_errorlog.txt','a') as erlog: # open error log and write to it
                erlog.write('Missing bias master at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'. Auto DR halted.\n')
            sys.exit() # exit the program since you can't calibrate files without a bias frame

        self.bias_h = self.bias_fits[0].header # split into header and data
        self.bias = self.bias_fits[0].data


        try:
            self.dark_fits = fits.open(self.path_to_cal+'dark_master.fit') 
            prnt(self.filename,'Successfully opened dark master %s' % self.path_to_cal+'dark_master.fit')
        except:
            prnt(self.filename,'Failed to open dark master %s' % self.path_to_cal+'dark_master.fit. Wrote to DR_errorlog.txt')
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Missing dark master at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'. Auto DR halted.\n')
            sys.exit()

        self.dark_h = self.dark_fits[0].header
        self.dark = self.dark_fits[0].data
        self.dxptime = self.dark_h['EXPTIME'] # store the exposure time for the dark master for scaling purposes
        exptime = light_h['EXPTIME'] # store light image exposure time


        try: # open filter-specific flat
            flat_fits = fits.open(self.path_to_cal+'flat_master_'+light_h['FILTER']+'.fit') 
            prnt(self.filename,'Successfully opened '+self.path_to_cal+'flat_master_'+light_h['FILTER']+'.fit')
        except:
            prnt(self.filename,'Failed to open flat master %s' % self.path_to_cal+'flat_master_'+light_h['FILTER']+'.fit. Wrote to DR_errorlog.txt')
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Missing '+light_h['FILTER']+'flat master at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'. Auto DR halted.\n')
            sys.exit()
        
        flat_h = flat_fits[0].header
        flat = flat_fits[0].data

        ## perform the actual data reduction
        if self.calibrate(light_h,light)==True:
            prnt(self.filename,'Calibrating image...' )

            bias_corrected_image = light - self.bias # subtract the bias
            dark_corrected_image = bias_corrected_image - (exptime/self.dxptime)*self.dark # scale the dark linearly w/ exptime and subtract
            final_image = dark_corrected_image / flat # divide by the flat field (already normalized)
            
            self.write_to_header(light_h)
            self.save_file(light_h, final_image, max(self.dates),self.filename)


        elif self.calibrate(light_h,light)=='OnlyDark': # auto dark
            prnt(self.filename,'Calibrating image...' )
            
            final_image = light / flat # divide by the flat field

            self.write_to_header(light_h)
            self.save_file(light_h, final_image, max(self.dates),self.filename)


        elif self.calibrate(light_h,light)=='Redundant':
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Attempted redundant calibration on '+self.filename+' at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'\n')
            prnt(self.filename,'Image already calibrated')
            self.save_file(light_h, light, max(self.dates),self.filename)
            

        elif self.calibrate(light_h,light)=='Binning':
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Image '+self.filename+' not 1x1 binning, rejected calibration at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'.')
            prnt(self.filename,'Image not 1x1 binning')

        elif self.calibrate(light_h,light)=='Temp':
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Image '+self.filename+' temp '+light_h['CCD-TEMP']+' degrees C, rejected calibration at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'.')
            prnt(self.filename,'Image taken at > -4 degrees C')

        elif self.calibrate(light_h,light)=='Size':
            with open('DR_errorlog.txt','a') as erlog:
                erlog.write('Image '+self.filename+' not full size, rejected calibration at '+strftime("%Y%m%d %H:%M GMT", gmtime())+'.')
            prnt(self.filename,'Image not full size')

        del flat_fits 

    def calibrated(self):
        return self.c

    def Source(self): # gathers source extraction data from .SRC file
        src = np.loadtxt(self.calibrated_path+self.filename.replace('.fits','.SRC'))
        objects = src[:,0:2] # pixel X,Y coordinates of the objects in question
        X_pos = src[:,0]
        Y_pos = src[:,1]
        A = src[:,8]
        B = src[:,9]
        theta = src[:,10]
        prnt(self.filename,'Gathered source data')
        return {'obj':objects,'X':X_pos,'Y':Y_pos,'A':A,'B':B,'theta':theta}

    def Convert(self): # converts obj list in pixel coordinate to RA-dec coordinates
        hdr = self.hdr
        w = wcs.WCS(hdr)
        objects = self.source['obj']
        world = w.wcs_pix2world(objects, 1)
        prnt(self.filename,'Converted coordinates to RA/Dec')
        return world

    def Photometry(self):
        ### perform aperture photometry
        hdr,img = self.hdr,self.img
        bkg = sep.Background(img)
        img_sub = img - bkg
        flux, fluxerr, flag = sep.sum_circle(img_sub, self.source['X'], self.source['Y'], self.aperture_size, err=bkg.globalrms, gain=1.0)
        magnitudes = -2.5*np.log(flux)
        median_i = np.median(magnitudes)


        ### retrieve magnitudes from catalog
        time = hdr['DATE-OBS']
        time = datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%f')
        filt = hdr['FILTER']
        objects = self.world


        v = Vizier(columns=['UCAC4','+_r','RAJ2000','DEJ2000','Bmag','Vmag','rmag'])
        output = {'id':[],'RA_C':[],'DEC_C':[],'RA_M':[],'DEC_M':[],'DIF':[],'MAG_R':[],'MAG_V':[],'MAG_B':[],'CMAG_R':[],'CMAG_V':[],'CMAG_B':[],'DATETIME':[],'IMGNAME':[]}
        cmags = []
        misfires = 0
        
        for n in range(len(objects)):
            catalog = 'UCAC4'
            result = v.query_region(coord.SkyCoord(ra=objects[n,0], dec=objects[n,1],
            unit=(u.degree, u.degree), frame='fk5'),radius='2s',catalog=catalog)
            try:
                result = result[0]
            except:
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

            ids = np.array(result['UCAC4'],str)
            ra = np.array(result['RAJ2000'],float)
            dec = np.array(result['DEJ2000'],float) # catalog RA and Dec
            dif = np.array(result['_r'],float)
            fluxtype = filt+'mag'
            if filt=='R':
                fluxtype = 'rmag'
            flux = np.array(result[fluxtype],float)

            for i in range(len(ids)):
                if dif[i] <= 2 and i==np.argmin(dif): # min residual value and less than 2 arcsec off
                    prnt(self.filename,'Star match in %s, mag %s, residual %s arcsec' % (catalog,flux[i],dif[i]))
                    output['id'].append(ids[i])
                    output['RA_C'].append(ra[i])
                    output['DEC_C'].append(dec[i])
                    output['RA_M'].append(objects[n,0])
                    output['DEC_M'].append(objects[n,1])
                    output['DIF'].append(dif[i])
                    output['DATETIME'].append(time)
                    output['IMGNAME'].append(self.filename)
                    cmags.append(flux[i])

                    
                else:
                    prnt(self.filename,'No star match within 2 arcseconds')
                    misfires += 1
                    continue


        prnt(self.filename,'Output %s stars' % len(output['id']))
        prnt(self.filename,'Output %s unique stars' % len(set(output['id'])))
        prnt(self.filename,'Missed %s objects' % misfires)
        prnt(self.filename,'Wrote magnitude data to sources.csv')
        print(' ')
        sleep(1)
        print("\033c")
        print('-'*termsize)

        cmags_nonan = [k for k in cmags if not math.isnan(float(k))]
        median_c = np.median(np.array(cmags_nonan))

        d = median_c - median_i
        magnitudes += d
        for n in range(len(objects)):
            for j in ['R','V','B']:
                magtype = 'MAG_'+j
                if j==filt:
                    output[magtype].append(magnitudes[n])
                else:
                    output[magtype].append('---')
                magtype = 'CMAG_'+j
                if j==filt:
                    output[magtype].append(cmags[n])
                else:
                    output[magtype].append('---')

        return output

    def Plot(self):
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
        name = self.filename.replace('Calibrated Images/','TestImages/')
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
   
    def Extract(self):
        self.source = self.Source()
        self.world = self.Convert()
        self.output = self.Photometry()
        self.writeData()
