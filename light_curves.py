# this is written to work in Python 2.7, as is the rest of the data pipeline, 
# the only major difference between this and Python 3.6 is the need to use the raw_input() function
# instead of just input(), as the two work differently in different versions

import pandas as pd
import numpy as np
import csv
import math
from datetime import datetime
from time import sleep
import matplotlib.pyplot as plt
import matplotlib
import sys
import os
import random


print("\033c")
rows, columns = os.popen('stty size', 'r').read().split()
termsize = int(columns)


def bold(msg):
    return u'\033[1m%s\033[0m' % msg

def header():
    print('-'*termsize)
    print(bold('Plotting Light Curves:').center(termsize))
    print('a part of The Guilford College Cline Observatory Data Analysis Pipeline Project'.center(termsize))
    print('-'*termsize)
    print('')

def printright(printout,clear=False,delay=False):
    if clear==True:
        print("\033c")
        header()
    print(' '*(termsize-len(printout)-8)+printout)
    if delay==True:
        sleep(1)

def options(head,items):
    print('')
    print('\t'+bold(head))
    print(('-'*(termsize-16)).center(termsize))
    for j in items:
        print('\t'+j)
    choice = raw_input('\tChoice: ')
    return choice

def print_slow(t,indent=1,speed=100):
    i = ' '*8*indent
    sys.stdout.write(i)
    #print(i,end='')
    for l in t:
        sys.stdout.write(l)
        sys.stdout.flush()
        sleep(random.random()*10.0/speed)
    print('')
    sleep(0.2)

def init():
    header()
    sleep(1)
    printright('[press ^C at any time to skip this intro]')
    try:
        print_slow('This program facilitates the following tasks:',indent=1)
        print_slow('Analysis of data collected and run through the data pipeline',indent=2)
        print_slow('Plotting of light curves for various optical filters and time spans',indent=2)
        print_slow('..........',indent=int((termsize/2-5)/8),speed=10)
        sleep(1)
        print("\033c")
        header()
        sleep(3)
    except KeyboardInterrupt:
        print("\033c")
        header()
        sleep(0.5)
    
def openData(filename):
    with open(filename) as csvFile:
        reader = csv.reader(csvFile)
        keys = next(reader)
    dictionary = dict()
    for i in keys:
        df = pd.read_csv(filename)
        dictionary[i]=np.array(df[i])
    return dictionary

init()
sources = openData('sources.csv')


while True:
    choice = options('Star Entry: ',
                ['[1] Star lookup by UCAC4 ID',
                '[2] Star lookup by coordinates',
                '[3] Bring up testing data',
                '[q] Quit the program'])

    if choice=='q':
        print("\033c")
        sys.exit()

    if choice=='3':
        choice='643-094180'
        printright('Choice Registered: Test Star '+choice,clear=True)
        indices = np.nonzero(sources['id']==choice)[0]
        printright('%s data points found' % len(indices),delay=True)

    if choice=='2':
        printright('Choice Registered: Coordinate Lookup',clear=True)
        coords = raw_input("\tRA and DEC coordinates in decimal decrees (format 'RA,DEC'): ")
        coords = coords.split(',')
        RA = float(coords[0])
        DEC = float(coords[1])
        difference = np.zeros(np.shape(sources['RA_M']))
        for j in range(len(sources['RA_M'])):
            RA_M = float(sources['RA_M'][j])
            DEC_M = float(sources['DEC_M'][j])
            difference[j] = np.sqrt((RA_M-RA)**2+(DEC_M-DEC)**2)

        indices = np.nonzero(difference<(2/60/60))[0]
        printright('%s data points found' % len(indices),delay=True)

    if choice=='1':
        choice = raw_input('\tUCAC4 ID: ')
        printright('Choice Registered: UCAC4 ID Star '+choice,clear=True)

        indices = np.nonzero(sources['id']==choice)[0]
        printright('%s data points found' % len(indices),delay=True)

    filt = options('Optical Filter Choice: ',
                ['[R] R band (red) fitler',
                '[V] V band (visual) filter',
                '[B] B band (blue) filter',
                '[a] View all filter data',
                '[q] Quit the program'])

    printright('Choice Registered: '+filt,clear=True,delay=True)
    if filt=='q':
        print("\033c")
        sys.exit()


    if not filt=='a':
        timeflag = options('Date/Time Entry: ',
                ['[1] Enter start/end manually',
                '[a] View all time data',
                '[q] Quit the program'])
                
        if timeflag=='q':
            print("\033c")
            sys.exit()

        if timeflag=='1':
            printright('Choice Registered: Manual Entry',clear=True,delay=True)
            print('')
            print('\t'+bold('Date/Time Entry: '))
            print(('-'*(termsize-16)).center(termsize))
            print('\tEnter dates as YYYY/MM/DD/HH/mm in GMT/24hr')
            start = raw_input("\tStart time: ")
            end = raw_input("\tEnd time: ")
            printright('Choice Registered: from %s to %s' % (start,end),clear=True)
            start = datetime.strptime(start,'%Y/%m/%d/%H/%M')
            end = datetime.strptime(end,'%Y/%m/%d/%H/%M')
        
        if timeflag=='a':
            start = datetime.strptime(np.amin(sources['DATETIME']), '%Y-%m-%d %H:%M:%S.%f')
            end = datetime.strptime(np.amax(sources['DATETIME']), '%Y-%m-%d %H:%M:%S.%f')
            printright('Choice Registered: All time data',clear=True,delay=True)

        mags,time = [],[]
        for x in indices:
            mags.append(sources['MAG_'+filt][x])
            time.append(datetime.strptime(sources['DATETIME'][x], '%Y-%m-%d %H:%M:%S.%f'))

        mags = [mags[x] for x in range(len(mags)) if isinstance(mags[x], float)]
        try:
            time = [time[x] for x in range(len(time)) if isinstance(mags[x], float)]
        except IndexError:
            print('\tNo %s magnitude data found' % filt)
            sleep(1.5)
            print('\tReturning to Star Lookup')
            sleep(2.5)
            print("\033c")
            sleep(0.6)
            header()
            sleep(0.6)
            continue

        saveflag = raw_input("\tSave plot as file? [y/n]: ")
        plt.figure(figsize=(10,8))
        plt.scatter(time, mags, c='k', marker='.', label=filt+' mag')
        plt.legend()
        plt.gca().invert_yaxis()

        plt.ticklabel_format(useOffset=False,axis='y')
        
        plt.xlim(start,end)
        plt.xlabel('Time')
        plt.ylabel('Magnitude')
        plt.xticks(rotation=50)
        plt.margins(0.2)
        plt.subplots_adjust(bottom=0.15)

        if choice=='C':
            plt.title('Star at %s, %s' % (RA,DEC))
            if saveflag=='y':
                filename = 'star_%s_%s.png' % (RA,DEC)
        else: 
            plt.title('UCAC4 %s' % choice)
            if saveflag=='y':
                filename = 'star_%s' % choice
        if saveflag=='y':
            plt.savefig(filename)
        plt.show()


    if filt=='a':
        timeflag = options('Date/Time Entry: ',
                ['[1] Enter start/end manually',
                '[a] View all time data',
                '[q] Quit the program'])

        if timeflag=='q':
            print("\033c")
            sys.exit()

        if timeflag=='1':
            printright('Choice Registered: Manual Entry',clear=True,delay=True)
            print('')
            print('\t'+bold('Date/Time Entry: '))
            print(('-'*(termsize-16)).center(termsize))
            print('\tEnter dates as YYYY/MM/DD/HH/mm in GMT/24hr')
            start = raw_input("\tStart time: ")
            end = raw_input("\tEnd time: ")
            printright('Choice Registered: from %s to %s' % (start,end),clear=True)
            start = datetime.strptime(start,'%Y/%m/%d/%H/%M')
            end = datetime.strptime(end,'%Y/%m/%d/%H/%M')
        if timeflag=='a':
            start = datetime.strptime(np.amin(sources['DATETIME']), '%Y-%m-%d %H:%M:%S.%f')
            end = datetime.strptime(np.amax(sources['DATETIME']), '%Y-%m-%d %H:%M:%S.%f')
            printright('Choice Registered: All time data',clear=True,delay=True)

        filt = ['R','V','B']
        mags,time = [[],[],[]],[[],[],[]]
        for x in indices:
            for j in range(len(filt)):
                if sources['MAG_'+filt[j]][x]=='---':
                    mags[j].append(sources['MAG_'+filt[j]][x])
                else:
                    mags[j].append(float(sources['MAG_'+filt[j]][x]))
                time[j].append(datetime.strptime(sources['DATETIME'][x], '%Y-%m-%d %H:%M:%S.%f'))

        nodata = 0
        for j in range(len(filt)):
            time[j] = [time[j][x] for x in range(len(time[j])) if isinstance(mags[j][x], float)]
            mags[j] = [mags[j][x] for x in range(len(mags[j])) if isinstance(mags[j][x], float)]
        

        saveflag = raw_input("\tSave plot as file? [y/n]: ")
        plt.figure(figsize=(10,8))
        plt.scatter(time[0], mags[0], c='r', marker='.', label='R mag')
        plt.scatter(time[1], mags[1], c='g', marker='.', label='V mag')
        plt.scatter(time[2], mags[2], c='b', marker='.', label='B mag')
        plt.legend()
        plt.gca().invert_yaxis()

        plt.ticklabel_format(useOffset=False,axis='y')
        
        plt.xlim(start,end)
        plt.xlabel('Time')
        plt.ylabel('Magnitude')
        plt.xticks(rotation=50)
        plt.margins(0.2)
        plt.subplots_adjust(bottom=0.15)

        if choice=='C':
            plt.title('Star at %s, %s' % (RA,DEC))
            if saveflag=='y':
                filename = 'star_%s_%s.png' % (RA,DEC)
        else: 
            plt.title('UCAC4 %s' % choice)
            if saveflag=='y':
                filename = 'star_%s' % choice
        if saveflag=='y':
            plt.savefig(filename)
        plt.show()

    if saveflag=='y':
        printright('Plot saved as %s' % filename,clear=True,delay=True)
    else:
        printright('Plot not saved',clear=True,delay=True)

    loopflag = options('What would you like to do next?',
                ['[1] Return to the beginning to plot another curve',
                '[q] Quit the program'])
    print("\033c")
    if loopflag=='q':
        sys.exit()

    if loopflag=='1':
        sleep(0.6)
        header()
        sleep(0.6)
        continue