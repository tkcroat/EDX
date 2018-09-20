# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 12:51:54 2016
Series of scripts to handle SEM-EDX files from NSS (either psmsa or emsa);  
guided by SEM_backfit_params, SEMEDXquantparams, 

@author: tkc
"""
import os, re, glob
import pandas as pd 
import numpy as np
from math import factorial # used by Savgol matrix
from io import StringIO
import datetime
from scipy import optimize
from PIL import Image, ImageDraw, ImageFont
import tkinter as tk  
import warnings
warnings.simplefilter(action = "ignore", category = RuntimeWarning)

def renamePSset(p_sname, str1, str2):
    ''' Find/replace rename for all 4 files .p_s, .psref, _pt1.psmsa and 
    change to .csi text removed bad space character problem or can be used 
    for any string replacement (i.e. misnamed file)
    '''
    newname=p_sname.replace(str1, str2)
    # Create new p_s file
    with open(newname,'w') as outfile:
        with open(p_sname,'r') as infile:
            mylines=list(infile)
            for i, line in enumerate(mylines):
                if 'psmsa' in line:
                    outfile.write(line.replace(str1,str2))
                elif 'psref' in line:
                    outfile.write(line.replace(str1,str2))
                else:
                    outfile.write(line)
    # delete original p_sfile
    os.remove(p_sname)

    # Rename all the *pt1.psmsa files
    pts=glob.glob("*%s*_pt*psmsa" % (p_sname.split('.')[0]))
    for i, fname in enumerate(pts):
        newname=fname.replace(str1,str2)
        os.rename(fname,newname)
    # Rename psref file
    ps=glob.glob("*%s*psref" % (p_sname.split('.')[0]))
    newname=ps[0].replace(str1, str2)
    os.rename(ps[0],newname)
    # rename any jpg files 
    jpgs=glob.glob("*%s*jpg" % (p_sname.split('.')[0]))
    for i, fname in enumerate(jpgs):        
        newname=fname.replace(str1, str2)
        os.rename(fname,newname)
    # now replace text within p_s file
    return 

def getelemmaps(thisrow, elements):
    '''Open selected element map images (NSS CSV style) for this NSS row (different available elements sets are selectable)
    returned as list of (elemname, maparray; similiar structure as for created Augermaps) '''
    availableelems=thisrow.iloc[0]['Elements']
    availableelems=[str(elem.strip()) for elem in availableelems.split(',')]
    if elements==['all']: # open and return all element maps
        elemlist=availableelems
    else: # passed value is normally a list of elements (not string with comma separated elements)
        elemlist=[str(el) for el in elements if el in availableelems]
    basename=thisrow.iloc[0]['Basename']
    if 'Grey' not in elemlist:
        elemlist.append('Grey') # also always grab secondary electron image
    elementmaps=[]
    for i, elem in enumerate(elemlist):
        fname=basename+'_'+elem+'.csv'
        try:
            myarr=np.loadtxt(fname, delimiter=',', skiprows=5) # direct load of image csv data into numpy array
            elementmaps.append([elem,myarr]) # package and return as elem name, element maps
        except:
            print("Couldn't open ", fname)
    return elementmaps

def getsinglemap(elementmaps,elem):
    '''Extract single numpy array map from set of elemental maps '''
    availelems=[]
    for i in range(0, len(elementmaps)):
        availelems.append(elementmaps[i][0])
    pos=availelems.index(elem)
    myarr=elementmaps[pos][1]
    return myarr
    
def batchNSSreader():
    '''Grab all params for all NSS x-ray image csv files (those produced from .si files by), pull details about each line, 
    renames element files (remove K, L, etc.)'''
    filelist=glob.glob('*Grey.csv')
    mycols=['Basename','SIname','Xres','Yres','Pixsize','Elements','Elemdetails']
    if os.path.isfile('NSScsvparams.csv'):
        NSScsvparams=pd.read_csv('NSScsvparams.csv', encoding='cp437')
    else:
        NSScsvparams=pd.DataFrame(columns=mycols) # Reconstructed each time?
    processed=np.ndarray.tolist(NSScsvparams.Basename.unique())
    for i, file in enumerate(filelist):
        base=file.split('_Grey.csv')[0]
        if base in processed: # skip csv file if already processed
            continue
        newcsvparams=pd.DataFrame(index=np.arange(0,1),columns=mycols) # single df row
        basestr='*'+base+'*'
        matchlist=glob.glob(basestr) # finds other element maps
        elemdetails=[] # new for each basefile
        elemlist=[] # comma separated element string
        for j, name in enumerate(matchlist): # loop through all elements 
            if 'Grey' in name:
                continue # not an element plot 
            tempstr=name.split(base)[1]
            newname=name.split(tempstr)[0]
            try:
                tempstr=tempstr.split('_')[1] # remove leading underscore
                tempstr=tempstr.split('.csv')[0] # remove file extension
            except: # Handles 100.csv image files in some pristine grains SI folders
                continue
            # Rename element map subfiles (but keep full info in elemdetails)
            elemdetails.append(tempstr)
            try:
                tempstr=tempstr.split(' ')[0] # has space only if no prior run
            except:
                pass
            elemlist.append(tempstr) # add bare element name
            # now rename by removing element detail
            if ' ' in name: # presence of space indicates not yet renamed (prior processing loop)
                newname=newname+'_'+tempstr+'.csv'
                os.rename(name,newname) # rename and remove subelem info (for ease of later finding image source file)
        newcsvparams.iloc[0]['Basename']=base
        newcsvparams.iloc[0]['Elements']=", ".join(elemlist)
        newcsvparams.iloc[0]['Elemdetails']=", ".join(elemdetails) # detail of elements (i.e. 'Mg K' instead of 'Mg')
        # find SIname, X and Y map resolutions and pixel size
        newcsvparams=readNSSparams(newcsvparams)# modify and return single-rowed df
        NSScsvparams=pd.concat([NSScsvparams, newcsvparams], ignore_index=True)
    NSScsvparams.to_csv('NSScsvparams.csv', index=False)
    return NSScsvparams

def readNSSparams(thisdfrow):
    '''Read NSS CSV extracted from SI header params from grey (others identical)
    passing single rowed df '''
    greyfile=thisdfrow.iloc[0]['Basename']+'_Grey.csv'
    if not os.path.isfile(greyfile):
        print('Data file ', greyfile, ' not found.')
        return thisdfrow
    with open(greyfile, 'r') as file:
        header=file.readlines()[0:5] # 5 header lines with filename, X, Y resolutions, pixel size, type
    # get SI name (can differ from CSV x-ray map image's saved name)
    tempstr=header[0].split('Image Name, ')[1]
    tempstr=tempstr.split(' Grey')[0]
    thisdfrow.iloc[0]['SIname']=tempstr
    # get X and Y resolutions
    tempstr=header[1].split('Pixels, ')[1]
    thisdfrow.iloc[0]['Xres']=int(tempstr)
    tempstr=header[2].split('Pixels, ')[1]
    thisdfrow.iloc[0]['Yres']=int(tempstr)
    tempstr=header[3].split('Size, ')[1]
    thisdfrow.iloc[0]['Pixsize']=tempstr # value and units (is this always in microns?)
    return thisdfrow
        
def csvbackup(filename):
    ''' version if df not loaded into memory... makes existing file into backup and checks backup
    list eliminating excessive ones (based on targetdate list)
    ''' 
    # if existing backup, rename
    now=datetime.datetime.now()
    mystr='*'+ filename +'*.csv'
    filelist=glob.glob(mystr)
    basename=filename+'.csv'
    if basename in filelist: #rename existing version 
        newname=filename+'_'+datetime.date.strftime(now, "%d%b%y")+'.csv'
        os.rename(basename, newname)
        filelist=glob.glob(mystr) # refresh filelist
    # targetdates gives ideal ages of list of backup files
    targetdates=[datetime.timedelta(120,0,0),datetime.timedelta(7,0,0)]
    
    dates=[] # list of file backup dates
    fileage=[] # age of backup
    for i,name in enumerate(filelist):        
        if '_' not in name: # unnecessary for this version 
            continue
        thisdate=name.split(filename+'_')[1]
        thisdate=thisdate.split('.csv')[0]
        thisdate=datetime.datetime.strptime(thisdate, "%d%b%y")
        age=now-thisdate
        fileage.append(age)
        dates.append([thisdate, age, name])
    dates.sort() # sort earliest to latest
    fileage.sort(reverse=True)
    # find list of files closest to ~4 mo old backup, 1 week old, and recent backup (2 days ish)
    #  keep at 4 mo and 1 week 
    keepindices=[]
    for i,thisage in enumerate(targetdates):
        # find closest entry to each
        if fileage: # skip if no such file names were found
            ind, age = min(enumerate(fileage), key=lambda x: abs(x[1]-thisage))
            keepindices.append(ind)
    if dates:
        keepindices.append(len(dates)-1) # always keep most recent backup
    
    # for list of lists, any way to just make list of element 1 
    for i, datelist in enumerate(dates):
        if i not in keepindices:
            os.remove(datelist[2])
    return
    
def replacelogentries(SEMfiles, Backfitparamslog, Integquantlog, Peakfitlog):
    '''After refitting misfits and checking new versions, replace those entries from prior saved logs with new ones '''
    # Read existing log files (with mostly good fits)
    try:
        Backlog=pd.read_csv('Backfitparamslog.csv', encoding='cp437')
        Integlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
        Peaklog=pd.read_csv('Peakfitlog.csv', encoding='cp437')
    except:
        print ('Prior background fit files were not saved!')
        return
    for index, row in SEMfiles.iterrows():
        thisname=SEMfiles.loc[index]['Filename']
        # remove existing entries
        Backlog=Backlog[Backlog['Filename']!=thisname]
        Integlog=Integlog[Integlog['Filename']!=thisname]
        Peaklog=Peaklog[Peaklog['Filename']!=thisname]
    # now add back new passed entries
    Backlog=pd.concat([Backlog,Backfitparamslog],ignore_index=True)
    Integlog=pd.concat([Integlog,Integquantlog],ignore_index=True)
    Peaklog=pd.concat([Peaklog,Peakfitlog],ignore_index=True)
    # Save modified files
    Backlog.to_csv('Backfitparamslog.csv', index=False)
    Integlog.to_csv('Integquantlog.csv', index=False)
    Peaklog.to_csv('Peakfitlog.csv', index=False)
    return Backfitparamslog, Integquantlog, Peakfitlog

def combinefiles(combinefilelist, savename):
    ''' Direct summation of counts from files in list; devise new filename and save directly '''
    EDXsum=openEDX(combinefilelist[0]) # open the first csv file as df
    # Assumed that these files have same energy structure 
    for i in range(1,len(combinefilelist)):
        EDXfile=openEDX(combinefilelist[i])
        mycols=['Counts', 'Backfit', 'Subdata', 'Gauss'] # all but derivative
        for i, colname in enumerate(mycols):
            EDXsum[colname]=EDXsum[colname]+EDXfile[colname]
    EDXsum=makesavgol(EDXsum) # recreate derivative
    EDXsum.to_csv(savename, index=False)
    # devise filename for saving (base)
    return

def makelogentry(df,nums, savename):
    '''constructs log entry for combined (summed) spectrum '''
    mycols=df.dtypes.index
    newrow=pd.DataFrame(index=np.arange(0,1),columns=mycols)
    newrow=df.iloc[0] # new values as series
    # change filenumber and filename
    newrow.Filenumber=nums
    newrow.Filename=savename
    # sum livetime, realtime
    newrow.Livetime=df.Livetime.sum()
    newrow.Realtime=df.Realtime.sum()
    # average detected, converted, stored, deadfraction
    newrow.Detected=df.Detected.mean()
    newrow.Converted=df.Converted.mean()
    newrow.Stored=df.Stored.mean()
    newrow.Deadfraction=df.Deadfraction.mean()
    return newrow # pass back the new entry as a series 
    
def combineEDX(df):
    ''' Direct sum of multiple spectra from same basename, point # (simplest comination logic without
    using logbook lastnumber style)
    also make new log entry  csv name is 1and2'''
    myarr=df.Basename.unique()
    for i, val in enumerate(myarr): # loop through basenames
        thisbase=df[df['Basename']==val]
        # usually only 1 point but occasionally more (only combine if same defined point)
        mypts=thisbase.Point.unique() # list of pts in
        mypts=np.ndarray.tolist(mypts) # convert to iterable list
        for j, pt in enumerate(mypts):  # then loop through points
            thissample=thisbase[thisbase['Point']==pt]
            if len(thissample>1): # combine via straight summation
                combinefilelist=thissample.Filename.unique()
                combinefilelist=np.ndarray.tolist(combinefilelist)
                # devise filename for combined file (basename + spnums +pt)
                filenumbers=thissample.Filenumber.unique()
                filenumbers=np.ndarray.tolist(filenumbers)
                nums="+".join(str(x) for x in filenumbers)
                savename=val+'('+nums+')_pt'+ str(pt)+'.csv'
                combinefiles(combinefilelist, savename) # sums spectra and saves directly
                # append entry to log file
                thisentry=makelogentry(thissample, nums, savename)
                # append series to log df 
                df=df.append(thisentry)
    csvbackup('EDXparamlog')
    df.to_csv('EDXparamlog.csv', index=False)
    return df # return log with all the new entries
    
def makepsimage(filename,SpatialAreas):
    '''Open psref as tif, overlay all rectangles and points from SpatialArea and directly save as jpg '''
    thisimage=Image.open(filename)
    draw=ImageDraw.Draw(thisimage) # single draw instance to label above image
    ttfont=ImageFont.truetype('arial.ttf', size=20)
    # TODO make sure font size adjustment is working
    # NSs box/area positions are slightly up and left of those direct from p_s file (use adjustment below)
    for index, row in SpatialAreas.iterrows():
        x1=SpatialAreas.iloc[index]['X1']
        y1=SpatialAreas.iloc[index]['Y1']
        x2=SpatialAreas.iloc[index]['X2']
        y2=SpatialAreas.iloc[index]['Y2']
        adj1=SpatialAreas.iloc[index]['Adj1'] # x scaling factor
        adj2=SpatialAreas.iloc[index]['Adj2'] # y scaling factor 
        if x2>0: # rectangular area region
            x1=x1*adj1/512
            x2=x2*adj1/512
            y1=y1*adj2/384
            y2=y2*adj2/384
            draw.rectangle((x1,y1,x2,y2), outline='red') # red rectangle at specified position
            message=str(index+1)
            draw.text((x2+2,y2+2),message, font=ttfont)
        else: # point 
            x1=x1*adj1/512
            y1=y1*adj2/384
            draw.rectangle((x1-4,y1-4,x1+4,y1+4), fill='red')
            message=str(index+1)
            draw.text((x1+2,y1+2),message,font=ttfont)
    jpgname=filename.split('.psref')[0]+'.jpg'
    thisimage.convert('RGB').save(jpgname)
    thisimage.close()
    return
	
def getpsinfo(psname):
    '''pass p_s file string and extract/return points info
    still unclear why image dimensions are listed at 96% of psref (i.e 491x368 instead of 512x384)	'''
    with open(psname, 'r') as file:
        filedata = file.read()
    pattern=re.compile(r'\.psmsa') # needed params are between each psmsa and next new line (\. is literal)
    matches=re.finditer(pattern,filedata)
    ptbreaks=[]
    # after parsed into these breaks contains psmsa name on first line (redundant) and has x,y, etc. params on second
    for m in matches:
        ptbreaks.append(m.start())
    dim=len(ptbreaks)
	# now make the dataframe of appropriate length
    mycols=['Filenumber','Filename','Areanumber','X1','Y1','X2','Y2','Adj1','Adj2','Areafract']
    SpatialAreas=pd.DataFrame(index=np.arange(0,dim), columns=mycols) # easier to use frame than dict
    ptbreaks.append(len(filedata)) # add EOF position
    ptitems=[] # list of strings holding data for each point parsed using .psmsa identified string
    for i in range(0,len(ptbreaks)-1):
        ptitems.append(filedata[ptbreaks[i]:ptbreaks[i+1]])
	# now extract needed numbers from each ptitem (truncate at /n and then just comma separated floats)
    for i, val in enumerate(ptitems): # extract and store info in df (which has same dimension)
        newval=val.split('\n')[1] # split data string at newline
        newval=newval.split('\n')[0]
        floatlist=[float(s) for s in newval.split(',')] # turn comma separated string into list of floats
        intlist=[int(float) for float in floatlist] # convert float list to int list
        SpatialAreas=SpatialAreas.set_value(i,'Areanumber',i+1) # pts are not zero-based
        SpatialAreas=SpatialAreas.set_value(i,'X1',intlist[0]) # pixel positions  in field
        SpatialAreas=SpatialAreas.set_value(i,'Y1',intlist[1])
        SpatialAreas=SpatialAreas.set_value(i,'Adj1',intlist[2])
        SpatialAreas=SpatialAreas.set_value(i,'Adj2',intlist[3])
        if len(intlist)==8:# indicates area spectrum
            SpatialAreas=SpatialAreas.set_value(i,'X2',intlist[4])
            SpatialAreas=SpatialAreas.set_value(i,'Y2',intlist[5])
            # calculate fractional area (relative to 512 x 384 full P&S image size)
            thisarea=(intlist[4]-intlist[0])*(intlist[4]-intlist[0])/(512*384)
            SpatialAreas=SpatialAreas.set_value(i,'Areafract',thisarea)
    return SpatialAreas

def processpointshoot():
    '''Finds psref and associated p_s file to get positions of points/rectangles  within psref file '''
    # for each psref in list, find associated p_s
    filelist=glob.glob('*.psref')
    # parse into chunks using SPOT or AREA
    mycols=['Filenumber','Filename','Areanumber','X1','Y1','X2','Y2','Adj1','Adj2','Areafract']
    SpatialAreasLog=pd.DataFrame(columns=mycols) # empty df log
    for i, filename in enumerate(filelist): # loop through all psref (tif) files	
        # find associated .p_s file
        psname=filename.split('.psref')[0]+'.p_s'
        # if p_s file exists, call getpsinfo to open and process
        if os.path.isfile(psname):
            SpatialAreas=getpsinfo(psname)
        else:
            print("Can't find p_s file associated with ", filename)
            continue # move to next psref and skip
        # Add filenumber and filename to each area number
        pattern=re.compile(r'(\(\d+\))') # find spectral number inside parens

        match=re.search(pattern, filename)
        if not match:
            print('Skipped: No spectral number in parenthesis in', filename)
            continue
        # basename=filename[0:match.start()] # pulls out base name
        specnum=int(filename[match.start()+1:match.end()-1]) # gets spectral number within parentheses
        for index,row in SpatialAreas.iterrows(): # copy common name to each row
            SpatialAreas=SpatialAreas.set_value(index,'Filenumber',specnum)
            SpatialAreas=SpatialAreas.set_value(index,'Filename',filename)
        # convert psref and superimpose spatial areas (directly saved as jpg with same filename)
        makepsimage(filename,SpatialAreas) # create jpg with P&S overlaid areas (directly saved)
        # write spatial areas from each to log
        SpatialAreasLog=pd.concat([SpatialAreasLog,SpatialAreas], ignore_index=True)
    SpatialAreasLog=SpatialAreasLog[mycols]
    csvbackup('Spatialareaslog')
    SpatialAreasLog.to_csv('Spatialareaslog.csv',index=False)
    return SpatialAreasLog
    
def returnoverlaps(fitranges):
    ''' Returns the energy ranges (as list of lists) where background fits overlap (for cross-fading
    of adjacent backgrounds'''
    overlaps=[]
    for i in range(1, len(fitranges)):
        range1=fitranges[i-1]        
        range2=fitranges[i]
        end1=int(range1.split('-')[1])
        start2=int(range2.split('-')[0])
        overlaps.append([start2,end1]) # append all overlaps even if ranges don't overlap
    return overlaps

def rangefromstring(x):
    result = []
    for part in x.split(','):
        if '-' in part:
            a, b = part.split('-')
            a, b = int(a), int(b)
            result.extend(range(a, b + 1))
        else:
            a = int(part)
            result.append(a)
    return result

def findminpoints(rangestrings, EDXfile):
    '''Pass index energy range (or ranges) and find/return minimum(s) in EDX file within these ranges 
    rangestring is index # range as string '''
    stringlists=[str(s) for s in rangestrings.split(',')]
    addlist=[]
    for i, val in enumerate(stringlists):
        thisminrange=rangefromstring(val)
        dfslice=EDXfile[min(thisminrange):max(thisminrange)] # slice based on this index range
        thismin=dfslice['Counts'].idxmin()
        addlist.append(int(thismin))
    return addlist # list of index #s of minimums within the ranges specified 
            
def unpackfitregs(df, EDXfile):
    ''' Loaded data frame has ev range, list of background regions and fit type
    unpack from dataframe into list of each 
    1) fitrange is total ev range (i.e.0-100), 2 )fitpts are index #s (or energy in eV) of regions commonly without peaks
    3) mandminpts - mandatory included index # of minimum vals between certain peaks (forced inclusion into background fit region)
    4) fittype (mostly parabola) and 5) threshold for derivative knockout (but threshold not applied to mandatory pts)
    add overlap range between adjacent fits''' 
    Fitregions=[]
    # TODO test to ensure that normal eV range corresponds to range of indices
    for i in range(0,len(df)):
        tempstr=df.iloc[i]['Backgroundregs']
        indexrange=rangefromstring(tempstr) # converts string describing range to actual range
        # This forces inclusion of lowest point in a given range as defined in minpoint cols
        # TODO set this up so that multiple strings can be used
        mandminlist=[] # empty default list of index # of mandatory included points for backfit
        if str(df.iloc[i]['Minpoint'])!='nan':
            mandminlist=findminpoints(df.iloc[i]['Minpoint'], EDXfile)
        Fitregions.append([df.iloc[i]['Fitrange'],indexrange, mandminlist, df.iloc[i]['Fittype'], df.iloc[i]['Threshold']])
    return Fitregions

def makesavgol(df):
    '''Perform python smooth-diff used to guide selection of background regions for SEM-EDX or TEM-EDX spectra
    ''' 
    df['Savgol']=0.0  # add/initialize col for 2nd deriv Sav-gol
    thisdim=len(df) 
    thisreg=df['Counts'] # convert to Series (keep these index)
    myarr=np.asarray(thisreg) # convert to numpy array
    window_size=11
    deriv=2 
    order=2 # order of savgol fit 
    rate=1
    order_range = range(order+1) # range object
    half_window = (window_size -1) // 2 # type int
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    # b is matrix 3 by window size
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv) # series as long as array
    # linalg.pinv gets pseudo-inverse of a matrix (window-sized series)
    # .A of any matrix returns it as ndarray object 
    
    # Pad the signal at the extremes with values taken from the signal itself
    firstvals = myarr[0] - np.abs(myarr[1:half_window+1][::-1] - myarr[0] )
    lastvals = myarr[-1] + np.abs(myarr[-half_window-1:-1][::-1] - myarr[-1])
    myarr= np.concatenate((firstvals, myarr, lastvals))
    # Now convolve input signal and sav-gol processing 1D array .. thisreg is numpy array w/ savgol results
    myarr=np.convolve( myarr, m[::-1], mode='valid')
    
    thisreg=pd.Series(myarr) # convert array to series 
    thisreg.loc[0:thisdim]=myarr # copies numpy array but keeps same indices
    df['Savgol']=thisreg # copy deriv column to dataframe
    return df #  returns savitsky-golay smooth diff over same full region 

# TESTING
# EDXfile=openEDX('GC6bands(1).csv')

def openEDX(EDXfileName):
    '''Open csv as dataframe if it exists or if not strip header from psmsa/emsa and import as dataframe ''' 
    csvname=str(EDXfileName.split('.')[0])+'.csv'
    try:
        EDXfile=pd.read_csv(csvname, encoding='cp437')
    except: # if csv doesn't exist, just open/strip psmsa
        with open(EDXfileName, 'r') as file:
            filedata = file.read()
        filedata =filedata.split('#SPECTRUM    :')[1]
        filedata =filedata.split('#ENDOFDATA   : ')[0]
        thisdata=StringIO(filedata)        
        EDXfile=pd.read_csv(thisdata)
        try:        
            EDXfile=EDXfile.drop(EDXfile.columns[[2]], axis=1) # drop erroneous 3rd column if present
        except:   
            print('') # ignore error if 3rd column not present
        EDXfile.columns=['Energy','Counts']    
    return EDXfile # should return data as pandas dataframe

def openorcreatelogbook(filelist):
    ''' Looks for existing csv or xls log file, load and pass in dict if present... if not found makes new one by calling makeblanklog ''' 
    logfile=glob.glob('*EDX_log*') # find ... Auger_logbook.
    if len(logfile)>0: # found logbook
        if len(logfile)>1:
            print('Error: multiple EDX logs found in folder... using', logfile[0])
        name=logfile[0]
        if '.xls' in name: # open log tab of existing excel file
            EDXlogbook=pd.read_excel(name, sheetname='Log')        
        if '.csv' in name: # open csv
            EDXlogbook=pd.read_csv(name)
        logbool=True
    else:
        print('No EDX logbook found... add sample names directly to paramlog')
        EDXlogbook=makeblanklog(filelist)        
        logbool=False
    return EDXlogbook, logbool

# If log not created during session, just use params file instead of sepsarate xls
def makeblanklog(filelist):
    ''' Make blank Excel log matching existing set of files (in cases where one was not created during data collection'''
    mycols=['Project', 'Basename', 'Filenumber', 'Lastnumber', 'Point', 'Filename', 'FilePath', 'Sample', 'Comments']    
    SEMlogbook=pd.DataFrame(columns=mycols) # blank df 
    # get project name from directory 
    fullpath=os.getcwd()
    pattern=re.compile(r'(\\)')
    match=re.finditer(pattern, fullpath)
    indices=[m.start(0) for m in match]
    projname=fullpath[indices[-1]+1:] # get project name (last folder) from full path
    for i, filename in enumerate(filelist):
        Samplelogrow=pd.DataFrame(index=np.arange(0,1), columns=mycols) # single df row for given file 
        pattern=re.compile(r'(\(\d+\))') # find spectral number inside parens
        match=re.search(pattern, filename)
        if match:
            basename=filename[0:match.start()] # pulls out base name
            specnum=int(filename[match.start()+1:match.end()-1]) # gets spectral number within parentheses
        else: # must have different naming convention (not auto-named)
            basename=filename.split('.')[0]
            specnum=1
        if '.psmsa' in filename: # emsa doesn't have point #
            ptnum=int(filename.split('_pt')[1].split('.')[0]) # gets pt number from 'name(10)_pt1.psmsa'
            Samplelogrow=Samplelogrow.set_value(0,'Point', ptnum)
        elif '.emsa' in filename: # no pt number for emsa ... set to 1
            ptnum=1
            Samplelogrow=Samplelogrow.set_value(0,'Point', ptnum)
        Samplelogrow=Samplelogrow.set_value(0,'Basename', basename)
        Samplelogrow=Samplelogrow.set_value(0,'Filenumber', specnum)

        Samplelogrow=Samplelogrow.set_value(0,'Filename', filename)
        Samplelogrow=Samplelogrow.set_value(0,'Project', projname)
        Samplelogrow=Samplelogrow.set_value(0,'FilePath', fullpath)
        SEMlogbook=pd.concat([SEMlogbook,Samplelogrow], ignore_index=True)
    SEMlogbook.sort_values(['Basename', 'Filenumber'])
    SEMlogbook=SEMlogbook[mycols] # reorder columns to standard
    print('Blank logbook created for project ', projname, '; Sample names and comments can be manually entered .')
    return SEMlogbook

def fitparabola(df, EDXfileName, fitrange):
    '''Pass appropriate chunk from Auger spectral dataframe, perform polynomial/parabola fit
    return chunk with backfit column added
    fitrange for error handling only'''
    xcol=df['Energy']
    ycol=df['Counts'] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
        A,B,C=np.polyfit(xcol, ycol, 2)
    except: # deal with common problems with linregress
        print('Fitting error from ', "{0:.2f}".format(df.Energy.min()),'to ',"{0:.2f}".format(df.Energy.max()), ' in file ', EDXfileName)
        print('df is empty =', df.empty) # indicate problem with empty frame probably due to thresholds
        print('fitrange =', fitrange) 
        fitparams=('n/a','n/a','n/a') # return all n/a
        return df, fitparams
    fitparams=(A, B, C) # tuple to return coeffs of 2nd order poly fit
    for index,row in df.iterrows(): # write this fit into this chunk of data (redundant?)
        xval=df.loc[index]['Energy']
        yval= A * xval**2+ B * xval + C
        df=df.set_value(index, 'Backfit', yval)
    return df, fitparams

def fitcubic(df, EDXfileName, fitrange):
    '''Pass appropriate chunk from EDX , perform polynomial/cubic fit
    return chunk with backfit column added
    fitrange for error handling only'''
    xcol=df['Energy']
    ycol=df['Counts'] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
        A,B,C,D=np.polyfit(xcol, ycol, 3)
    except: # deal with common problems with linregress
        print('Fitting error from ', "{0:.2f}".format(df.Energy.min()),'to ',"{0:.2f}".format(df.Energy.max()), ' in file ', EDXfileName)
        print('df is empty =', df.empty) # indicate problem with empty frame probably due to thresholds
        print('fitrange =', fitrange) 
        fitparams=('n/a','n/a','n/a') # return all n/a
        return df, fitparams
    fitparams=(A, B, C, D) # tuple to return coeffs of 2nd order poly fit
    for index,row in df.iterrows(): # write this fit into this chunk of data (redundant?)
        xval=df.loc[index]['Energy']
        yval= A * xval**3+ B * xval**2 + C*xval + D
        df=df.set_value(index, 'Backfit', yval)
    return df, fitparams
    
def findfitregion(df, fitregion, mandminpts, threshold, fitrange, EDXfileName):
    '''df is normally EDXfile Passing single list of allowable index #s for background fits (no duplicates) 
    remove those with high from list of allowable indices any that show high smoothed-derivatives (i.e. not good for background fitting 
    fitrange and EDXfilename -- error handling only'''
    Backfitdf=df.loc[[x for x in fitregion]] 
    # these are loaded from EDX_backfit_regions.csv
    Backfitdf=Backfitdf.dropna(subset=['Counts']) # drops above (set to na by ix)
    # now additionally filter out those with derivative above threshold value
    Backfitslice=Backfitdf[(Backfitdf['Savgol']<threshold) & (Backfitdf['Savgol']>-threshold)]
    if Backfitslice.empty==True:
        print('Threshold too low for ', fitrange, ' in ', EDXfileName)
        try:        
            while len(Backfitslice)<4: # reslice until at least 3 points appear for fitting this region
                threshold=threshold+1 # incrementally raise threshold
                Backfitslice=Backfitdf[(Backfitdf['Savgol']<threshold) & (Backfitdf['Savgol']>-threshold)]
            print ('Threshold reduced to ', str(threshold))
        except KeyboardInterrupt: # probably not necessary
            print('Threshold while loop interrupted!')
    # Add in the mandatory minimum fitting points (no threshold applied).. .list of ints
    for i, val in enumerate(mandminpts):
        # test if this index # is included in Backfitslice
        if val not in Backfitslice.index:
            Backfitslice=Backfitslice.append(df.loc[val]) # adds in mandatory backfit values from EDXfile df
    # pass back this set of backfitpts (stored in backfitlog for each region)
    backptslist=list(Backfitslice.index.values)
    return Backfitslice, backptslist

def findelemregions(Elements, EDXquantparams):
    ''' Takes element string and returns standard Elemdata for each elem symbol containing params 
    needed for peak finding and quant 
    tuple for integ peak is symbol, ideal peak index #, and integ kfactor
    don't apply energy shifts here... apply later when doing integrate''' 
    Elemdata=[]
    for i, elem in enumerate(Elements):
        try:
            # find row in AESquantparams for this element
            thiselemdata=EDXquantparams[(EDXquantparams['element']==elem)]
            thiselemdata=thiselemdata.squeeze() # series with this elements params
            
            # integ peak position value is relative to negpeak in smooth-diff (i.e. -5 is 5 eV below ideal negpeak)
            idealindex=int((thiselemdata.energy+.01)*100) # ideal index value of EDX-EDX peak from energy in keV
            kfact=thiselemdata.kfactor # typical sensitivity k-factor associated with element for integration
            errkfact=thiselemdata.errkfact 
            mass=thiselemdata.mass
            maxshift=int(thiselemdata.maxshift) # on indices so must be int
            # full peak width in keV from EDXquantparams (usually 0.15keV or 15 channels at 0.1eV/chan)
            # integration width in channels for direct integration for this element
            width=int(((thiselemdata.fullwidth*100)-1)/2) 
            # total # of channels in AESquantparams but include n-1/2 channels on either side of peak center (usually width is 8 channels)
            #Elemdata is a list (of length number of elements) containing length5 tuples
            elemtuple=(elem, idealindex, maxshift, width, kfact, errkfact, mass) # add tuple with info for this element
            Elemdata.append(elemtuple) # now contains proper limits on fitting regions 
        except:
            print('Quant parameters not properly loaded for', elem)
    return Elemdata
    
def fitgauss(df, halfwidth, elem, EDXfileName, savegauss=True):
    ''' Gaussian fit of direct peaks (pass EDXfile just around peaks region
    no need to save Gaussian fit, just return width and other params 
    integwidth pass from AESquantparams value'''
    # Remove any nan values from peak region (shouldn't be any though)
    df=df.dropna(subset=['Subdata']) # remove nan entries from peak
    # Estimate initial Gaussian parameters from data
    xc=df['Subdata'].idxmax() # estimate center based on peak max index
    xc=df.loc[xc]['Energy'] # associated energy value near center
    peakarea=df['Subdata'].sum()  # decent area estimate
    y0=0 #
    width=0.01*(2*halfwidth+1) # full width estimate in keV from half-width in channels
    params0=[xc,width,peakarea,y0] # initial params list (first guess at gaussian params)
    
    xcol=df['Energy']
    ycol=df['Subdata']
    xcol=xcol.as_matrix() # convert both to numpy matrices
    ycol=ycol.as_matrix()
    
    # define standard gaussian funct (xc, width, area and yoffset are init params)
    gaussian=lambda params, x: params[3]+params[2]/(params[1]*np.sqrt(2*np.pi))*np.exp(-((x-params[0])**2/(2*params[1]**2)))
    
    # thisgauss= gaussian(params0,xcol) 
    errfunc=lambda p, xcol, ycol: ycol- gaussian(p,xcol) # lambda error funct definition
    # sigma2FWHM = lambda sigma: sigma * sqrt(2 * log(2)) * 2 / sqrt(2) # convert Gaussian widths to FWHM?
    
    try:
        fitparams, cov, infodict, mesg, ier =optimize.leastsq(errfunc,params0,args=(xcol,ycol),full_output=True)    
        ss_err=(infodict['fvec']**2).sum()
        ss_tot=((ycol-ycol.mean())**2).sum()
        rsquared=1-(ss_err/ss_tot)
        
    except: # fitting problem 
        print('Gaussian fitting error for', elem, ' peak in file ', EDXfileName)
        fitparams=('n/a','n/a','n/a','n/a') # return all n/a
        rsquared='n/a'
        ier='n/a'
        return df, fitparams, rsquared, ier
    if savegauss==True:
        df['Gauss']='' # add col for gaussian fit
        for index,row in df.iterrows():
            xval=df.loc[index]['Energy']
            yval=fitparams[3]+fitparams[2]/(fitparams[1]*np.sqrt(2*np.pi))*np.exp(-((xval-fitparams[0])**2/(2*fitparams[1]**2)))
            df.set_value(index,'Gauss',yval)
    return df, fitparams, rsquared, ier

def fitpeaks(EDXfile, Elemdata, logmatch, savegauss, peakcols, integcols):
    ''' Gaussian fit of major peaks in single spectrum, shift is list of energy shifts of negpeak (same order as Eledata (opens source spectrum as EDXfile, 
    fits peak backgrounds above and below using Elemdata, also saves linear fit params to logdataframe with position/amplitude/etc;
    desired elements out of data range are skipped (in prior findindices function)
    # Saving of gaussian fits of peaks could be stored as separate csv if this was ever desired... probably not
    ''' 
    EDXfileName=logmatch.Filename # only used for error reporting 
    # Create temp df to hold and pass linear fit data    
    Peakfits=pd.DataFrame(columns=peakcols) # blank df for this spectrum's peak fits
    Integresults=pd.DataFrame(columns=integcols) # blank df for this spectrum's integration results
    
    # fit all elemental peaks with gaussian, determine shift and perform integration (incl. error)
    for i, (elem, idealindex, maxshift, halfwidth, kfact, errkfact, mass) in enumerate(Elemdata):
        Peakfitrow=pd.DataFrame(index=np.arange(0,1),columns=peakcols) # single dataframe row for this
        Integresultrow=pd.DataFrame(index=np.arange(0,1),columns=integcols) # blank df row
        # linear fit below this elem's peak (shifts and adjustments already made)
        # use 10 more channels than those used for integration for gaussian fits
        fitregion=EDXfile[idealindex-halfwidth-5:idealindex+halfwidth+6]
        if fitregion.empty==True: # skip if no data present (peak out of range problem)
            continue
        # Gaussian fit of subtracted data peaks > 50 cnts
        if fitregion['Subdata'].max()>50: # add flag to skip gaussian
            fitregion, fitparams, rsquared, ier = fitgauss(fitregion, halfwidth, elem, EDXfileName, savegauss=True)
            if savegauss==True: # save Gaussian peaks as separate column
                if 'Gauss' not in EDXfile.dtypes.index: # add col if not already present                
                    EDXfile['Gauss']='' # add blank col for gaussian fit if not present
                # copy gaussian fit to Augerfile... fitregion only modified in new Gauss peak fit column
                EDXfile.loc[fitregion.index,fitregion.columns]=fitregion
                      # determination of peak shift 
        # If gaussian fit is successful set center integration channel to index nearest xc
        # ier flag of 1,2,3,4 if fit succeeds but rsquared threshold is better
            if rsquared!='n/a': # somewhat successful gaussian fit 
                if rsquared>0.4:
                    xc=fitparams[0] # center of gaussian fit in keV
                    centerindex=int((xc+.01)*100)
                    shift= centerindex- idealindex # energy shift in channels
                    if abs(shift)>maxshift: # maxshift is element specific maximum move of integration window
                        # common problem with weak peaks... only use print for troubleshoot
                        # print('Warning: Gaussian shift of ', str(shift), ' channels indicated for ', elem, ' in ', EDXfileName)
                        if shift>0: # keep peak shift the same but only allow 3 channel shift in integration window
                            centerindex=idealindex+maxshift # set to max shift
                        else:
                            centerindex=idealindex-maxshift
            # TODO Maybe a better way of setting maximum allowable shift        
                else: 
                    # common problem with mass fitting so skip print report
                    # print('Low quality gaussian fit for ', elem, ' in ', EDXfileName)
                    centerindex=idealindex # already stores index number of central peak (ideal - sm-diff shift value)
                    shift='n/a'
            # Write gaussian fit params to peakfit (eventually copied to peakfitlog)
                    
            else: # Fit attempted but failed result
                print ('Fit attempted but result failed for ', elem, ' in ', EDXfileName)
                fitparams=['n/a','n/a','n/a','n/a']            
                rsquared='n/a'
                
        else: # indication of failed Gaussian fit (use prior knowledge of peak position)
            # common problem with weak peaks... only use print for troubleshoot
            # print('Skip gaussian fit of tiny ', elem, ' peak in ', EDXfileName)
            # set center integration channel to value passed by integpeak 
            # this is ideal energy value but adjusted by shift found using smooth-diff quant method
            centerindex=idealindex # already stores index number of central peak (ideal - sm-diff shift value)
            shift='n/a'
            fitparams=['n/a','n/a','n/a','n/a']            
            rsquared='n/a'
        # Perform integration over peak center channel + integwidth on either side 
        EDXpeak=EDXfile[centerindex-halfwidth:centerindex+halfwidth+1]
        integcounts=EDXpeak['Subdata'].sum() # get counts sum 
        backgroundcnts=EDXpeak['Backfit'].sum() # sum counts over identical width in background fit
        # Used for peak significance i.e. typically 2 sigma of background integration over identical width
        # full integ width is 1.2*FWHM but integwidth here is closest integer half-width

        # end of element loop
        Peakfitrow.loc[0]['Element']=elem
        Peakfitrow.loc[0]['Xc']=fitparams[0]
        Peakfitrow.loc[0]['Width']=fitparams[1]
        Peakfitrow.loc[0]['Peakarea']=fitparams[2]
        Peakfitrow.loc[0]['Y0']=fitparams[3]
        Peakfitrow.loc[0]['Rsquared']=rsquared        
        Peakfits=pd.concat([Peakfits, Peakfitrow], ignore_index=True) # copy peak rows individually to df
  
        # Copy integration results for this peak into df row
        Integresultrow.iloc[0]['Element']=elem
        Integresultrow.iloc[0]['Energy']=centerindex # index of center as determined by fitting (if successful)
        Integresultrow.iloc[0]['Shift']=shift # energy shift from ideal in channels (0.01 eV)
        Integresultrow.iloc[0]['Rawcounts']=EDXpeak['Counts'].sum() 
        Integresultrow.iloc[0]['Backcounts']=backgroundcnts
        Integresultrow.iloc[0]['Subtractedcounts']=integcounts
        # Adjusted counts must be determined later for pathological overlaps
        # 2 sigma err due to counting statistics
        Integresultrow.iloc[0]['% err']=round(2/np.sqrt(integcounts),3)
        Integresultrow.iloc[0]['Significance']=round(integcounts/(np.sqrt(backgroundcnts)),3)
		# TODO add 2/sqrt(n) calc of associated percent error (also can calculate later)
        Integresultrow.iloc[0]['Correctedcounts']=integcounts*kfact/mass
        # Calculated combined error for 2sig counting stats + loaded k-factor error
        comberr=np.sqrt(errkfact**2+(2/np.sqrt(integcounts))**2)
        # calculate error in Correctedcounts for given elemental peak
        Integresultrow.iloc[0]['Errcorrcnts']=(integcounts*kfact/mass)*comberr
        Integresultrow.iloc[0]['Kfact']=kfact
        Integresultrow.iloc[0]['Fullwidth']=2*halfwidth
        Integresults=pd.concat([Integresults,Integresultrow], ignore_index=True)
        
    # assign params that are common to this spectrum (all elemental peaks)
    for index,row in Peakfits.iterrows(): 
        Peakfits.loc[index]['Filenumber']=logmatch.Filenumber   
        Peakfits.loc[index]['Basename']=logmatch.Basename
        Peakfits.loc[index]['Filename']=logmatch.Filename
        Peakfits.loc[index]['Point']=logmatch.Point
        Peakfits.loc[index]['Filepath']=logmatch.FilePath
        Peakfits.loc[index]['Sample']=logmatch.Sample
        Peakfits.loc[index]['Comments']=logmatch.Comments
    for index,row in Integresults.iterrows(): # assign
        Integresults.loc[index]['Filenumber']=logmatch.Filenumber   
        Integresults.loc[index]['Filename']=logmatch.Filename
        Integresults.loc[index]['Basename']=logmatch.Basename
        Integresults.loc[index]['Point']=logmatch.Point
        Integresults.loc[index]['Filepath']=logmatch.FilePath
        Integresults.loc[index]['Sample']=logmatch.Sample
        Integresults.loc[index]['Comments']=logmatch.Comments
    Peakfits=Peakfits[peakcols] # put back in original order
    Integresults=Integresults[integcols] # put back in original order
    return EDXfile, Peakfits, Integresults # df with direct peak fitting info for all areas/ all elements

def crossfadeback(EDXfile, Backfitparams, overlaps):
    ''' For overlapping regions (defined as index #s in EDX_backfit_regions), blend the adjacent parabolic fits 
    from separate parabolic fits of lower and high range across the boundary region
    can also handle parabolic to cubic'''    
    # For each find overlapping range
    for i, fitrange in enumerate(overlaps):
        [start,end]=fitrange
        if Backfitparams.iloc[i]['Fittype']=='parabola' and Backfitparams.iloc[i+1]['Fittype']=='parabola':
            A0=Backfitparams.iloc[i]['A'] # grab fit params from lower and upper regions
            B0=Backfitparams.iloc[i]['B']
            C0=Backfitparams.iloc[i]['C']
            A1=Backfitparams.iloc[i+1]['A']
            B1=Backfitparams.iloc[i+1]['B']
            C1=Backfitparams.iloc[i+1]['C']
            thisrange=abs(end-start) # total for crossfading
            for j in range(start, end):
                xval=EDXfile.iloc[j]['Energy']
                yval=(1-(j-start)/thisrange)*(C0+B0*xval+A0*xval**2)+((j-start)/thisrange)*(C1+B1*xval+A1*xval**2)
                EDXfile=EDXfile.set_value(j, 'Backfit',yval)
        if Backfitparams.iloc[i]['Fittype']=='parabola' and Backfitparams.iloc[i+1]['Fittype']=='cubic':
            A0=Backfitparams.iloc[i]['A'] # grab fit params from lower and upper regions
            B0=Backfitparams.iloc[i]['B']
            C0=Backfitparams.iloc[i]['C']
            A1=Backfitparams.iloc[i+1]['A']
            B1=Backfitparams.iloc[i+1]['B']
            C1=Backfitparams.iloc[i+1]['C']
            D1=Backfitparams.iloc[i+1]['D']
            thisrange=abs(end-start) # total for crossfading
            for j in range(start, end):
                xval=EDXfile.iloc[j]['Energy']
                yval=(1-(j-start)/thisrange)*(C0+B0*xval+A0*xval**2)+((j-start)/thisrange)*(D1+C1*xval+B1*xval**2+A1*xval**3)
                EDXfile=EDXfile.set_value(j, 'Backfit',yval)
        if Backfitparams.iloc[i]['Fittype']=='cubic' and Backfitparams.iloc[i+1]['Fittype']=='parabola':
            A0=Backfitparams.iloc[i]['A'] # grab fit params from lower and upper regions
            B0=Backfitparams.iloc[i]['B']
            C0=Backfitparams.iloc[i]['C']
            D0=Backfitparams.iloc[i]['D']
            A1=Backfitparams.iloc[i+1]['A']
            B1=Backfitparams.iloc[i+1]['B']
            C1=Backfitparams.iloc[i+1]['C']
            thisrange=abs(end-start) # total for crossfading
            for j in range(start, end):
                xval=EDXfile.iloc[j]['Energy']
                yval=(1-(j-start)/thisrange)*(D0+C0*xval+B0*xval**2+A0*xval**3)+((j-start)/thisrange)*(C1+B1*xval+A1*xval**2)
                EDXfile=EDXfile.set_value(j, 'Backfit',yval)
    return EDXfile
''' 
i, [fitrange, fitpts, mandminpts, fittype, threshold]=0,Fitregions[0]
'''        
def fitbackgrounds(EDXfile, Fitregions, logmatch, fitcols):
    ''' Background fit for each direct peak(opens source spectrum as EDXfile, 
    fits peak backgrounds above and below using Elemdata, saves background to source csv (overwrites existing fits), 
    also saves linear fit params to logdataframe with position/amplitude/etc;
    Fitregions stores total ev range, background regions, fit type and thresholdfor deriv knockout '''
    # Create temp df to hold and pass linear fit data
    EDXfileName=logmatch.Filename # 
    Backfitparams=pd.DataFrame(columns=fitcols) # empty df to hold all rows from this spectrum
        # all fit regions modify fit region boundaries for this spectrum based on smooth-differentiated peak (2nd deriv, Savgol (poly=2, pts=11))
        # global shifts from smdifpeaks and local shift based on smoothed 2nd derivative 
        # already incorporated into Elemdata values (lower1,2 and upper1,2 fully adjusted)
    
    # loop through and fit all peaks for each element in this spatial area    
    fitranges=[]        
    for i, [fitrange, fitpts, mandminpts, fittype, threshold] in enumerate(Fitregions):
        # get fitranges for cross-fading of backgrounds in overlapping regions
        fitranges.append(fitrange)
        # create new df row for each fitted range
        Backfitparamrow=pd.DataFrame(index=np.arange(0,1),columns=fitcols)
        # modify fit region for this spectrum (eliminate those with high derivative )
        # Threshold level defined by many attempted fits to actual data
        Thisbackfit, backptslist=findfitregion(EDXfile, fitpts, mandminpts, threshold, fitrange, EDXfileName)
        # problem ... this function also can knock out the mandatory min included pts
        # Force counts to zero near origin
        for index,row in Thisbackfit.iterrows():
            if index < 5:
                Thisbackfit=Thisbackfit.set_value(index,'Counts',0)
        # now do parabolic fit over this region (return df with backfit col)
        if fittype=='parabola':
            Thisbackfit, fitparams = fitparabola(Thisbackfit, EDXfileName, fitrange)
            # unpack polynomial fit parameters 
            A=fitparams[0]
            B=fitparams[1]
            C=fitparams[2]    
            # now copy this function over entire range of fit
            lower=int(fitrange.split('-')[0])
            upper=int(fitrange.split('-')[1])
            # TODO these are essentially index numbers (since energy in eV is nearly same as index range)    
            if A!='n/a': # test for successful fit (all fitparams set to n/a)
                try:
                    for i in range(lower,upper): 
                        xval=EDXfile.iloc[i]['Energy']
                        EDXfile.set_value(i,'Backfit',A * xval**2 + B * xval + C) # just set values directly from fit results
                except:
                    print('Problem finding energy value for range ', str(lower),' to ', str(upper),'in ',EDXfileName)
                    print('Fitrange =', fitrange,'A=',A,'B=',B,'C=',C)
            # now store values forthis df row (slower but easy)
            Backfitparamrow.iloc[0]['Fitrange'] = fitrange
            Backfitparamrow.iloc[0]['Fittype'] = fittype
            Backfitparamrow.iloc[0]['A'] = A # parabolic fit params
            Backfitparamrow.iloc[0]['B'] = B
            Backfitparamrow.iloc[0]['C'] = C
            Backfitparamrow.iloc[0]['Backfitpts'] = backptslist # store list of index #s actually used
        # TODO test and incorporate other fit types
        if fittype=='cubic':
            Thisbackfit, fitparams = fitcubic(Thisbackfit, EDXfileName, fitrange)
            # unpack polynomial fit parameters 
            A=fitparams[0] # highest power first
            B=fitparams[1]
            C=fitparams[2]
            D=fitparams[3]
            # now copy this function over entire range of fit
            lower=int(fitrange.split('-')[0])
            upper=int(fitrange.split('-')[1])
            # TODO these are essentially index numbers (since energy in eV is nearly same as index range)    
            if A!='n/a': # test for successful fit (all fitparams set to n/a)
                try:
                    for i in range(lower,upper): 
                        xval=EDXfile.iloc[i]['Energy']
                        EDXfile.set_value(i,'Backfit',A * xval**3 + B * xval**2 + C * xval + D) # just set values directly from fit results
                except:
                    print('Problem finding energy value for range ', str(lower),' to ', str(upper),'in ',EDXfileName)
                    print('Fitrange =', fitrange,'A=',A,'B=',B,'C=',C,'D=',D)
            # now store values forthis df row (slower but easy)
            Backfitparamrow.iloc[0]['Fitrange'] = fitrange
            Backfitparamrow.iloc[0]['Fittype'] = fittype
            Backfitparamrow.iloc[0]['A'] = A # parabolic fit params
            Backfitparamrow.iloc[0]['B'] = B
            Backfitparamrow.iloc[0]['C'] = C
            Backfitparamrow.iloc[0]['D'] = D
            Backfitparamrow.iloc[0]['Backfitpts'] = backptslist # store list of index #s actually used        
        # Set value for subtracted spectral data (nan if nothing in backfit)
        EDXfile['Subdata']=EDXfile['Counts']-EDXfile['Backfit']
        # concatentate single row with log 
        Backfitparams=pd.concat([Backfitparams, Backfitparamrow], ignore_index=True)
        # END OF BACKGROUND FITTING LOOP FOR SINGLE SPECTRUM
    # Cross-fade of adjacent background fit regions (fitranges)
    overlaps=returnoverlaps(fitranges) # returns overlapping ranges as list of [start,end] from string ranges
    # len overlaps one less than 
    EDXfile=crossfadeback(EDXfile, Backfitparams, overlaps)
    
    # assign params that are common to all areas/all peaks into rows of df (copied from original log)
    for index, row in Backfitparams.iterrows():
        Backfitparams.loc[index]['Basename']=logmatch.Basename
        Backfitparams.loc[index]['Filenumber']=logmatch.Filenumber
        Backfitparams.iloc[index]['Filename']=logmatch.Filename
        Backfitparams.iloc[index]['FilePath']=logmatch.FilePath
        Backfitparams.iloc[index]['Sample']=logmatch.Sample
        Backfitparams.iloc[index]['Point']=logmatch.Point
        Backfitparams.iloc[index]['Comments']=logmatch.Comments
        Backfitparams.iloc[index]['Date']=logmatch.Date
        Backfitparams.iloc[index]['Beamkv']=logmatch.Beamkv
        Backfitparams.iloc[index]['Livetime']=logmatch.Livetime
        Backfitparams.iloc[index]['Timeconst']=logmatch.Timeconst
        Backfitparams.iloc[index]['Deadfraction']=logmatch.Deadfraction
    Backfitparams=Backfitparams[fitcols] # put back in original order
    return EDXfile, Backfitparams # df with direct peak fitting info for all areas/ all elements
        
def batchEDXquant(EDXfiles, Fitregionsdf, EDXquantparams, Elements, Backfitlog, Integlog, Peakfitlog, **kwargs):
    ''' Batch quantification of all peaks in Elements list 
    returns df with peak positions, amplitudes, width, energy shift, etc.
    kwargs: redo_backfit-- refit of backgrounds (default false)
			redo_integration-- default true; does integration over after changing backfits (i.e. with interactive refitter)
			clear_old_backfits - if true overwrite all prior background fits
            savegauss-  save gaussian fitting column to csv data file 
			'''  
    fitcols=['Basename', 'Filenumber', 'Filename', 'FilePath', 'Sample', 'Comments', 'Date', 'Point', 'Beamkv', 
    'Livetime','Timeconst','Deadfraction','Fitrange', 'Backfitpts', 'Fittype', 'A', 'B', 'C', 'D', 'Rval', 'Pval', 'Stderr']    
    # blank df for gaussian peak fit results
    peakcols=['Basename', 'Filenumber', 'Point','Filename', 'Filepath', 'Sample', 'Comments', 'Element',
    'Xc', 'Width', 'Peakarea', 'Y0','Rsquared'] # for gaussian peak fits 
    integcols=['Basename','Filenumber', 'Point', 'Filename', 'Filepath', 'Sample', 'Comments', 'Element', 'Energy', 'Shift', 'Rawcounts', 
    'Backcounts', 'Subtractedcounts', 'Adjcounts', '% err', 'Significance' , 'Correctedcounts', 'Errcorrcnts','Kfact','Fullwidth'] # for integration results     
    # TODO any problems if older existing backfitlogs don't match above list
    if kwargs.get('redo_backfit',False):
        # redo background fits even if already processed
        backprocessed=[] # no files to skip
        # Remove prior entries (assumed to be inferior)
        filelist=np.ndarray.tolist(EDXfiles.Filename.unique())
        Backfitlog=Backfitlog[~Backfitlog['Filename'].isin(filelist)]
    else:
        backprocessed=np.ndarray.tolist(Backfitlog.Filename.unique())
    if not kwargs.get('redo_integration', True):
        # Skip redo of integration (defaults true)
        integprocessed=np.ndarray.tolist(Integlog.Filename.unique())
    else:
        integprocessed=[]
        # Remove prior entries (assumed to be inferior)
        filelist=np.ndarray.tolist(EDXfiles.Filename.unique())
        Integlog=Integlog[~Integlog['Filename'].isin(filelist)]
        Peakfitlog=Peakfitlog[~Peakfitlog['Filename'].isin(filelist)]
    # get global values for Elemdata for standard spectra
    Elemdata=findelemregions(Elements, EDXquantparams)
    # Eliminate already processed files from below quant loop (also prevents rerun of excluded files in EDXparamlog.csv)

    for i in range(0,len(EDXfiles)):
        # get ith row from parameters log for subset of selected spe files (i.e. from EDXfiles)
        logmatch=EDXfiles.iloc[i] #contains row with filename and all other parameters from a given spectra 
        logmatch=logmatch.squeeze() # convert/flatten to Series
        # load Auger spe file of interest here
        EDXfileName=logmatch.Filename  # get Auger filename from Series
        EDXfile=openEDX(EDXfileName) # read csv or raw psmsa
        beamkv=int(logmatch.Beamkv) # for spectrum truncation 
        if 'Counts' not in EDXfile:
            print('Counts not present in file ', EDXfileName)
            continue # skip to next spectrum (shouldn't happen)
        # zero out low energy cols for background fitting
        for index, row in EDXfile.iterrows():
            if index <= 5: # 
                EDXfile=EDXfile.set_value(index,'Counts',0)
        if 'Backfit' not in EDXfile: # add this background fit column if not present
            EDXfile['Backfit']=np.nan
        # optional overwrite of all prior backgrounds in csv file
        # Default false for quant rerun .. keep custom backgrounds
        overwrite=kwargs.get('clear_old_backfits', False)
        if overwrite==True: # clear all prior background, subdata and peak fitting regions
            EDXfile['Backfit']=np.nan
            EDXfile['Subdata']=np.nan
            EDXfile['Gauss']=np.nan
            print('Overwrite of backgrounds from', EDXfileName)
        # Sav-gol 2nd deriv column used to guide selection of fitting regions
        if 'Savgol' not in EDXfile: # returns df with this Savgol column added
            EDXfile=makesavgol(EDXfile)  # FUNCT pass full spectrum for given area (saved below)   
        if 'Subdata' not in EDXfile: # add col for subtracted peak data
            EDXfile['Subdata']=np.nan
        # Unpack fitregions from dataframe into list of data ev ranges, 
        # For each ev range we have total fit range, regions which it holding background and fit type

        if EDXfileName not in backprocessed: # skips background refitting if existing one is good
            Fitregions=unpackfitregs(Fitregionsdf, EDXfile) # turn dataframe into list of lists
            # Linear background fits over 4 or 5 energy regions (while dodging major peaks)
            EDXfile, Backfitparams = fitbackgrounds(EDXfile, Fitregions, logmatch, fitcols)
            print('New background fitting for ', EDXfileName)
        else:
            print('Using existing background fit for', EDXfileName)
        savegauss=kwargs.get('savegauss', True) # optional save of gaussian fitting col
        # TODO handling of odd spectra using Elemdatamod (looks up positions in this EDXfile)
        # use global Elemdata, fit major peaks in Subdata with gaussians (and determine shift)
        # Fits peaks, determines shifts and does integration and corrected counts calculations for each element
        if EDXfileName not in integprocessed: 
            EDXfile, Peakfits, Integresults = fitpeaks(EDXfile, Elemdata, logmatch, savegauss, peakcols, integcols)
            print('New peak fits and integrations for ', EDXfileName)
        else:
            print('Using existing integrations and peak fits for', EDXfileName)

        # truncate at beam kV setting 
        beamkv=int(logmatch.Beamkv) # for spectrum truncation 
        EDXfile=EDXfile[EDXfile['Energy']<beamkv]
        # Concat result rows from this spectrum to longer master log (if new values generated)
        if EDXfileName not in backprocessed:
            Backfitlog=pd.concat([Backfitlog,Backfitparams],ignore_index=True)
        Peakfitlog=pd.concat([Peakfitlog,Peakfits],ignore_index=True)
        if EDXfileName not in integprocessed:
            Integlog=pd.concat([Integlog,Integresults],ignore_index=True)
        
        # direct save of modified auger csv with new linear background fits (after all areas processed)
        if '.psmsa' in EDXfileName:
            EDXfileName=EDXfileName.split('.psmsa')[0]+'.csv' # save as .csv not .psmsa
            EDXfile.to_csv(EDXfileName, index=False)
            print(EDXfileName,' saved.')
        elif '.emsa' in EDXfileName:
            EDXfileName=EDXfileName.split('.emsa')[0]+'.csv' # save as .csv not .psmsa
            EDXfile.to_csv(EDXfileName, index=False)
            print(EDXfileName,' saved.')
        elif '.csv' in EDXfileName: # csv probably sum of two psmsa.. just overwrite
            EDXfile.to_csv(EDXfileName, index=False)
            print(EDXfileName,' saved.')
        else:
            print('EDXfile not saved ... file extension is not emsa or psmsa')

    # now sort and reorder master log files for backfits, peakfits, integration of all peaks all spectra   
    Backfitlog=Backfitlog.sort_values(['Basename','Filenumber'], ascending=True)
    for i, col in enumerate(fitcols): # handling older files possibly w/ different colname set
        if col not in Backfitlog.columns:
            Backfitlog[col]=''
    Backfitlog=Backfitlog[fitcols] # Put back in original order
    Peakfitlog=Peakfitlog.sort_values(['Basename','Filenumber'], ascending=True)
    for i, col in enumerate(peakcols): # handling older files possibly w/ different colname set
        if col not in Peakfitlog.columns:
            Peakfitlog[col]=''
    Peakfitlog=Peakfitlog[peakcols] # Put back in original order    
    Integlog=Integlog.sort_values(['Basename','Filenumber'], ascending=True)
    for i, col in enumerate(integcols): # handling older files possibly w/ different colname set
        if col not in Integlog.columns:
            Integlog[col]=''
    Integlog=Integlog[integcols] # Put back in original order
    print('Backfitlog, Integlog, Peakfitlog not saved... manual save required.')
    return Backfitlog, Peakfitlog, Integlog

def EDXquant(EDXfile, Fitregionsdf, EDXquantparams, Elements, **kwargs):
    ''' Single spectrum version of quant loop (normally only used on refitting of failed background fits '''
    
    # create empty dataframe for storing/passing linear fit params (same structure as in fitbackgrounds)
    mycols=['Basename', 'Filenumber', 'Filename', 'FilePath', 'Sample', 'Comments', 'Date', 'Point', 'Beamkv', 
    'Livetime','Timeconst','Deadfraction','Fitrange', 'Backfitpts', 'Fittype', 'A', 'B', 'C', 'D', 'Rval', 'Pval', 'Stderr']    
    # black df for gaussian peak fit results
    mycols2=['Basename', 'Filenumber', 'Point','Filename', 'Filepath', 'Sample', 'Comments', 'Element',
    'Xc', 'Width', 'Peakarea', 'Y0','Rsquared'] # for gaussian peak fits 
    # Now a blank frame for integrated quant results
    mycols3=['Basename','Filenumber', 'Point', 'Filename', 'Filepath', 'Sample', 'Comments', 'Element', 'Energy', 'Shift', 'Rawcounts', 
    'Backcounts', 'Subtractedcounts', 'Adjcounts', '% err', 'Significance' , 'Correctedcounts', 'Errcorrcnts',] # for integration results     
    # Get global values for Elemdata for standard spectra
    Elemdata=findelemregions(Elements, EDXquantparams)
    logmatch=EDXfile.squeeze() # flatten single chosen spectrum to Series
    EDXfileName=logmatch.Filename  # get Auger filename from Series
    EDXfile=openEDX(EDXfileName) # read csv or raw psmsa
    beamkv=int(logmatch.Beamkv) # for spectrum truncation 
    # On re-run can skip many checks present in first run 
    try:
        EDXfile['Backfit']=np.nan
        EDXfile['Subdata']=np.nan
        EDXfile['Gauss']=np.nan
        EDXfile['Subdata']=np.nan
    except:
        print('Problem clearing prior fits.')
    # unpack fit regions from loaded dataframe
    Fitregions=unpackfitregs(Fitregionsdf, EDXfile)
    # main background fitting (guided by allowable points in EDX_backfit_regions)
    EDXfile, Backfits = fitbackgrounds(EDXfile, Fitregions, logmatch)
    savegauss=kwargs.get('savegauss', True) # optional save of gaussian fitting col

    # use global Elemdata, fit major peaks in Subdata with gaussians (and determine shift)
    # Fits peaks, determines shifts and does integration and corrected counts calculations for each element
    EDXfile, Peakfits, Integresults = fitpeaks(EDXfile, Elemdata, logmatch, savegauss)
    # truncate at beam kV setting 
    beamkv=int(logmatch.Beamkv) # for spectrum truncation 
    EDXfile=EDXfile[EDXfile['Energy']<beamkv]
    # Concat result rows from this spectrum to longer master log
    
    # direct save of modified auger csv with new linear background fits (after all areas processed)
    if '.psmsa' in EDXfileName:
        EDXfileName=EDXfileName.split('.psmsa')[0]+'.csv' # save as .csv not .psmsa
        EDXfile.to_csv(EDXfileName, index=False)
    elif '.emsa' in EDXfileName:
        EDXfileName=EDXfileName.split('.emsa')[0]+'.csv' # save as .csv not .psmsa
        EDXfile.to_csv(EDXfileName, index=False)
    elif '.csv' in EDXfileName: # csv probably sum of two psmsa.. just overwrite
        EDXfile.to_csv(EDXfileName, index=False)
    else:
        print('EDXfile not saved ... file extension is not emsa or psmsa')

    # now sort and reorder master log files for backfits, peakfits, integration of all peaks all spectra  
    Backfits=Backfits.sort_values(['Basename','Filenumber'], ascending=True)
    Backfits=Backfits[mycols] # put back in original order
    Peakfits=Peakfits.sort_values(['Basename','Filenumber'], ascending=True)
    Peakfits=Peakfits[mycols2] # put back in original order    
    Integresults=Integresults.sort_values(['Basename','Filenumber'], ascending=True)
    Integresults=Integresults[mycols3] # put back in original order
    return Backfits, Peakfits, Integresults
	
def convertdate(datestring):
    '''Figure out the date format and convert to datetime  '''    
    thismatch=re.search(r'(\d+)-(\w{3})-(\d+)',datestring)
    if type(thismatch)!='NoneType': # found it
        date=datetime.datetime.strptime(datestring,'%d-%b-%Y')
        date=datetime.date.strftime(date,'%m/%d/%y') # returns date as string in common format
        return date
    else:
        print ('Problem converting date format')
        return datestring  # just return in same format as original
    return 
        
def stripparams(header, df, filename):
    '''Generalized parameter finding from header string;  pass string header and template DF with 
    column names, start and end strings and data type  '''
    # use EDXparam.csv as template for data rows
    collist=df.dtypes.index
    for i, col in enumerate(collist):
        start=df.iloc[1][col]
        # skip columns without search string (those not present in header)
        if str(start)=='nan':
            continue
        try:
            tempstr=header.split(start)[1]
            tempstr=tempstr.split('\n#')[0] # end of string is always next line starting with #
            tempstr=tempstr.strip()
            # convert if necessary 
            thistype=df.iloc[2][col] # expected data type as entered in param template
            if thistype=='int':
                tempstr=int(tempstr)
            elif thistype=='float':
                tempstr=float(tempstr)
            # special treatment for dates
            if col=='Date':
                tempstr=convertdate(tempstr)                  
            df=df.set_value(0,col, tempstr)
        except:
            print('Error extracting ', col, 'param data from ', filename)
    # calculate deadtime fraction 
    try:
        # realtime=float(df.iloc[0]['Realtime'])
        # livetime=df.iloc[0]['Livetime']
        deadfract=(df.iloc[0]['Realtime']-df.iloc[0]['Livetime'])/df.iloc[0]['Realtime']
        df=df.set_value(0,'Deadfraction', deadfract)
    except:
        print('Problem calculating deadtime fraction for ', filename)
    df=df.iloc[[0]] # eliminate string searches and pass back single completed df row
    return df

def getfromname(EDXparamrow,filename):
    ''' Extract base name, spectral number and pt # from filename and current data path 
    write all to param log '''
    # default autonumbered naming scheme for Noran emsa
    pattern=re.compile(r'(\(\d+\))') # find spectral number inside parens
    match=re.search(pattern, filename)
    path=os.getcwd() # current data directory
    if match:
        try:
            basename=filename[0:match.start()] # pulls out base name
            specnum=int(filename[match.start()+1:match.end()-1]) # gets spectral number within parentheses
            if '.psmsa' in filename:
                ptnum=int(filename.split('_pt')[1].split('.')[0]) # gets pt number from 'name(10)_pt1.emsa'
            elif '.emsa' in filename:
                ptnum=1
            EDXparamrow=EDXparamrow.set_value(0,'Point',ptnum)
            EDXparamrow=EDXparamrow.set_value(0,'Basename',basename)
            EDXparamrow=EDXparamrow.set_value(0,'Filenumber',specnum)
            EDXparamrow=EDXparamrow.set_value(0,'Filename',filename) # use definite filename instead of title from header
            EDXparamrow=EDXparamrow.set_value(0,'FilePath',path)
            return EDXparamrow
        except:
            print('Problem extracting names from logfile for ', filename)
    else: # unknown naming scheme.. individually named?
        # check for spectrum number e.g. sp2 at end
        teststr=filename.split('.')[0]
        if 'sp' in teststr: 
            match=re.search(r'sp\d+', teststr)
            if match:
                specnum=int(teststr[match.start()+2:]) # ending would be spectrum number
                basename=teststr[0:match.start()]
                ptnum=1
        else: # no interesting pattern found
            basename=teststr # just remove extension
            specnum=1
            ptnum=1
    # write all to df row for return
    EDXparamrow=EDXparamrow.set_value(0,'Point',ptnum)
    EDXparamrow=EDXparamrow.set_value(0,'Basename',basename)
    EDXparamrow=EDXparamrow.set_value(0,'Filenumber',specnum)
    EDXparamrow=EDXparamrow.set_value(0,'Filename',filename) # use definite filename instead of title from header
    EDXparamrow=EDXparamrow.set_value(0,'FilePath',path)
    # Default naming scheme for NORAN point and shoot psmsa
    # other major format is (basename)_sp(specnum)_pt(ptnum)
    ''' TODO finish this generalized extractor
    pattern=re.compile(r'(_sp)|(pt)|\.') # match both sp, pt and . if present
    match=re.finditer(pattern, filename)
    if match:
        thesebreaks=[]    
        for m in match:
            thesebreaks.append(m.start())
        if len(thesebreaks)==2: # basename_sp2.psmsa
            basename=filename[0:thesebreaks[0]]
            specnum=filename[thesebreaks[0]+2:thesebreaks[1]]
            ptnum=1
        if len(thesebreaks)==3:
            basename=filename[0:thesebreaks[0]]
            specnum=filename[thesebreaks[0]+2:thesebreaks[1]]
            ptnum=1
    '''
    return EDXparamrow
    
def getfromlog(EDXparamrow, EDXlogbook):
    ''' Find associated entry in EDX Excel log book  '''
    # need to match base name & spectral number (currently assumes same data entry for all pts within single point & shoot)
    basename=EDXparamrow.iloc[0]['Basename']
    specnum=EDXparamrow.iloc[0]['Filenumber']
    ptnum=EDXparamrow.iloc[0]['Point']
    EDXlogbook=EDXlogbook.replace(np.nan,'')
    mask=EDXlogbook['Basename'].str.contains(basename, case=False)
    match=EDXlogbook.loc[mask]    
    match=match[(match['Filenumber']==specnum)&(match['Point']==ptnum)]
    if len(match)>1:
        print('Problem finding unique log entry for ', basename, specnum)
    if len(match)==0:
        print('No log entry found for ', basename, specnum)
    try:
        EDXparamrow=EDXparamrow.set_value(0, 'Project', match.iloc[0]['Project']) # copy project name
        EDXparamrow=EDXparamrow.set_value(0, 'Sample', match.iloc[0]['Sample']) # copy sample name
        EDXparamrow=EDXparamrow.set_value(0, 'Comments', match.iloc[0]['Comments']) # copy comments
    except:
        print('Problem finding log entry for ', basename, specnum, ptnum)
    return EDXparamrow

def getparams(filelist, reprocess=False):
    '''Main loop that opens files, finds associated parameters from header and combines with logbook and sample name info '''
    # TODO check for consistency between dataframes
    # Check for existing EDX logbook, which may correlate file numbers with sample names
    EDXlogbook, logbool=openorcreatelogbook(filelist)
    if logbool:
        # prints possible errors from existing EDXlogbook to console
        checklogfile(filelist, EDXlogbook)
    try:
        paramtemplate=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\SEMparams.csv') 
    except:
        print('Paramtemplate for finding header params not found.')
        return
    csvbackup('EDXparamlog') # check for existing file, renames as backup and manages other backups
    mycols=['Project','Basename','Filenumber','Filename','FilePath','Sample','Point',
            'Phase','Comments','Date','Time','Beamkv','Livetime','Realtime',
            'Detected','Converted','Stored','Timeconst','Deadfraction']
    # look for existing SEMparam file (for open append)
    if reprocess==True: # create blank
        EDXparamlog=pd.DataFrame(columns=mycols)
        processed=[] # blank list
    else:
        if os.path.isfile('EDXparamlog.csv'): # keep any prior versions as it may have more complete comments/sample names, etc
            EDXparamlog=pd.read_csv('EDXparamlog.csv', encoding='cp437')
            # ensure it has same columns (add if necessary)
            for i, col in enumerate(mycols):
                if col not in EDXparamlog.columns:
                    EDXparamlog[col]=''
            processed=np.ndarray.tolist(EDXparamlog.Filename.unique()) # list of already processed files
        else:
            EDXparamlog=pd.DataFrame(columns=mycols) # asks for no reprocessing, but prior processing not done
            processed=[] # blank list
    for i,filename in enumerate(filelist):
        if reprocess==False and filename in processed:
            continue # skip already-processed files
        else:
            with open(filename, 'r') as file:
                filedata = file.read()
            header=filedata.split('#SPECTRUM    :')[0]
            # all params extracted and placed in single dataframe row
            EDXparamrow=stripparams(header, paramtemplate, filename)
            EDXparamrow=getfromname(EDXparamrow,filename)  # get basename, spectral number and pt # from name
            
            # now find matching sample info, project names from Excel logbook if one exists
            EDXparamrow=getfromlog(EDXparamrow, EDXlogbook)
            EDXparamlog=pd.concat([EDXparamlog,EDXparamrow], ignore_index=True)
    EDXparamlog=EDXparamlog.sort_values(['Basename','Filenumber','Point'], ascending=True) # sort by basename then filenumber
    EDXparamlog=EDXparamlog[mycols] # put back in original order
    EDXparamlog=EDXparamlog.reset_index(drop=True) # reset the index
    EDXparamlog.to_csv('EDXparamlog.csv',index=False)
    return EDXparamlog
    
def checklogfile(filelist, SEMlogbook):
    ''' Checks the user Auger logbook Excel for consistency with the actual data file list from directory
    prints out filenumbers that have a problem to console''' 
    # TODO modify this to allow for different base names
    psmsa=[] # list of filenumbers of psmsa files in directory
    pattern=re.compile(r'(\(\d+\))') # find spectral number inside parens (standard default naming by Noran EDX software)

    for i, name in enumerate(filelist): # deals with 3 cases (psmsa, sem or map)
        match=re.search(pattern, name)
        if match:        
            specnum=match.group(0).split('(')[1]
            specnum=int(specnum.split(')')[0]) # spectral autoincrement number is set in 
            psmsa.append(specnum) # add to psmsa files in dir list     
    logpsmsacombine=[] # combineable psmsa files from excel logbook
    logpsmsa=[] # other file from Excel logbook (sem, map or single psmsa)  
    alllog=[]
    combinelist=SEMlogbook[(SEMlogbook['Lastnumber']>0)] # get file ranges to combine 
    tempdf=SEMlogbook.replace(np.nan, 0)    
    singlelist=tempdf[(tempdf['Lastnumber']==0)] # these are ones not to combine (most for psmsa)
    for i in range(0,len(singlelist)):
        try:        
            num=int(singlelist.iloc[i]['Filenumber'])
            logpsmsa.append(num)
        except:
            print('Error: file numbers must be integers')
    for i in range(0,len(combinelist)):
        try:
            first=int(combinelist.iloc[i]['Filenumber'])
            last=int(combinelist.iloc[i]['Lastnumber'])
            for i in range(first,last+1):            
                logpsmsacombine.append(i)
        except:
            print('Error: file numbers must be integers')
    alllog=logpsmsacombine+logpsmsa
    # alldir=psmsa+semmap
    for i,val in enumerate(logpsmsacombine): # ensure all are in psmsa list
        if val not in psmsa:
            print('File combination error in Excel logfile.  File # ',val,' is not a combinable psmsa file.')
    missingdata=[i for i in alllog if i not in psmsa] # list comprehension for missing data file (but present in excel logbook)
    missingentry=[i for i in psmsa if i not in alllog] # list comprehension for data file with missing log entry (could also convert to sets and compare)
    for i, val in enumerate(missingdata):
        print ('Data file number ', val, ' mentioned in logbook but missing from directory')
    for i, val in enumerate(missingentry):
        print ('Data file number ', val, ' present in directory but missing from logbook')
    # check for duplicate entries in logbook
    myset=set([x for x in alllog if alllog.count(x) > 1])
    for i in myset:
        print('Duplicate entry for file number', i, ' in Excel logbook.')
    return

def savegoodfits(bf, peak, integ, goodfits, overwriteall=False):
    ''' Save subset of good fits to csv (currently by filenumber, overwriteall=False replaces
    those in goodfits but leave other existing entries unchanged '''
    bf=bf[bf['Filenumber'].isin(goodfits)]
    integ=integ[integ['Filenumber'].isin(goodfits)]
    peak=peak[peak['Filenumber'].isin(goodfits)]
    if overwriteall==False: # keeps prior run if not in this good set
        try:
            integ2=pd.read_csv('Integquantlog.csv', encoding='cp437')
            integ2=integ2[~integ2['Filenumber'].isin(goodfits)] # replace those in good list
            integ=pd.concat([integ,integ2],ignore_index=True)
            integ=integ.sort_values(['Filenumber'])
        except:
            print('Existing integquantlog not found in ', os.getcwd())
        try:
            bf2=pd.read_csv('Backfitparamslog.csv', encoding='cp437')
            bf2=bf2[~bf2['Filenumber'].isin(goodfits)]
            bf=pd.concat([bf,bf2],ignore_index=True)
            bf=bf.sort_values(['Filenumber'])
        except:
            print('Existing backfitparamslog not found in ', os.getcwd())
        try:
            peak2=pd.read_csv('Peakfitlog.csv', encoding='cp437')
            peak2=peak2[~peak2['Filenumber'].isin(goodfits)]
        except:
            print('Existing peakfitlog not found in ', os.getcwd())
    peak.to_csv('Peakfitlog.csv', index=False)
    bf.to_csv('Backfitparamslog.csv', index=False)
    integ.to_csv('Integquantlog.csv', index=False)
    return

def loadprocessfiles():
    '''Loads/returns common EDX processing files from working directory 
    files can be excluded if marked in EDX log comments file '''
    fitcols=['Basename', 'Filenumber', 'Filename', 'FilePath', 'Sample', 'Comments', 
        'Date', 'Point', 'Beamkv', 'Livetime','Timeconst','Deadfraction','Fitrange', 
        'Backfitpts', 'Fittype', 'A', 'B', 'C', 'D', 'Rval', 'Pval', 'Stderr']    
    peakcols=['Basename', 'Filenumber', 'Point','Filename', 'Filepath', 'Sample', 
        'Comments', 'Element','Xc', 'Width', 'Peakarea', 'Y0','Rsquared'] 
    integcols=['Basename','Filenumber', 'Point', 'Filename', 'Filepath', 'Sample', 
        'Comments', 'Element', 'Energy', 'Shift', 'Rawcounts', 'Backcounts', 
        'Subtractedcounts', 'Adjcounts', '% err', 'Significance' , 'Correctedcounts', 
        'Errcorrcnts','Kfact','Fullwidth']
    paramcols=['Project','Basename','Filenumber','Filename','FilePath','Sample','Point',
            'Phase','Comments','Date','Time','Beamkv','Livetime','Realtime',
            'Detected','Converted','Stored','Timeconst','Deadfraction']
    if os.path.exists('EDXparamlog.csv'):
        EDXlog=pd.read_csv('EDXparamlog.csv', encoding='cp437')
        start=len(EDXlog)
        EDXlog['Comments']=EDXlog['Comments'].replace(np.nan,'')
        EDXlog=EDXlog[~EDXlog['Comments'].str.contains("exclude",na=False, case=False)]
        if start-len(EDXlog)!=0:
            print('Dropped',str(int(start-len(EDXlog))), 'excluded spectral files.')
    else:
        files=glob.glob('*paramlog.csv')
        if len(files)==1:
            print('Loaded params file', files[0])
            EDXlog=pd.read_csv(files[0], encoding='cp437')
        else:
            print("Couldn't find EDX params file in existing folder.")
            # load blanks to avoid error but cd probably needed
            EDXlog=pd.DataFrame(columns=paramcols) 
            Integlog=pd.DataFrame(columns=integcols)
            Backfitlog=pd.DataFrame(columns=fitcols)
            Peakfitlog=pd.DataFrame(columns=peakcols)
    if os.path.exists('Peakfitlog.csv'):
        Peakfitlog=pd.read_csv('Peakfitlog.csv', encoding='cp437')
    else:
        Peakfitlog=pd.DataFrame(columns=peakcols)
    if os.path.exists('Backfitparamslog.csv'):
        Backfitlog=pd.read_csv('Backfitparamslog.csv', encoding='cp437')
    else:
        Backfitlog=pd.DataFrame(columns=fitcols)
    if os.path.exists('Integquantlog.csv'):
        Integlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
    else:
        Integlog=pd.DataFrame(columns=integcols)
    # Print TEM or SEM to console based on beam kV
    if EDXlog['Beamkv'].max()>30:
        print(EDXlog['Beamkv'].max(),'keV TEM spectra loaded.')
        EDXquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\TEMquantparams.csv', encoding='utf-8')
        Interferences=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\TEM_interferences.csv', encoding='utf-8')
    else:
        print(EDXlog['Beamkv'].max(),'keV SEM spectra and quant params loaded.')
        EDXquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\SEMquantparams.csv', encoding='utf-8')
        Interferences=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\SEM_interferences.csv', encoding='utf-8')
    return EDXlog, Backfitlog, Integlog, Peakfitlog, EDXquantparams, Interferences


def loadcomps():
    ''' Find and load most recent EDXcomp and EDXsummary from cwd '''
    # load EDXcomp file
    filelist=glob.glob('EDXcomp*') 
    EDXcomp=findload(filelist)
    filelist=glob.glob('EDXsumm*') 
    EDXsumm=findload(filelist)
    # Return list of prior real and excluded elements using EDXsumm(ary)
    if len(EDXsumm)>0:
        mycols=EDXsumm.columns.tolist()
        Elements=mycols[:mycols.index('Total')]
        Elemexcl=mycols[mycols.index('Total'):]
        Elements=[i.split('%')[1] for i in Elements if '%' in i]
        Elemexcl=[i.split('%')[1] for i in Elemexcl if '%' in i]
    else:
        Elements=[]
        Elemexcl=[]
    return EDXcomp, EDXsumm, Elements, Elemexcl
    
    
def findload(filelist):
    ''' Find and load most recent file matching glob pattern '''
    if len(filelist)>1:
        dates=[s.split('_')[1].split('.')[0] for s in filelist]
        try:
            dates=[datetime.datetime.strptime(val, "%d%b%y") for val in dates]
            datepos=dates.index(max(dates)) # position of newest date (using max)
            newfile=filelist[datepos]
        except:
            print('File date comparison failed... using first one')
            newfile=filelist[0]
        df=pd.read_csv(newfile)
    elif len(filelist)==1:
        df=pd.read_csv(filelist[0])
    else:
        df=pd.DataFrame()
    return df

def pickelemsGUI(EDXquantparams):
    ''' Quick method of interactively selecting elements for plotting 
    has some hard-coded presets that can be changed using preset dictionaries below
    only elements with info in quant params csv files are selectable
    Note.. only tkinter variables exist after root.destroy
    '''
    # Enter default dictionaries for preset buttons (two currently available)
    preset1={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1,'FeL':1}
    preset2={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1,'O':1}
    # All available elemenst are those with entries in edxquantparams.csv
    elems=np.ndarray.tolist(EDXquantparams.element.unique()) 
    # Subset of elements selected (on) by default
    elemdict={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1}
    root = tk.Tk()
    varlist=[] # list of tkinter IntVars
    for i, col in enumerate(elems): # set up string variables
        varlist.append(tk.IntVar())
        val=elemdict.get(col,0) # set to 1 or 0 based on above default dictionary
        varlist[i].set(val) # set default value based on elemdict
        
    tk.Label(root, text='Select elements for plotting or quant').grid(row=0,column=0)
    
    def choose1():
        ''' Have available preset defaults and adjust checkbox values '''
        # preset1={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1,'Fe2':1}
        # Still have to pass these through as tkinter ints
        for i, col in enumerate(elems): # set up string variables
            val=preset1.get(col,0) # set to 1 or 0 based on above default dictionary
            varlist[i].set(val) # set default value based on elemdict
        root.destroy()

    def choose2():
        ''' Have available preset defaults and adjust checkbox values '''
        # preset2={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1,'Fe2':1}
        # Still have to pass these through as tkinter ints
        for i, col in enumerate(elems): # set up string variables
            val=preset2.get(col,0) # set to 1 or 0 based on above default dictionary
            varlist[i].set(val) # set default value based on elemdict
        root.destroy()
        
    for i, col in enumerate(elems):
        # choose row, col grid position (starting row 1)
        thisrow=i%3+1 # three column setup
        thiscol=i//3
        ent=tk.Checkbutton(root, text=elems[i], variable=varlist[i])
        ent.grid(row=thisrow, column=thiscol)
    # Add preset 1 button (defined above)
    els=list(preset1)
    mystr=', '.join(els)
    c=tk.Button(root, text=mystr, command=choose1)
    lastrow=len(elems)%3+2
    c.grid(row=lastrow, column=0)
    # Add preset 2 button

    els=list(preset2)
    mystr=', '.join(els)
    d=tk.Button(root, text=mystr, command=choose2)
    lastrow=len(elems)%3+3
    d.grid(row=lastrow, column=0)
    # add done button
    e=tk.Button(root, text='done')
    e.bind("<Button-1>", lambda event: root.destroy())
    lastrow=len(elems)%3+4
    e.grid(row=lastrow, column=0)

    root.mainloop()

    elemlist=[] # list of strings with plot number and x or y
    for i, val in enumerate(varlist): # result in normal string, not tkinter StringVar
        if val.get()==1:
            elemlist.append(elems[i]) # add element if box is checked 
    return elemlist