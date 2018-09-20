# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 12:10:57 2017

@author: tkc
"""
import numpy as np
import pandas as pd
import re, os, glob

def getstagefromapf(df, matchstr, apffilename):
    '''Pass NSScsvparams and get matching stage positions from .apf file and 
    add to NSS params file (uses matching string in case
    NSS file contains multiple arrays '''
    if 'Xpos' not in df:
        df['Xpos']=np.nan
    if 'Ypos' not in df:
        df['Ypos']=np.nan
    tempdf=df[df['SIname'].str.contains(matchstr)]
    with open(apffilename, "r") as infile:
        wholefile = infile.read()
        # find Pt: instances
        pattern=re.compile(r'Pt: ')
        matches=re.finditer(pattern,wholefile) # 
        ptbreaks=[] # finds string start position of each "Pt: xpos ypos" line 
        ptitems=[]
        for m in matches:
            ptbreaks.append(m.start())
        for i in range(0,len(ptbreaks)-1):
            ptitems.append(wholefile[ptbreaks[i]:ptbreaks[i+1]])
        ptitems.append(wholefile[ptbreaks[-1]:]) # get final position
        if len(ptitems)!=len(tempdf):
            print('Mismatch between # pts in apf file and associated number in NSScsvparams.')
            return df # return unmodified
        for i, val in enumerate(ptitems):
            xval=float(ptitems[i].split(' ')[1].strip())
            yval=float(ptitems[i].split(' ')[2].strip())
            # set using index in original df of ith filtered row
            df=df.set_value(tempdf.index[i], 'Xpos', xval) 
            df=df.set_value(tempdf.index[i], 'Ypos', yval) 
    return df

def renamexrayfiles():
    '''Find csv xray files in directory and find replace common elements  '''
    filelist=glob.glob('* *.csv') # x-ray files containing space char
    elemdict={'Si K':'Si','Mg K':'Mg','In L':'In','Ti K':'Ti','S K':'S','O K':'O','Na K':'Na','Ga L':'Ga',
             'Ca K':'Ca', 'Pt M':'Pt','Fe L':'FeL','Fe K':'FeK'}
    for i, name in enumerate(filelist):
        # use regex? 
        pattern=re.compile(r'_([a-z]+) (\w+).csv', flags=re.IGNORECASE)
        
        match=re.search(pattern, name)
        match.group(0)
        start=name.split('_')[0]
        mid=name.split('_')[1]
        mid=mid.split('.csv')[0]
        mid=elemdict.get(mid,mid) # replace if found
        newname=start+'_'
        
    

    newlist=filelist # set of new names
        
    for i, key in enumerate(elemdict):
        
        newval=elemdict.get(key,'')
        
        newlist=[str.replace(key)]
        
        
def findelems(Stagepositions):
    ''' Check xray files and see which elements are available... assume same element subset for all '''
    xrayfile=Stagepositions.iloc[0]['Xrayfile']
    tempstr='*'+xrayfile+'*'
    filelist=glob.glob(tempstr)
    elems=[s.replace(xrayfile+'_','').replace('.csv','').strip() for s in filelist]
    elemdict={'Si K':'Si','Mg K':'Mg','In L':'In','Ti K':'Ti','S K':'S','O K':'O','Na K':'Na','Ga L':'Ga',
                 'Ca K':'Ca', 'Pt M':'Pt','Fe L':'FeL','Fe K':'FeK'}
    elems=[s.replace('.csv','') for s in elems]
    
def assemblexraymap(Stagepositions, elem, pixsizes=[512,384]):
    ''' Assemble larger SEM xray map based on shifts found using image montage 
    inputs are csv files for each element extracted using NSS software (no hack of .si yet)
    xpix and ypix from stagepos are upper left insertion position for new arra
    x-ray maps are normally at different resolution than images... adjust for this '''
    arraylist=[] # load all the image arrays
    for index, row in Stagepositions.iterrows():
        # csv filename for ith np array for elem
        fname=Stagepositions.loc[index]['Xrayfile']+'_'+elem+'.csv' 
        if os.path.isfile(fname):
            # direct load of image csv data into numpy array
            arraylist.append(np.loadtxt(fname, delimiter=',', skiprows=5)) 
        else:
            print('Problem loading array ', fname)    
    # adjust xpix and ypix for resolution of x-ray image.. gray images are 512 x 384
    height, width=arraylist[0].shape
    Stagepositions['Xpix']=Stagepositions['Xpix']*(width/pixsizes[0])
    Stagepositions['Ypix']=Stagepositions['Ypix']*(height/pixsizes[1])
    # need to round first for proper handling of negative decimals
    Stagepositions['Xpix']=Stagepositions['Xpix'].round(decimals=0)
    Stagepositions['Ypix']=Stagepositions['Ypix'].round(decimals=0)
    Stagepositions['Xpix']=Stagepositions['Xpix'].astype(int)
    Stagepositions['Ypix']=Stagepositions['Ypix'].astype(int)

    # Calculate total array size (Xpix and Ypix are upper left of each array)
    ymax=int(Stagepositions.Ypix.max())+height
    ymin=int(Stagepositions.Ypix.min())
    xmin=int(Stagepositions.Xpix.min())
    xmax=int(Stagepositions.Xpix.max())+width
    # adjust for different resolutions of 
    xsize=xmax-xmin
    ysize=ymax-ymin
    # make blank 2D array of appropriate size 
    comboarray=np.empty([ysize, xsize]) # row-column indexing- reverse of cartesian x,y
    if ymin<0:
        Stagepositions['Ypix']=Stagepositions['Ypix']-ymin
    for i, thisarr in enumerate(arraylist):
        xpos=Stagepositions.iloc[i]['Xpix']
        ypos=Stagepositions.iloc[i]['Ypix']
        for j in range(0,height):
            for k in range(0, width):
                comboarray[j+ypos,k+xpos]=thisarr[j,k]
    return comboarray