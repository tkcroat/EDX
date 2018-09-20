# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:16:43 2016
TEM-EDX processing functions
@author: tkc
"""

import pandas as pd
import numpy as np
import os, sys, re, glob
if 'C:\\Users\\tkc\\Documents\\Python_Scripts' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts')
    
#%%
    
datapath='H:\Research_data\Stardust\C2010W\TEM'    

df=TEMintegcomp

def findfilename(df, TEMparamlog):
    '''Match info from origin quant log with actual filenames '''
    for index, row in df.iterrows():
        thisbasename=df.loc[index]['Basename']
        thisfilenum=df.loc[index]['Filenumber']
        match=TEMparamlog[(TEMparamlog['Basename']==thisbasename) & (TEMparamlog['Filenumber']==thisfilenum)]
        if not match.empty:
            fname=match.iloc[0]['Filename']
            df=df.set_value(index,'Filename',fname)
        else:
            print('Filename not found for ', thisbasename)
    return df
        
        
    
    
def processoldEDX(excelname, datapath):
    '''Import batch processed results from origin and convert to consistent format with python generated data '''
    if not datapath.endswith('\\'):
        datapath=datapath+'\\'
    df=pd.read_excel(excelname, sheetname='Batch')
    delcols=[int(i) for i in range(18,len(df.columns))] #
    delcols.extend([14,15])
    df.drop(df.columns[delcols], axis=1, inplace=True) # remove all undesirable columns
    mycols=['Filename','Element','Energy','Type','Rawcounts','Low back','High back','Backcounts','Subtractedcounts','Adjcounts','Lc','Ld','Lq',  'Correctedcounts','% err','ID']
    df.columns=mycols # rename to be consistent with other python results 
    # TODO find correct filename from filelist (match on known params)
    df=df.dropna(subset=['Filename']) # drop empty last row
    # TODO calculate shift for each element (energy - ideal energy)
    df['Shift']=''
    df['Significance']=df['Subtractedcounts']/(np.sqrt(df['Backcounts']))
    df['Basename']=''
    df['Sample']=''
    df['Comments']=''
    df['Filenumber']=0
    df['Filepath']=datapath # passed as arg
    df['Errcorrcnts']=''
    df['Point']=1 # only single point for manually collected TEM-EDX
    mycols2=['Basename','Filenumber', 'Point', 'Filename', 'Filepath', 'Sample', 'Comments', 'Element', 'Energy', 'Shift', 'Rawcounts', 
    'Backcounts', 'Subtractedcounts', 'Adjcounts', '% err', 'Significance' , 'Correctedcounts', 'Errcorrcnts',] 
    for index, row in df.iterrows():
        name=df.loc[index]['Filename']
        try:
            basename=name.split('sp')[0]
            tempval=name.split('sp')[1]
            filenum=int(tempval.split('.')[0])
            df=df.set_value(index,'Basename', basename)
            df=df.set_value(index,'Sample', basename) # in most cases same as basename for old TEM spectra
            df=df.set_value(index,'Filenumber', filenum)
        except:
            print("Couldn't find basename, filenumber for ", name)
    df=df[mycols2] # put in standard order
    return df


def generateTEMlog(filelist):
    ''' '''
filelist=TEMintegcomp.Spectrum.unique()
filelist=np.ndarray.tolist(filelist)
# TODO make a matching log (maybe empty maybe pulled from summary) from above batch quant file