# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:43:40 2016

@author: tkc
"""
import pandas as pd
import numpy as np
import os, sys, re, glob
if 'C:\\Users\\tkc\\Documents\\Python_Scripts' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts')
import SEMEDX_quant_functions as SEMquant
import SEMEDX_batch_import_functions as SEMimport
import TEM_EDX_functions as TEM

#%% 
# Load and process TEM-EDX compositions
TEMparamlog=pd.read_csv('TEMparamlog.csv', encoding='cp437')

C2010Wcnts_origin=pd.read_csv('C2010W_SEM_cnts_log_originlab.csv',encoding='cp437')
TEMquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\SEMEDX\\TEMquantparams.csv', encoding='utf-8')

# Read in underlying TEM emsa files from source data directroy 
excelname='C:\\Users\\tkc\\Desktop\\Research\\Stardust_craters\\EDXS\\C2010W_FIBTEM_EDXquant.xlsx'
datapath='H:\Research_data\Stardust\C2010W\TEM'
paramtemplate=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\SEMparams.csv') # shoudl be same as for SEM-EDX

TEMlogbook=SEMimport.openorcreatelogbook(filelist) # blank logbook need for EMSA param input

os.chdir('H:\Research_data\Stardust\C2010W\TEM') # switch to TEM-EDX data directory (or copy into console)
datapath=os.getcwd()
filelist=glob.glob('*.emsa') 

TEMparamlog= SEMimport.getparams(filelist, paramtemplate, TEMlogbook) # get all file params from headers
TEMparamlog.to_csv('TEMEDXparamlog.csv', index=False)

TEMinteglog=TEM.processoldEDX(excelname, datapath) # get origin batch generated quant results
TEMinteglog=findfilename(TEMintegcomp, TEMparamlog) # get correct data file emsa name

# Makes count adjustments and recalculate errors (only acting on peaks with pathological overlap/ interferences)
# values are in Adjcounts column (Subtractedcounts unaffected)
Interferences=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\SEMEDX\\TEM_interferences.csv', encoding='utf-8')
TEMinteglog=SEMquant.recalcadjbatch(TEMinteglog, Interferences, TEMquantparams) 
Integquantlog=SEMquant.recalcadjbatch(Integquantlog, Interferences, TEMquantparams) 
TEMinteglog=recalcadjbatch(TEMinteglog, Interferences, TEMquantparams) 

Integquantlog=SEMquant.recalccorrcounts(Integquantlog, TEMquantparams) # recalculcate corrected counts based on current kfactors
TEMinteglog.to_csv('TEMinteglog_origin.csv', index=False)

# Now can calculate compositions from TEM-EDX 
Elements=['Ca', 'Fe', 'Si', 'Mg', 'S','Ni'] # Fe2 is usually Fel line
Elements=['Ca', 'Fe', 'Si', 'Mg', 'O']

# Now calculate compositions 
C2010TEMcomp=SEMquant.calccomp(TEMfiles, Integquantlog, Elements, sigthreshold=2)
C2010TEMcomp.to_csv('C2010WTEM_comp_origin.csv', index=False)

# recalculation of 
#%% OTHER STUFF
# TODO  Choose subset of spectra and combine into averaged spectrum 
TEMfiles=TEMparamlog[(TEMparamlog['Basename']=='f1-2c2') & TEMparamlog['Filenumber'].isin([3,7,9,10])]
TEMfiles=SEMimport.combineEDX(TEMfiles)
