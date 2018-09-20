# -*- coding: utf-8 -*-
"""
Spyder Editor
SEM_batch_conversion script
Extracts important header info into parameter log, designed to read out pertinent header information from all emsa files within a folder.
No need to convert psmsa into csv ... just always strip header when opening
Output into single log file for import into Excel or elsewhere
"""

#%% Load modules
import glob, sys, os # already run with functions 
import pandas as pd
import numpy as np
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX')
import EDX_import_functions as EDXimport
import EDX_quant_functions as EDXquant
import EDX_plot_functions as EDXplot
import EDX_refit_tk_gui as EDXrf
import EDX_quantplotter_tk_gui as EDXqpl
#%% 
# datapath = filedialog.askdirectorypwd
# initialdir="H:\\Research_data", title = "choose data directory")
filelist=glob.glob('*.psmsa')+glob.glob('*.emsa')  # psmsa option

#%% Main file processing loop for emsa or psmsa parameter extraction

# Create parameters log for all SEM-EDX files (autosaved with prior backup) using parameter template
# Checks for existing EDXlogbook correlating filenames w/ sample 
EDXlog= EDXimport.getparams(filelist)
EDXlog= EDXimport.getparams(filelist, reprocess=True) # alt version that reacquires params from existing files
EDXlog.to_csv('EDXparamlog.csv',index=False)

# Creation of jpg images with points/areas superimposed (from .psref and .p_s files).. jpgs directly saved
# returns df with spatial areas (automatically saved w/ backup)
SpatialAreasLog=EDXimport.processpointshoot()
#%% 
# Combine files with same basename/point name (autosaves altered EDXlog with backup)
EDXlog=EDXimport.combineEDX(EDXlog)

#%%  Automated background fitting of SEM-EDX spectra 
# can drop or exclude files here if desired (filter of EDXlog)

# Various ways of slicing up above full parameters log list
EDXfiles=EDXlog
EDXfiles=EDXfiles[0:10][:] # grab first ten rows 
EDXfiles=EDXfiles.iloc[[0]] # select single row
EDXfiles=EDXfiles[EDXfiles['Filenumber'].str.contains("\+",na=False)] # choose only summed files
EDXfiles=EDXfiles[~EDXfiles['Comments'].str.contains("exclude",na=False, case=False)] # choose only summed files
EDXfiles=EDXfiles[EDXfiles['Timeconst']>12500] # backfits fail with small timeconst

#%% Reload of existing files (if reprocessing data) from working directory
EDXlog, Backfitlog, Integlog, Peakfitlog, EDXquantparams, Interferences=EDXimport.loadprocessfiles()

#%%
Elements=EDXimport.pickelemsGUI(EDXquantparams) # interactive element selection
Elements=['S','C','Ca','O','Cr', 'FeL','Fe','Mg','Al','Si'] # meteorites
Elements=['S','C','Ca','O','Cr', 'FeL','Fe','Mg','Al','Si'] # pristine SiC
Elements=['S','C','Ca','O','Cr', 'FeL','Fe','Mg','Al','Si','PtM','PtL','PtL2','Ga','GaL'] # meteorites +FIB artifact
Elements=['N','C','O','FeL','Fe','S','Ca','Mg','Al','Si','Ti'] # refractory analogs
Elements=np.ndarray.tolist(Integlog.Element.unique())# gets prior used element set
Elements.append('PtL2')
Elements.extend(['GaL','PtM', 'Ga','PtL','PtL2'])

# Load energy ranges without peaks for background fitting (various options and can also create custom version)
Fitregionsdf=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEM_backfit_regions.csv', encoding='utf-8')
Fitregionsdf=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEM_backfit_regions_alt.csv', encoding='utf-8')
# Version for pristine grains on graphene
Fitregionsdf=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEM_backfit_regions_pristine.csv', encoding='utf-8')
# TEM version
Fitregionsdf=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\TEM_backfit_regions.csv', encoding='utf-8')
Fitregionsdf=pd.read_csv('SEM_backfit_regions_alt.csv', encoding='utf-8') # local version

# If any modifications were made during quant of this data, load local version stored with data
EDXquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEMquantparams.csv', encoding='utf-8')
EDXquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\TEMquantparams.csv', encoding='utf-8')

#%%
# Run main quant loop (not autosaved so use to_csv save below after checks)
kwargs={}
Backfitlog, Peakfitlog, Integlog= EDXimport.batchEDXquant(EDXlog, Fitregionsdf, EDXquantparams, Elements, Backfitlog, Integlog, Peakfitlog, **kwargs)

# optional kwargs for above command
kwargs.update({'redo_backfit':True}) # default false for redo, redo of integration but not of background fits; no effect on new spectra
kwargs.update({'redo_integration':False}) # defaults true (false allows skip of existing integrations and gauss peak fits
# if quant rerun w/o changing backfits (i.e. after custom mods) skip clear of backfits
kwargs.update({'clear_old_backfits':True}) # default false option to not overwrite all backgrounds in csv files (defaults True)
kwargs.update({'savegauss':False}) # optional save of gaussian fit column into spectrum's csv file; default true

# Find/ Replace subset of files (processed in alternate manner) from above log files.. refit of failed fits  
Backfitlog.to_csv('Backfitparamslog.csv', index=False)
Peakfitlog.to_csv('Peakfitlog.csv', index=False)
Integlog.to_csv('Integquantlog.csv', index=False)

# After successful refit of subset of files, find/replace entries in original logbooks (saves after finishing)
Backfitlog, Peakfitlog, Integlog = EDXimport.replacelogentries(EDXlog, Backfitlog, Peakfitlog, Integlog)
#%% Run interactive EDXrefitter (if any plots, backfit points, etc. are bad)
EDXrf.launch_refitter()

EDXqpl.launch_plotter(os.getcwd())
# Redo integlog, peakfits if any backfits were changed (first reload saved changes from file)
EDXlog, Backfitlog, Integlog, Peakfitlog, EDXquantparams, Interferences=EDXimport.loadprocessfiles()
kwargs={'newback':False,'overwrite':False} # do not refit or overwrite backgrounds... use ones made with interactive refitter
Backfitlog, Peakfitlog, Integlog= EDXimport.batchEDXquant(EDXlog, Fitregionsdf, EDXquantparams, Elements, **kwargs)
# Manual save of peakfitlog and integlog are needed 
Peakfitlog.to_csv('Peakfitlog.csv', index=False)
Integlog.to_csv('Integquantlog.csv', index=False)
#%% PLOTTING to check quality of background fits, peaks, etc.

EDXfiles=EDXlog[0:5] # Selecting subsets of all SEM files 

# Plot counts and background over specified energy range
pkwargs={}
pkwargs.update({'xrange':'0.3-10'}) # optional x range for plot (default is 0-10? )
pkwargs.update({'backfitdf':Backfitlog}) # optional plotting of points used to create background fit 
pkwargs.update({'backfitpts':False}) # skip background pts but include fits 
pkwargs.update({'yrange':[-500,3000]}) # optional y range for plot.. defaults to data range
pkwargs.update({'plotelems':['O','Mg','S','Si', 'Ca', 'Fe', 'FeL']}) # list of elements to label on plots
pkwargs.update({'plotelems':['O','Mg','Si', 'Fe']}) 
pkwargs.update({'PDFname':'counts_report_9Jan18.pdf'}) # alt save name (defaults to countsback_report.pdf)
pkwargs.update({'savgol':True}) # include savgol differentiated plot (default False)
EDXplot.reportcounts(EDXfiles, EDXquantparams, **pkwargs)
EDXplot.reportcounts(EDXlog, EDXquantparams, **pkwargs)

# plot report with subtracted counts and optionally gaussian peak fits (if they exist)
EDXplot.reportSEMpeaks(EDXfiles, plotelems, SEMquantparams, addgauss=True, PDFname='peak_report.pdf') 
# TODO Place center of integration on plot for significant peaks 

# plot subtracted data around major elements including corrected counts 
EDXplot.reportsubdatamajor(EDXfiles, Integquantlog, PDFname='Subcounts_major_report.pdf')

reportcountspeakfits(EDXfiles, Fitregionsdf, plotrange, plotelems, SEMquantparams)
# Now proceed to EDX_quant_main for interference adjustments, \\osition calcs, etc.

# Renaming of troublesome p_s and psmsa files (i.e. containing blanks)
psfiles=glob.glob('*.p_s')
badpsfiles=[i for i in psfiles if '\xa0' in i]
for i, psfile in enumerate(badpsfiles):
    EDXimport.renamePSset(psfile, '\xa0', '_')

train=pd.read_csv('Backfit_training.csv')

