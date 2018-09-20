# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 15:07:07 2016

@author: tkc
"""

import pandas as pd
import os, sys, glob
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX')
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Utilities' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Utilities')
import EDX_import_functions as EDXimport
import EDX_plot_functions as EDXplot
import EDX_quant_functions as EDXquant
import EDX_quantplotter_tk_gui as EDXqpl

import pandas_utilities as pdutil

#%% 
os.chdir('c:\\Temp\\TEM')
os.chdir('c:\\Temp\\EDX')
# Loading of prior parameter logs, background fits, integration quant results, peakfits 
EDXlog, Backfitlog, Integlog, Peakfitlog, EDXquantparams, Interferences=EDXimport.loadprocessfiles() # load all existing from current w/d
EDXcomp=pd.read_csv('EDXcomp_26Jan18.csv', encoding='cp437')
# load existing EDX quant and summary files (on rerun)
EDXcomp, EDXsumm, Elements, Elemexcl=EDXimport.loadcomps()

# Alt/manually load custom SEM or TEM EDX quant parameters (kfactors)... defaults loaded above
EDXquantparams=pd.read_csv('EDXquantparams.csv', encoding='utf-8') # reading local version
Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
EDXlog=pd.read_csv('EDXparamlog.csv',encoding='cp437')
#%% ADJUSTMENTS FOR PATHOLOGICAL PEAK OVERLAPS (if present)
# Reload list of pathological interferences (correct version usually auto-loaded above)
Interferences=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\TEM_interferences.csv', encoding='utf-8')
Interferences=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEM_interferences.csv', encoding='utf-8')
# Calculate and return adjusted counts column to account for pathological overlap in case of interferences
# this also recalculates corrcnts and associated error
Integlog=EDXquant.recalcadjbatch(Integlog, Interferences, EDXquantparams) 

# Now recalculate corrected counts (using adjcounts if present)
# Recalculate corrected counts/errors after k-factor changes
Integlog=EDXquant.recalccorrcounts(Integlog, EDXquantparams)
Integlog.to_csv('Integquantlog.csv', index=False) # manually save after changes made

# Interactive spectral viewer and quant summary
EDXqpl.launch_plotter(os.getcwd())

# Tk gui interface w/ single EDX window plotter
EDXplot.EDXplot_gui(EDXlog, Elements, Backfitlog, EDXquantparams)
#%% Proceed to compositional calculations
# Select elements of interest
Elements=np.ndarray.tolist(Integlog.Element.unique())# gets prior used element set

Elements=EDXimport.pickelemsGUI(EDXquantparams) # use GUI checkbox
Elements=['Ca', 'Fe', 'FeL', 'Si', 'Mg', 'S'] # Fe2 is Fe-L line
Elements=['Ca', 'Fe', 'Si', 'Mg', 'Al', 'S','O','Ni']
Elements=['Ca', 'Fe', 'Si', 'Mg', 'Al', 'S','Ni']
Elements=['Fe', 'Ni', 'S']
Elements.remove('Pbl')
Elements.append('Na')
Elements.append('Cu')
Elemexcl=['Ga','GaL', 'PtM','Al','C','Cu','Cr','Na','FeL']
Elemexcl=['FeL','Cr']

Elemexcl=['Fe','Cr']
Elemexcl.append('S')
Elemexcl.extend(['Ga','GaL','PtM'])
Elemexcl.remove('C')
real=[i for i in Elements if i not in Elemexcl]
# Get elements from prior quant
Elements=[i.replace('%','') for i in EDXcomp.columns if i.startswith('%')]
Elements2=[i.replace('%','') for i in EDXcomp2.columns if i.startswith('%')]

# Calculate compositions for selected elements
EDXcomp=EDXquant.calccomp(EDXlog, Integlog, Elements, sigthreshold=2)

# Same as calccomp but returns subtracted counts for each elem
EDXcomp=EDXquant.getcountscomp(EDXlog, Integlog, Elements, sigthreshold=2)
tempcomp=EDXquant.getcountscomp(EDXlog, Integlog, Elements, sigthreshold=2)

EDXcomp.to_csv('EDXcomp_12Apr18.csv',index=False, float_format='%.1f')
EDXcomp=pd.read_csv('EDXcompfull_28Mar18.csv', encoding='cp437')
EDXcomp=pd.read_csv('EDXcomp_5Jan18.csv', encoding='cp437')
EDXcomp2=pd.read_csv('EDXcomp_13Nov17.csv', encoding='cp437')
EDXcomp['Comments']=EDXcomp['Comments'].replace(np.nan,'')

# Quick output statistical description of at.% to console
EDXquant.describecomps(EDXcomp, Elements)

# Print filenames of outliers in trace/minor elements
# check these in pdf report to see if real or erroneous
np.ndarray.tolist(EDXcomp[EDXcomp['%Fe']>10].Filename.unique())
# plotter to check these files
EDXplot.EDXplot_gui(EDXlog, Elements, Backfitlog, EDXquantparams)

# Abbreviated compositional summary (at.%)
EDXsumm=EDXquant.compsummary(EDXcomp, Elements, Elemexcl)

EDXsumm.to_csv('EDXsummary_12Apr18.csv', index=False, float_format='%.1f')

EDXsumm=pd.read_csv('EDXsummary_26Jan18.csv')

# Sync selected cols between two csv log files (i.e. after manual entries)
pdutil.synclogs('EDXparamlog.csv','EDXsummary_26Jan18.csv',['Phase','Comments'])
pdutil.synclogs('EDXcomp_12Apr18.csv','EDXsummary_12Apr18.csv',['Phase','Comments','Sample'])
pdutil.synclogs('EDXparamlog.csv','EDXsummary_12Apr18.csv',['Sample','Phase','Comments'])

# Plots to assist in setting k-factors (plot of correctedcounts w/ outliers return)
compdata, outliers=elemcompareplot(Integlog, 'PtL', 'PtM',thresh=0.1, errbars='xy')

# Print each matching composition + summary description
kwargs={}
kwargs.update({'string':'area5*xt4|area5*px1'}) # optional substring search for filename
kwargs.update({'string':'area8*xt*'})
kwargs.update({'%Si':'<40'})
tempcomps=EDXquant.printcomps(EDXcomp, Elements)
tempcomps=EDXquant.printcomps(tempcomp, Elements, **kwargs)
tempcomps=EDXquant.printcomps(EDXcomp, Elements, **kwargs)
tempcomps.to_csv('tempcomps.csv', index=False)

# Drop elements from EDXcomp (those never present)
EDXcomp=EDXquant.dropelement(EDXcomp, 'K')

# Retrieve set of elements present (superceded by loadcomps?)
Elements=EDXquant.getelements(EDXcomp)

# Alter element selection before renormalization
Elements=EDXquant.changeelemsGUI(EDXquantparams, Elements)
# Renormalize EDXcomp after changing included elements
EDXcomp2=calccomp(EDXlog, Integlog, Elements, sigthreshold=2)
EDXcomp.to_csv('EDXcomp_metalsbasis.csv', index=False)

# Get TEM image info (generally same folder as EMSA)
filelist=glob.glob('*.dm3')
TEMimagelog=EDXquant.getdm3params(filelist)
TEMimagelog=getdm3params(filelist)
TEMimagelog.to_csv('TEMimagelog.csv',index=False)
TEMimagelog=pd.read_csv('TEMimagelog.csv', encoding='cp437')

# Load diffraction data/ phase informatio
difflog=pd.read_excel('C:\\Users\\tkc\\Desktop\\Research\\Miscellaneous\\Fayetteville_diff_indexing.xlsx', sheetname='Summary')
difflog=difflog.dropna(subset=['a'])
# Make names consistent between diff and EDX data
EDXquant.crosschecknames(TEMimagelog,EDXcomp, difflog)
EDXquant.crosschecknames(dm3log,EDXlog, pd.DataFrame()) # if no diffraction data 
# Get subset with identified phase from diff
EDXcomp=EDXquant.findphase(EDXcomp, difflog)
# Quick plotting of compositional subsets 

# Plot
EDXfiles=EDXlog[0:5] 

kwargs={}
kwargs.update({'xrange':'0.3-10'}) # optional x range for plot (default is 0-10? )
kwargs.update({'plotbackfits':False}) # default true,  plot of background fits cols
kwargs.update({'plotbackpts':Backfitlog}) # Optional plotting of points used to create background fit 
kwargs.update({'yrange':[-500,3000]}) # optional y range for plot.. defaults to data range
kwargs.update({'plotelems':['O','Mg','S','Si', 'Ca', 'Fe', 'FeL']}) # list of elements to label on plots
kwargs.update({'PDFname':'counts_report_30Jan18.pdf'}) # alt save name (defaults to countsback_report.pdf)
kwargs.update({'savgol':True}) # include savgol differentiated plot (default False)
EDXplot.reportcounts(EDXlog, EDXquantparams, **kwargs)


# TODO spectral plot and compositional plot

# Make a standard ternary plot
kwargs={}
kwargs.update({'colorgroup':'Phase'}) # groupby on column name and plot by color
kwargs.update({'colorgroup':'Comments'})
kwargs.update({'colorgroup':'Sample'}) 
kwargs.update({'symbolsize':100})
kwargs.update({'symboltype':'v'})
kwargs.update({'fontsize':30}) # font size for axis labels 
kwargs.update({'title':'SEM-EDX'})
kwargs.update({'plotunknown':False}) # include compositions with undetermined phase (default true)
ternelems=['Fe', 'Si', 'Mg+Ca']
ternelems=['Fe', 'Si+Al', 'Mg+Ca']
ternelems=['Fe', 'S', 'Mg+Ca']
EDXslice=EDXcomp[EDXcomp['Filename'].str.contains('wide',case=False)]
EDXslice=EDXcomp[~EDXcomp['Comments'].str.contains('Pt cap',case=False)]

EDXplot.plotternary(EDXcomp, ternelems, **kwargs)
EDXplot.plotternary(EDXslice, ternelems, **kwargs)

EDXcomp=EDXcomp[EDXcomp['Phase']!='mix'] # drop certain phases 
for index, row in EDXcomp.iterrows():
    if row.Phase=='FeS':
        EDXcomp=EDXcomp.set_value(index,'Phase','FeNiS')

# Principal Components Analysis 

# Sort and drop worst of the duplicates
outputduplicates(C2010WEDXcomp)

scattercompplot(C2010Woldquant, C2010WEDXcomp, Elements, basis=False)

# COMPARING COMPOSITIONS CALCULATED IN DIFFERENT WAYS 
compdata, outliers=EDXplot.scattercompplot(EDXcompFe, EDXcompFe2, Elements, joinlist=['Filename'], thresh=0.1, basis=True, errbars='none')
compdata, outliers=scattercompplot(EDXcompFe, EDXcompFe2, Elements, joinlist=['Filename'], thresh=0.1, basis=True, errbars='none')

outliers=outliers[outliers['Element']=='Fe']
outliers.to_csv('outliers.csv', index=False)
negoutliers=outliers[outliers['Resid']<0]
posoutliers=outliers[outliers['Resid']>0]

negoutlier.to_csv('negoutlier.csv', index=False)
posoutlier.to_csv('posoutliers.csv', index=False)

# plot of outlier spectra 
# Plot counts and background over specified energy range
plotrange='0-6.0' # energy range in keV 
plotrange='0.1-7.0'
yrange=[-500,1500]
yrange=[] # defaults to plot range 
plotelems=['Mg','S','Si', 'Ca', 'Fe', 'Fe2'] # list of elements to label on plots

# Legacy version
negoutlier=pd.merge(negoutliers, EDXparamlog, how='inner', on=['Filename'], suffixes=('','_b'))
posoutlier=pd.merge(posoutliers, EDXparamlog, how='inner', on=['Filename'], suffixes=('','_b'))
compdatafull=pd.merge(compdata, EDXparamlog, how='inner', on=['Filename'], suffixes=('','_b'))
compdatafull.to_csv('compdata_full.csv',index=False)

# saving various files 
C2010WEDXcomp.to_csv('C2010W_SEMcomp_match.csv', index=False)

# Assemble larger dataset (paramlogs and integlogs) from all in subfolders
paramloglist=glob.glob('./**/EDXparamlog.csv', recursive=True)
integloglist=glob.glob('./**/Integquantlog.csv', recursive=True)
Masterparams, Masterinteglog=EDXquant.assembledataset(paramloglist,integloglist)
Masterparams.to_csv('Masterparamslog.csv', index=False)

Masterparams=pd.read_csv('Masterparamslog.csv', encoding='cp437')
Masterparams.Timeconst.hist() # histogram of timeconst values
Masterparams.Beamkv.hist()# histogram of beamkv

# now filter and process multi-project dataset
thismask=Masterparams['FilePath'].str.contains('south', case=False, na=False)
Masterparams=Masterparams.loc[~thismask] # non-projectile refractory mix craters

My7600set=Masterparams[(Masterparams['Beamkv']==10) & (Masterparams['Timeconst']==7600)]

# Getting compositional datasets 
datsouth15kv=EDXquant.getcountscomp(Southfoilparamlog, Masterinteglog, Elements, sigthreshold=2)
datnorth10kv7600=EDXquant.getcountscomp(north10kv7600, Masterinteglog, Elements, sigthreshold=2)

datsouth15kv.to_csv('Southfoil_15kv_compcntdata.csv', index=False)

# Example of getting 10 most Fe rich crater spectra (for subsequent individual plot reports)
Ferich10kv=datnorth10kv.sort_values(['Fecnts'], ascending=False)
Ferich10kv=Ferich10kv.head(10)

# Scatter plotting 

# Batch plot reports of underlying spectra 
plotrange='0-7.0' # energy range in keV 
plotrange='0.1-7.0'
yrange=[-500,5000]
yrange=[] # defaults to plot range 
plotelems=['O','Mg','Al','Si','Ti','Fe', 'Fe2'] # list of elements to label on plots
reportcounts(Ferich10kv, plotrange, yrange, plotelems, EDXquantparams, PDFname='Ferich10kV_counts_report.pdf')
