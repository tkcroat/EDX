# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 13:59:23 2017

@author: tkc
"""

#%% Alternate uncommon processing functions

# Separately reload fits and quant results (if already run)
EDXlog=pd.read_csv('EDXlog.csv', encoding='cp437')
Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
Backfitparamslog=pd.read_csv('Backfitparamslog.csv', encoding='cp437')
Peakfitlog=pd.read_csv('Peakfitlog.csv', encoding='cp437')

# Calculate SEM-EDX compositions (include if Ld or greater) for older files (originlab mostly)
# Calccomposition acts on older count files (calccomp on python imported & fit files)
C2010WEDXcomp=EDXquant.calccomposition(C2010WEDXcnts, EDXquantparams, Elements, sigthreshold=2)
C2010WEDXcomp=calccomposition(C2010WEDXcnts, EDXquantparams, Elements, sigthreshold=2)
SEMlog_origin=pd.read_csv('C2010W_SEMlog_originlab.csv', encoding='cp437')
C2010WEDXcnts=pd.read_csv('C2010W_SEMcounts.csv', encoding='cp437')
C2010Woldquant=pd.read_csv('C2010W_EDXquant.csv', encoding='cp437')
C2010Wsemcomp=pd.read_csv('C2010W_SEMcomp.csv', encoding='cp437')
C2010WTEMcomp=pd.read_csv('C2010W_FIBTEMcomp.csv', encoding='cp437')

# Section for processing Brendan TEM data 
# get sample name from filename using BH name convention 
EDXlog['Sample']=EDXlog['Filename'].str.split('_').str[1]
EDXlog['Sample']=EDXlog['Sample'].str.split('_').str[0]
EDXlog.to_csv('EDXparamlog.csv', index=False)
print(",".join(np.ndarray.tolist(EDXlog['Sample'].unique())))
crosschecknames(dm3log,EDXlog, pd.DataFrame()) # Check for naming convention problems

dm3log['Filename']=dm3log['Filename'].str.split('_').str[0]
np.ndarray.tolist(dm3log.Filename.unique())

#%% Old method of dealing with failed fits (not in batch mode) 
# better way now is to use tk interactive GUIrefitter
badfits=[3,8,10,12,13]
goodfits=[i for i in EDXfiles.Filenumber.tolist() if i not in badfits]
# Saves only good fits and removes bad ones
EDXimport.savegoodfits(Backfitparamslog, Peakfitlog, Integquantlog, goodfits, overwriteall=True)
EDXfiles=EDXfiles[EDXfiles['Filenumber'].isin(goodfits)]

# Fix bad background fits via exclusion of stray elements (i.e. Pt)
EDXfile=EDXlog[EDXlog['Filenumber']==3] # choose single bad spectrum
Backfits, Peakfits, Integ = EDXimport.SEMquant(EDXfile, Fitregionsdf, EDXquantparams, Elements, **kwargs)
interfere=EDXimport.pickelemsGUI(EDXquantparams)
kwargs={'intelems':interfere} # list of interfering elements
# TODO some better means of dealing with interfering peaks
# TODO get replace log entries for single file working
Backfitlog, Peakfitlog, Integlog = EDXimport.replacelogentries(EDXfiles, Backfitlog, Peakfitlog, Integlog)

# Deal separately with bad fits... normally with alternate background regions 
EDXplot.plotderivthresh(EDXfiles, Backfitparamslog, Fitregionsdf, plotrange) # mostly used to refine fitting process
