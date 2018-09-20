# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 12:51:54 2016

@author: tkc
"""
import os, re, glob
import pandas as pd 
import numpy as np
from math import factorial # used by Savgol matrix
from io import StringIO
import datetime
from scipy import optimize

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

def unpackfitregs(df):
    ''' Loaded data frame has ev range, list of background regions and fit type
    unpack from dataframe into list of each 
    1) fitrange is total ev range (i.e.0-100), 2 )backgroundregions are index #s (or energy in eV) of regions commonly without peaks
    3) fittype (mostly parabola) and 4) threshold for derivative knockout
    add overlap range between adjacent fits''' 
    Fitregions=[]
    overlapregs=[] # list of lists containing adjacent overlapping index ranges
    # TODO test to ensure that normal eV range corresponds to range of indices
    for i in range(0,len(df)):
        tempstr=df.iloc[i]['Backgroundregs']
        indexrange=rangefromstring(tempstr) # converts string describing range to actual range
        Fitregions.append([df.iloc[i]['Fitrange'],indexrange,df.iloc[i]['Fittype'], df.iloc[i]['Threshold']])
    return Fitregions

def makesavgol(df):
    '''Perform python smooth-diff used to guide selection of background regions for SEM-EDX spectra
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
    
def openSEM(SEMfileName):
    '''Open csv as dataframe if it exists or if not strip header from psmsa/emsa and import as dataframe ''' 
    csvname=str(SEMfileName.split('.')[0])+'.csv'
    try:
        SEMfile=pd.read_csv(csvname, encoding='cp437')
    except: # if csv doesn't exist, just open/strip psmsa
        with open(SEMfileName, 'r') as file:
            filedata = file.read()
        filedata =filedata.split('#SPECTRUM    :')[1]
        filedata =filedata.split('#ENDOFDATA   : ')[0]
        thisdata=StringIO(filedata)        
        SEMfile=pd.read_csv(thisdata)
        try:        
            SEMfile=SEMfile.drop(SEMfile.columns[[2]], axis=1) # drop erroneous 3rd column if present
        except:   
            print('') # ignore error if 3rd column not present
        SEMfile.columns=['Energy','Counts']    
    return SEMfile # should return data as pandas dataframe

def fitparabola(df, SEMfileName):
    '''Pass appropriate chunk from Auger spectral dataframe, perform polynomial/parabola fit
    return chunk with backfit column added '''
    xcol=df['Energy']
    ycol=df['Counts'] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
        A,B,C=np.polyfit(xcol, ycol, 2)
    except: # deal with common problems with linregress
        print('Fitting error from ', "{0:.2f}".format(df.Energy.min()),'to ',"{0:.2f}".format(df.Energy.max()), ' in file ', SEMfileName)
        fitparams=('n/a','n/a','n/a') # return all n/a
        return df, fitparams
    fitparams=(A, B, C) # tuple to return coeffs of 2nd order poly fit
    for index,row in df.iterrows(): # write this fit into this chunk of data (redundant?)
        xval=df.loc[index]['Energy']
        yval= A * xval**2+ B * xval + C
        df=df.set_value(index, 'Backfit', yval)
    return df, fitparams

def findfitregion(df, fitregion, threshold):
    '''Passing single list of allowable index #s for background fits (no duplicates) 
    remove those with high from list of allowable indices any that show high smoothed-derivatives (i.e. not good for background fitting '''
    Backfitdf=df.ix[[x for x in fitregion]] # filter out those not in allowable background ranges
    # these are loaded from SEM_backfit_regions.csv
    Backfitdf=Backfitdf.dropna(subset=['Counts']) # drops above (set to na by ix)
    # now additionally filter out those with derivative above threshold value
    Backfitdf=Backfitdf[(Backfitdf['Savgol']<threshold) & (Backfitdf['Savgol']>-threshold)]
    return Backfitdf

def findelemregions(Elements, SEMquantparams):
    ''' Takes element string and returns standard Elemdata for each elem symbol containing params 
    needed for peak finding and quant 
    tuple for integ peak is symbol, ideal peak index #, and integ kfactor
    don't apply energy shifts here... apply later when doing integrate''' 
    Elemdata=[]
    try:
        for i, elem in enumerate(Elements):
            # find row in AESquantparams for this element
            thiselemdata=SEMquantparams[(SEMquantparams['element']==elem)]
            thiselemdata=thiselemdata.squeeze() # series with this elements params
            
            # integ peak position value is relative to negpeak in smooth-diff (i.e. -5 is 5 eV below ideal negpeak)
            idealindex=int((thiselemdata.energy+.01)*100) # ideal index value of SEM-EDX peak from energy in keV
            kfact=thiselemdata.kfactor # typical sensitivity k-factor associated with element for integration
            errkfact=thiselemdata.errkfact 
            mass=thiselemdata.mass
            # full peak width in keV from SEMquantparams (usually 0.15keV or 15 channels at 0.1eV/chan)
            width=int(((thiselemdata.fullwidth*100)-1)/2) # integration width in channels for direct integration for this element
            # total # of channels in AESquantparams but include n-1/2 channels on either side of peak center (usually width is 8 channels)
            
            #Elemdata is a list (of length number of elements) containing length5 tuples
            elemtuple=(elem, idealindex, width, kfact, errkfact, mass) # add tuple with info for this element
            Elemdata.append(elemtuple) # now contains proper limits on fitting regions 
    except:
        print('Quant parameters are not properly loaded.')
    return Elemdata 
    
def fitgauss(df, halfwidth, elem, SEMfileName, savegauss=True):
    ''' Gaussian fit of direct peaks (pass SEMfile just around peaks region
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
        print('Gaussian fitting error for', elem, ' peak in file ', SEMfileName)
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

def fitpeaks(SEMfile, Elemdata, logmatch, savegauss=True):
    ''' Gaussian fit of major peaks in single spectrum, shift is list of energy shifts of negpeak (same order as Eledata (opens source spectrum as SEMfile, 
    fits peak backgrounds above and below using Elemdata, also saves linear fit params to logdataframe with position/amplitude/etc;
    desired elements out of data range are skipped (in prior findindices function)
    # Saving of gaussian fits of peaks could be stored as separate csv if this was ever desired... probably not
    ''' 
    SEMfileName=logmatch.Filename # only used for error reporting 
    # Create temp df to hold and pass linear fit data    
    mycols=['Basename', 'Filenumber', 'Point','Filename', 'Filepath', 'Sample', 'Comments', 'Element',
    'Xc', 'Width', 'Peakarea', 'Y0','Rsquared'] # for gaussian peak fits 
    mycols2=['Basename','Filenumber', 'Point', 'Filename', 'Filepath', 'Sample', 'Comments', 'Element', 'Energy', 'Shift', 'Rawcounts', 
    'Backcounts', 'Subtractedcounts', 'Adj counts', '% err', 'Significance' , 'Basis', 'Errbasis',]  # for integration results 
    Peakfits=pd.DataFrame(columns=mycols) # blank df for this spectrum's peak fits
    Integresults=pd.DataFrame(columns=mycols2) # blank df for this spectrum's integration results
    
    # fit all elemental peaks with gaussian, determine shift and perform integration (incl. error)
    for i, (elem, idealindex, halfwidth, kfact, errkfact, mass) in enumerate(Elemdata):
        Peakfitrow=pd.DataFrame(index=np.arange(0,1),columns=mycols) # single dataframe row for this
        Integresultrow=pd.DataFrame(index=np.arange(0,1),columns=mycols2) # blank df row
        # linear fit below this elem's peak (shifts and adjustments already made)
        # use 10 more channels than those used for integration for gaussian fits
        fitregion=SEMfile[idealindex-halfwidth-5:idealindex+halfwidth+6]
        if fitregion.empty==True: # skip if no data present (peak out of range problem)
            continue
        # Gaussian fit of subtracted data peaks > 50 cnts
        if fitregion['Subdata'].max()>50:
            fitregion, fitparams, rsquared, ier = fitgauss(fitregion, halfwidth, elem, SEMfileName, savegauss=True)
            if savegauss==True: # save Gaussian peaks as separate column
                if 'Gauss' not in SEMfile.dtypes.index: # add col if not already present                
                    SEMfile['Gauss']='' # add blank col for gaussian fit if not present
                # copy gaussian fit to Augerfile... fitregion only modified in new Gauss peak fit column
                SEMfile.loc[fitregion.index,fitregion.columns]=fitregion
                      # determination of peak shift 
        # If gaussian fit is successful set center integration channel to index nearest xc
        # ier flag of 1,2,3,4 if fit succeeds but rsquared threshold is better
            if rsquared!='n/a': # somewhat successful gaussian fit 
                if rsquared>0.4:
                    xc=fitparams[0] # center of gaussian fit in keV
                    centerindex=int((xc+.01)*100)
                    shift= centerindex- idealindex # energy shift in channels
                    if abs(shift)>3:
                        print('Warning: Gaussian shift of ', str(shift), ' channels indicated for ', elem, ' in ', SEMfileName)
                        if shift>0: # keep peak shift the same but only allow 3 channel shift in integration window
                            centerindex=idealindex+3 # set to max shift
                        else:
                            centerindex=idealindex-3
            # TODO Maybe a better way of setting maximum allowable shift        
                else: 
                    print('Low quality gaussian fit for ', elem, ' in ', SEMfileName)
                    centerindex=idealindex # already stores index number of central peak (ideal - sm-diff shift value)
                    shift='n/a'
            # Write gaussian fit params to peakfit (eventually copied to peakfitlog)
                    
            else: # Fit attempted but failed result
                print ('Fit attempted but result failed for ', elem, ' in ', SEMfileName)
                fitparams=['n/a','n/a','n/a','n/a']            
                rsquared='n/a'
                
        else: # indication of failed Gaussian fit (use prior knowledge of peak position)
            print('Skip gaussian fit of tiny ', elem, ' peak in ', SEMfileName)
            # set center integration channel to value passed by integpeak 
            # this is ideal energy value but adjusted by shift found using smooth-diff quant method
            centerindex=idealindex # already stores index number of central peak (ideal - sm-diff shift value)
            shift='n/a'
            fitparams=['n/a','n/a','n/a','n/a']            
            rsquared='n/a'
        # Perform integration over peak center channel + integwidth on either side 
        SEMpeak=SEMfile[centerindex-halfwidth:centerindex+halfwidth+1]
        integcounts=SEMpeak['Subdata'].sum() # get counts sum 
        backgroundcnts=SEMpeak['Backfit'].sum() # sum counts over identical width in background fit
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
        Integresultrow.iloc[0]['Rawcounts']=SEMpeak['Counts'].sum() 
        Integresultrow.iloc[0]['Backcounts']=backgroundcnts
        Integresultrow.iloc[0]['Subtractedcounts']=integcounts
        # Adjusted counts must be determined later for pathological overlaps
        # 2 sigma err due to counting statistics
        Integresultrow.iloc[0]['% err']=round(2/np.sqrt(integcounts),3)
        Integresultrow.iloc[0]['Significance']=round(integcounts/(np.sqrt(backgroundcnts)),3)
		# TODO add 2/sqrt(n) calc of associated percent error (also can calculate later)
        Integresultrow.iloc[0]['Basis']=integcounts*kfact/mass
        # Calculated combined error for 2sig counting stats + loaded k-factor error
        comberr=np.sqrt(errkfact**2+(2/np.sqrt(integcounts))**2)
        # calculate error in basis for given elemental peak
        Integresultrow.iloc[0]['Errbasis']=(integcounts*kfact/mass)*comberr
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
    Peakfits=Peakfits[mycols] # put back in original order
    Integresults=Integresults[mycols2] # put back in original order
    return SEMfile, Peakfits, Integresults # df with direct peak fitting info for all areas/ all elements
    
def fitbackgrounds(SEMfile, Fitregions, logmatch, oddspectrum=False):
    ''' Background fit for each direct peak(opens source spectrum as SEMfile, 
    fits peak backgrounds above and below using Elemdata, saves background to source csv (overwrites existing fits), 
    also saves linear fit params to logdataframe with position/amplitude/etc;
    Fitregions stores total ev range, background regions, fit type and thresholdfor deriv knockout '''
    # Create temp df to hold and pass linear fit data
    SEMfileName=logmatch.Filename # 
    mycols=['Basename', 'Filenumber', 'Filename', 'FilePath', 'Sample', 'Comments', 'Date', 'Point', 'Beamkv', 
    'Livetime','Timeconst','Deadfraction','Fitrange', 'Fittype', 'A', 'B', 'C', 'D', 'Rval', 'Pval', 'Stderr']
    Backfitparams=pd.DataFrame(columns=mycols) # empty df to hold all rows from this spectrum
        # all fit regions modify fit region boundaries for this spectrum based on smooth-differentiated peak (2nd deriv, Savgol (poly=2, pts=11))
        # global shifts from smdifpeaks and local shift based on smoothed 2nd derivative 
        # already incorporated into Elemdata values (lower1,2 and upper1,2 fully adjusted)
    
    # loop through and fit all peaks for each element in this spatial area            
    for i, [fitrange, fitregs, fittype, threshold] in enumerate(Fitregions):
        # create new df row for each fitted range
        Backfitparamrow=pd.DataFrame(index=np.arange(0,1),columns=mycols)
        # modify fit region for this spectrum (eliminate those with high derivative )
        # Threshold level defined by many attempted fits to actual data
        Thisbackfit=findfitregion(SEMfile, fitregs, threshold)
        # Force counts to zero near origin
        for index,row in Thisbackfit.iterrows():
            if index < 5:
                Thisbackfit=Thisbackfit.set_value(index,'Counts',0)
        # now do parabolic fit over this region (return df with backfit col)
        if fittype=='parabola':
            Thisbackfit, fitparams = fitparabola(Thisbackfit, SEMfileName)
            # unpack polynomial fit parameters 
            A=fitparams[0]
            B=fitparams[1]
            C=fitparams[2]    
            # now copy this function over entire range of fit
            lower=int(fitrange.split('-')[0])
            upper=int(fitrange.split('-')[1])
            # TODO these are essentially index numbers (since energy in eV is nearly same as index range)    
            if A!='n/a': # test for successful fit (all fitparams set to n/a)
                for i in range(lower,upper): 
                    xval=SEMfile.iloc[i]['Energy']
                    SEMfile.set_value(i,'Backfit',A * xval**2 + B * xval + C) # just set values directly from fit results
            # now store values forthis df row (slower but easy)
            Backfitparamrow.iloc[0]['Fitrange'] = fitrange
            Backfitparamrow.iloc[0]['Fittype'] = fittype
            Backfitparamrow.iloc[0]['A'] = A # parabolic fit params
            Backfitparamrow.iloc[0]['B'] = B
            Backfitparamrow.iloc[0]['C'] = C
        # TODO test and incorporate other fit types
        # Set value for subtracted spectral data (nan if nothing in backfit)
        SEMfile['Subdata']=SEMfile['Counts']-SEMfile['Backfit']
        # concatentate single row with log 
        Backfitparams=pd.concat([Backfitparams, Backfitparamrow], ignore_index=True)
        # END OF BACKGROUND FITTING LOOP FOR SINGLE SPECTRUM
    # create subtracted peak for entire spectrum
    
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
        
    Backfitparams=Backfitparams[mycols] # put back in original order
    return SEMfile, Backfitparams # df with direct peak fitting info for all areas/ all elements

        
def batchSEMquant(SEMfiles, Fitregionsdf, SEMquantparams, Elements, overwrite=True, savegauss=True):
    ''' Batch quantification of all peaks in Elements list 
    returns df with peak positions, amplitudes, width, energy shift, etc. '''   
    # create empty dataframe for storing/passing linear fit params (same structure as in fitbackgrounds)
    mycols=['Basename', 'Filenumber', 'Filename', 'FilePath', 'Sample', 'Comments', 'Date', 'Point', 'Beamkv', 
    'Livetime','Timeconst','Deadfraction','Fitrange', 'Fittype', 'A', 'B', 'C', 'D', 'Rval', 'Pval', 'Stderr']    
    Backfitparamslog=pd.DataFrame(columns=mycols) # empty df to hold all fits/ all spectra 
    # black df for gaussian peak fit results
    mycols2=['Basename', 'Filenumber', 'Point','Filename', 'Filepath', 'Sample', 'Comments', 'Element',
    'Xc', 'Width', 'Peakarea', 'Y0','Rsquared'] # for gaussian peak fits 
    Peakfitlog=pd.DataFrame(columns=mycols2)
    # Now a blank frame for integrated quant results
    mycols3=['Basename','Filenumber', 'Point', 'Filename', 'Filepath', 'Sample', 'Comments', 'Element', 'Energy', 'Shift', 'Rawcounts', 
    'Backcounts', 'Subtractedcounts', 'Adj counts', '% err', 'Significance' , 'Basis', 'Errbasis',] # for integration results     
    Integquantlog=pd.DataFrame(columns=mycols3) 
    # Get global values for Elemdata for standard spectra
    Elemdata=findelemregions(Elements, SEMquantparams)
    
    for i in range(0,len(SEMfiles)):
        # get ith row from parameters log for subset of selected spe files (i.e. from SEMfiles)
        logmatch=SEMfiles.iloc[i] #contains row with filename and all other parameters from a given spectra 
        logmatch=logmatch.squeeze() # convert/flatten to Series
        # load Auger spe file of interest here
        SEMfileName=logmatch.Filename  # get Auger filename from Series
        SEMfile=openSEM(SEMfileName) # read csv or raw psmsa
            
        # Check for necessary columns in SEMfile
        if 'Counts' not in SEMfile:
            print('Counts not present in file ', SEMfileName)
            continue # skip to next spectrum (this shouldn't happen)
        for index, row in SEMfile.iterrows():
            if index <= 5: # 
                SEMfile=SEMfile.set_value(index,'Counts',0)
        if 'Backfit' not in SEMfile: # add this background fit column if not presentfs
            SEMfile['Backfit']=np.nan
        if overwrite==True: # clear all prior background, subdata and peak fitting regions
            SEMfile['Backfit']=np.nan
            SEMfile['Subdata']=np.nan
            SEMfile['Gauss']=np.nan
        # Sav-gol 2nd deriv column used to guide selection of fitting regions
        if 'Savgol' not in SEMfile: # returns df with this Savgol column added
            SEMfile=makesavgol(SEMfile)  # FUNCT pass full spectrum for given area (saved below)   
        if 'Subdata' not in SEMfile: # add col for subtracted peak data
            SEMfile['Subdata']=np.nan
        # Unpack fitregions from dataframe into list of data ev ranges, 
        # For each ev range we have total fit range, regions which it holding background and fit type
        Fitregions=unpackfitregs(Fitregionsdf) # turn dataframe into list of lists
        # Linear background fits over 4 or 5 energy regions (while dodging major peaks)
        SEMfile, Backfitparams, = fitbackgrounds(SEMfile, Fitregions, logmatch)
        # TODO handling of odd spectra using Elemdatamod (looks up positions in this SEMfile)
        # use global Elemdata, fit major peaks in Subdata with gaussians (and determine shift)
        # Fits peaks, determines shifts and does integration and basis calculations for each element
        SEMfile, Peakfits, Integresults = fitpeaks(SEMfile, Elemdata, logmatch, savegauss=True)

        # concat result rows from this spectrum to longer master log
        Backfitparamslog=pd.concat([Backfitparamslog,Backfitparams],ignore_index=True)
        Peakfitlog=pd.concat([Peakfitlog,Peakfits],ignore_index=True)
        Integquantlog=pd.concat([Integquantlog,Integresults],ignore_index=True)
        
        # direct save of modified auger csv with new linear background fits (after all areas processed)
        if '.psmsa' in SEMfileName:
            SEMfileName=SEMfileName.split('.psmsa')[0]+'.csv' # save as .csv not .psmsa
        SEMfile.to_csv(SEMfileName, index=False)
    # now sort and reorder master log files for backfits, peakfits, integration of all peaks all spectra
    Backfitparamslog=Backfitparamslog.sort_values(['Basename','Filenumber'], ascending=True)
    Backfitparamslog=Backfitparamslog[mycols] # put back in original order
    Peakfitlog=Peakfitlog.sort_values(['Basename','Filenumber'], ascending=True)
    Peakfitlog=Peakfitlog[mycols2] # put back in original order    
    Integquantlog=Integquantlog.sort_values(['Basename','Filenumber'], ascending=True)
    Integquantlog=Integquantlog[mycols3] # put back in original order    
    return Backfitparamslog, Peakfitlog, Integquantlog
	
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
    
def getfromlog(SEMparamrow, SEMlogbook):
    ''' Find associated entry in SEM Excel log book  '''
    # need to match base name & spectral number (currently assumes same data entry for all pts within single point & shoot)
    basename=SEMparamrow.iloc[0]['Basename']
    specnum=SEMparamrow.iloc[0]['Filenumber']
    ptnum=SEMparamrow.iloc[0]['Point']
    mask=SEMlogbook['Basename'].str.contains(basename, case=False)
    match=SEMlogbook.loc[mask]    
    match=match[(match['Filenumber']==specnum)&(match['Point']==ptnum)]
    if len(match)!=1:
        print('Problem finding unique log entry for ', basename, specnum)
    SEMparamrow=SEMparamrow.set_value(0, 'Project', match.iloc[0]['Project']) # copy project name
    SEMparamrow=SEMparamrow.set_value(0, 'Sample', match.iloc[0]['Sample']) # copy sample name
    SEMparamrow=SEMparamrow.set_value(0, 'Comments', match.iloc[0]['Comments']) # copy comments
    return SEMparamrow
    
def stripparams(header, df, filename):
    '''Generalized parameter finding from header info:  pass string header and template DF with column names, start and end strings and data type  '''
    # use SEMparam.csv as template for data rows
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
		

def getfromname(SEMparamrow,filename):
    ''' Extract base name, spectral number and pt # from filename and current data path 
    write all to param log '''
    pattern=re.compile(r'(\(\d+\))') # find spectral number inside parens
    match=re.search(pattern, filename)
    path=os.getcwd() # current data directory 
    try:
        basename=filename[0:match.start()] # pulls out base name
        specnum=int(filename[match.start()+1:match.end()-1]) # gets spectral number within parentheses
        ptnum=int(filename.split('_pt')[1].split('.')[0]) # gets pt number from 'name(10)_pt1.emsa'
        SEMparamrow=SEMparamrow.set_value(0,'Basename',basename)
        SEMparamrow=SEMparamrow.set_value(0,'Filenumber',specnum)
        SEMparamrow=SEMparamrow.set_value(0,'Filename',filename) # use definite filename instead of title from header
        SEMparamrow=SEMparamrow.set_value(0,'Point',ptnum)
        SEMparamrow=SEMparamrow.set_value(0,'FilePath',path) 
    except:
        print('Problem extracting names from logfile for ', filename)
    return SEMparamrow


def getparams(filelist, paramtemplate, SEMlogbook):
    '''Main loop that opens files, finds associated parameters from header and combines with logbook and sample name info '''
    # TODO check for consistency between dataframes
    mycols=['Project','Basename','Filenumber','Filename','FilePath','Sample','Point','Comments','Date','Time','Beamkv','Livetime','Realtime','Detected','Converted','Stored','Timeconst','Deadfraction']
    SEMparamlog=pd.DataFrame(columns=mycols)
    
    for i,filename in enumerate(filelist):
        with open(filename, 'r') as file:
            filedata = file.read()
        header=filedata.split('#SPECTRUM    :')[0]
        # all params extracted and placed in single dataframe row
        SEMparamrow=stripparams(header, paramtemplate, filename)
        SEMparamrow=getfromname(SEMparamrow,filename)  # get basename, spectral number and pt # from name
        
        # now find matching sample info, project names from Excel logbook
        SEMparamrow=getfromlog(SEMparamrow, SEMlogbook)
        
        SEMparamlog=pd.concat([SEMparamlog,SEMparamrow], ignore_index=True)
    SEMparamlog=SEMparamlog.sort_values(['Basename','Filenumber','Point'], ascending=True) # sort by basename then filenumber
    SEMparamlog=SEMparamlog[mycols] # put back in original order
    SEMparamlog=SEMparamlog.reset_index(drop=True) # reset the index
    return SEMparamlog

def checklogfile(filelist, SEMlogbook):
    ''' Checks the user Auger logbook Excel for consistency with the actual data file list from directory
    prints out filenumbers that have a problem to console''' 
    # TODO modify this to allow for different base names
    psmsa=[] # list of filenumbers of psmsa files in directory
    pattern=re.compile(r'(\(\d+\))') # find spectral number inside parens

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
    
def openorcreatelogbook(filelist):
    ''' Looks for existing csv or xls log file ... if not found makes new one by calling makeblanklog ''' 
    logfile=glob.glob('*SEMEDX_log*') # find ... Auger_logbook.
    if len(logfile)==1: # found logbook
        name=logfile[0]
        if '.xls' in name: # open log tab of existing excel file
            SEMlogbook=pd.read_excel(name, sheetname='Log')        
        if '.csv' in name: # open csv
            SEMlogbook=pd.read_csv(name)
    elif len(logfile)==0:
        SEMlogbook=makeblanklog(filelist)
    else:
        print('Error: There must be only one SEMEDX log in folder.')
        return 
    return SEMlogbook
            
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
        basename=filename[0:match.start()] # pulls out base name
        specnum=int(filename[match.start()+1:match.end()-1]) # gets spectral number within parentheses
        ptnum=int(filename.split('_pt')[1].split('.')[0]) # gets pt number from 'name(10)_pt1.emsa'
        Samplelogrow=Samplelogrow.set_value(0,'Basename', basename)
        Samplelogrow=Samplelogrow.set_value(0,'Filenumber', specnum)
        Samplelogrow=Samplelogrow.set_value(0,'Point', ptnum)
        Samplelogrow=Samplelogrow.set_value(0,'Filename', filename)
        Samplelogrow=Samplelogrow.set_value(0,'Project', projname)
        Samplelogrow=Samplelogrow.set_value(0,'FilePath', fullpath)
        SEMlogbook=pd.concat([SEMlogbook,Samplelogrow], ignore_index=True)
    csvname=projname+'SEMEDX_logbook.csv'
    SEMlogbook.sort_values(['Basename', 'Filenumber'])
    SEMlogbook=SEMlogbook[mycols] # reorder columns to standard
    SEMlogbook.to_csv(csvname, index=False)
    print('Blank logbook created for project ', projname, '; Sample names and comments can be manually entered .')
    return SEMlogbook

#%% 


#%% OLD STUFF Function to compute smoothed-deriv of counts column 
# http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
'''
def smooth_diff(emsafile):
    import numpy as np
    from math import factorial
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except:
        ValueError, msg:
    raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
    raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
    raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
'''