# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:30:11 2016

@author: tkc
"""
import re, sys
import pandas as pd
import numpy as np
from collections import defaultdict
import tkinter as tk
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Utilities' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Utilities')
#%% 

# Testing     row=EDXsumm.loc[index]


def compsummary(EDXcomp, Elements, Elemexcl):
    ''' Compositional summary that keeps at % and identifying fields only '''
    mycols=['Filename', 'Sample', 'Comments', 'SEMbasis', 'Phase']
    # Handle quant including the excluded elements
    missing=[i for i in Elemexcl if i not in Elements]
    if len(missing)>0:
        print(','.join(missing),' excluded elements missing from Element list')
    missing=[i for i in Elements if i not in EDXcomp.columns]
    if len(missing)>0:
        print(','.join(missing),' elements not present in EDXcomp... removed')
        for i, elem in enumerate(missing):
            Elements.remove(elem)
    real=[i for i in Elements if i not in Elemexcl]
    # order with real elems first and excluded second
    for i, elem in enumerate(real):
        mycols.append('%'+elem)
    mycols.append('Total')
    if 'Total' not in EDXcomp.columns:
        EDXcomp['Total']=np.nan
    for i, elem in enumerate(Elemexcl):
        mycols.append('%'+elem)
    for index, row in EDXcomp.iterrows():
        elemsumm=0.0
        # Compute at.% including all elems
        for i, elem in enumerate(Elements):
            elemsumm+=row['%'+elem]
        for i, elem in enumerate(Elements):
            EDXcomp=EDXcomp.set_value(index,'%'+elem, 
                100*EDXcomp.loc[index]['%'+elem]/elemsumm)
        # Redo the sum only for included (real) elems
        elemsumm=0.0
        for i, elem in enumerate(real):
            elemsumm+=row['%'+elem]
        for i, elem in enumerate(real):
            EDXcomp=EDXcomp.set_value(index,'%'+elem, 
                    100*EDXcomp.loc[index]['%'+elem]/elemsumm)
        EDXcomp=EDXcomp.set_value(index,'Total', 
                    elemsumm)
    EDXsumm=EDXcomp[mycols]
    return EDXsumm

def calcadjcnts(Integquantlog, elem1, elem2, adjcntratio, adjerr, kfacts):
    '''Performs count, corrected count and error adjustment using pair of elements; elem1 is quant element in question
    elem2 is other interfering peak, adjcntratio and adjerr from SEM_interferences 
    kfacts is list of params for this elem from SEMquantparams ''' 
    filelist=np.ndarray.tolist(Integquantlog.Filename.unique())
    kfactor, errkfact, mass=kfacts # unpack kfactor main params
    for i, fname in enumerate(filelist):
        match1=Integquantlog[(Integquantlog['Filename']==fname) & (Integquantlog['Element']==elem1)]
        match2=Integquantlog[(Integquantlog['Filename']==fname) & (Integquantlog['Element']==elem2)]
        if len(match2)!=1 or len(match1)!=1:
            print('Problem finding ', elem1,' and/or', elem2, 'for ', fname)
            continue
        elem2cnts=match2.iloc[0]['Subtractedcounts']
        elem1cnts=match1.iloc[0]['Subtractedcounts']
        err1=match1.iloc[0]['% err'] # fractional/relative error for element in question
        err2=match2.iloc[0]['% err'] # fractional/relative error for interfering element
        term2fracterr=np.sqrt(err2**2+adjerr**2) # fractional error in term 2
        term2err=term2fracterr*elem2cnts*adjcntratio # absolute error in term2 correction
        adjcnts=elem1cnts-elem2cnts*adjcntratio # adjusted counts ratio for element 1
        newerr=np.sqrt(err1**2+term2err**2) # absolute error in adjusted counts
        newadjerr=newerr/adjcnts # New fractional error in elem 1
        match1=match1.set_value(match1.index[0], 'Adjcounts', max(adjcnts,0))
        if adjcnts>0:
            match1=match1.set_value(match1.index[0], '% err', newadjerr) # reset error to include that due to interference
        # now recalculate corrected counts and associated error
        backcnts=match1.iloc[0]['Backcounts']
        newcorrcnts=adjcnts*kfactor/mass # standard value
        if adjcnts<2*np.sqrt(backcnts): # 2sigma of background as lower limit
            newcorrcnts=2*np.sqrt(backcnts)*kfactor/mass # set to lower limit
            print ('2sigma of background limiting value used for ', elem1, fname)
            match1=match1.set_value(match1.index[0],'Significance',0) # set to zero as marker of limiting value
        match1=match1.set_value(match1.index[0],'Correctedcounts',newcorrcnts)
        # find combined 2 sigma error (counts and k-factor) as percent error
        comberr=np.sqrt(errkfact**2+newadjerr**2) # combine the fractional errors 
        match1=match1.set_value(match1.index[0],'Errcorrcnts',newcorrcnts*comberr)
        Integquantlog.loc[match1.index,match1.columns]=match1 # copies altered row back to main log
    return Integquantlog
        
def recalcadjbatch(df, Interferences, SEMquantparams):
    '''Calculate adjusted counts (from unresolvable interferences), then recalculate corrected counts with updated k-factor/ k-fact error 
    and mass result stored in corrcnts column and used for subsequent compositional determinations
    can change AESquantresults and recalc at any time ''' 
    if 'Correctedcounts' not in df:
        df['Correctedcounts']=0.0 # new column for adjusted amplitude Correctedcounts (if not already present)
    if 'Errcorrcnts' not in df:
        df['Errcorrcnts']=0.0 # new column for error
    # loop for each element, mask df, get appropriate k-factor & mass
    df=df.reset_index(drop=True) # go ahead and reset index
    elemlist=np.ndarray.tolist(df.Element.unique()) # list of unique elements from df
    for index, row in Interferences.iterrows():
        elem1=Interferences.loc[index]['Elem1']
        elem2=Interferences.loc[index]['Elem2']
        if elem1 in elemlist and elem2 in elemlist: # do the subtraction
            adjratio=Interferences.loc[index]['Adjcntratio']
            adjerr=Interferences.loc[index]['Percenterr'] # error as percentage for this adjustment
            matchelem=SEMquantparams[SEMquantparams['element']==elem1]
            if not matchelem.empty:
                kfacts=[matchelem.iloc[0]['kfactor'],matchelem.iloc[0]['errkfact'],matchelem.iloc[0]['mass']] # kfactor and mass for this element/peak
                df= calcadjcnts(df, elem1, elem2, adjratio, adjerr, kfacts) # makes correction for all in list, incl new corrcnts and error estimates
    return df

def assembledataset(paramloglist, integloglist):
    '''Construct master paramlog, integlog, backfitlog and peakslog for a list of directories '''
    mycols=['Project','Basename','Filenumber','Filename','FilePath','Sample','Point','Comments','Date','Time','Beamkv','Livetime','Realtime','Detected','Converted','Stored','Timeconst','Deadfraction']
    Masterparamlog=pd.DataFrame(columns=mycols)
    mycols2=['Basename','Filenumber', 'Point', 'Filename', 'Filepath', 'Sample', 'Comments', 'Element', 'Energy', 'Shift', 'Rawcounts', 
    'Backcounts', 'Subtractedcounts', 'Adjcounts', '% err', 'Significance' , 'Correctedcounts', 'Errcorrcnts',] # for integration results     
    Masterinteglog=pd.DataFrame(columns=mycols2) # empty frame   
    for i, logfile in enumerate(paramloglist):
        thisparam=pd.read_csv(logfile)
        Masterparamlog=pd.concat([Masterparamlog,thisparam], ignore_index=True)
    for i, logfile in enumerate(integloglist):
        thisinteg=pd.read_csv(logfile)
        Masterinteglog=pd.concat([Masterinteglog,thisinteg], ignore_index=True)
    return Masterparamlog, Masterinteglog
	
def organizecolumns(df1,mycols):
    ''' Pass df and template (list of desired columns in desired order) and return reorganized newdf
    '''
    cols1=df1.columns.tolist()
    newdf=df1 # avoids modification of passed df
    uniquelist=[i for i in cols1 if i not in mycols]
    for i,colname in enumerate(uniquelist): # remove cols from df1 that are absent from df2
        # newdf.drop(colname, axis=1, inplace=True) # this modifies both passed and returned dfs
        newdf=newdf.drop(colname, axis=1)
    newdf=newdf[mycols] # reorder columns based on template df
    return newdf
    
def thresholdSEM(thisval,thisthresh,threshold):
    '''Pass value, threshold smdifpeak info on single peak; return df if ratio above threshold, return empty frame if below '''
    if threshold=='Lq':
        keeplist=['Lq']
    if threshold=='Ld':
        keeplist=['Lq','Ld']
    if threshold=='Lc':
        keeplist=['Lq','Ld','Lc']
    if thisthresh in keeplist:
        value=thisval
    else:
        value=0
    return value # return empty df if amplitude below threshold

def outputduplicates(df):
    '''Prints out names of samples with duplicated entries into console'''
    tempdf=df.duplicated(['Sample']) # series with 2nd of duplicated entries as True
    for i in range(0,len(tempdf)):
        if tempdf[i]==True:
            sample=df.iloc[i]['Sample']
            print('Duplicated spectra for sample: ', sample)
    return

def parseelemlist(elemlist):
    '''Find and separate multielement peaks to be averaged (e.g. Fe2 & Fe) from longer string of element peaks
    e.g. splits "Mg Fe Fe2 Si" into "Mg Si" and "{Fe,[Fe,Fe2]} dictionary'''
    # Strip numbers from strings within list  
    newlist=[re.match('\D+',i).group(0) for i in elemlist]
    
    # find duplicated peaks (multiple peaks per element)
    Multielem = defaultdict(list)
    for i, item in enumerate(newlist):
        Multielem[item].append(i)
    Multielem = {k:v for k,v in Multielem.items() if len(v)>1} # dictionary with duplicated item and list with indices
    
    duplist=list(Multielem.values()) # get list  
    duplist=[item for sublist in duplist for item in sublist] # single list with positions of duplicated elements
    
    # now alter multipeak elements list to give dict with element and then list of peak for that element    
    for key,value in Multielem.items():
        templist=value # dictionary value is list of elem peak index positions
        peaklist=[]
        for i, index in enumerate(templist): # create new list with original elem peak from index positions
            peaklist.append(elemlist[index])
        # now replace list of index positions with elempeak names
        Multielem.update({key:peaklist}) # key will be multipeak element string i.e. "Fe"
    # finally construct new single elements list with multipeak ones removed (handle each separately)
    newelemlist=[]
    for i in range(0,len(elemlist)):
        if i not in duplist:
            newelemlist.append(elemlist[i])
    return newelemlist, Multielem
    
    
def parseelem2(elemlist, Multielem):
    ''' After multielement peaks removed, also move secondary peaks used as primary to dict (handle separately)
    e.g. splits "S Mg Fe2 Si" into "S Mg Si" and "{Fe,[Fe2]} dictionary; same structure and df output 
    for averaging of Fe, Fe2, or straight Fe2 or straight Fe'''
    
    newelemlist=[]
    for i, elem in enumerate(elemlist):
        if re.search(r'\d',elem): # has number
            templist=[] # peakIDs added as list (of length 1)
            templist.append(elem) # list containing single string (keeps identical data structure)
            newkey=re.search(r'\D+','Fe2').group(0)
            Multielem.update({newkey:templist}) # add to existing dictionary for separate handling
        else:
            newelemlist.append(elemlist[i]) # just copy over 
    return newelemlist, Multielem # return altered element list and multielem dictionary    

def getcountscomp(df, Integquantlog, elemlist, sigthreshold=2):
    '''Retrieves counts for elements listed in elemlist; uses subset of files in df and finds associated elemental quant
    data from accompanying integquantlog (usually fuller set);  returns counts, corrected counts and at % 
    normal usage though is for looking at count ratios (i.e. determining pathological overlaps)
    '''
    # thresholds=getelemthresholds(elemlist, SEMquantparams) # Get list of sigma levels for significance/inclusion 
    # thresholds for both single and multipeak
    # two element lists needed (elements with one peak and elements with compositions averaged from two peaks i.e. Fe2, Fe3)
    df=df.reset_index(drop=True)
    df['SEMbasis']=0.0 # resets to zero if already present from calcamplitude
    mycols=['Filenumber', 'Project', 'Basename', 'Point', 'Filename', 'FilePath', 'Sample', 'Comments','SEMbasis','Beamkv','Timeconst']
    for i, elem in enumerate(elemlist):  # add columns for basis
        cntsname=elem+'cnts' # at % columns named %S, %Mg, etc.        
        df[cntsname]=0.0
        df[elem]=0.0 # add col for each element to spelist
        mycols.append(cntsname)
        mycols.append(elem)
    for i, elem in enumerate(elemlist):  # now add at.% columns (e.g. %S, %Mg)
        colname='%'+elem # at % columns named %S, %Mg, etc.
        errname='err%'+elem        
        mycols.append(colname)  # add to column list template
        mycols.append(errname)
        df[colname]=0.0
        df[errname]=0.0
    for i in range(0,len(df)): # loop through each spectrum (match basename, filenumber and point number for psmsa)
        basename=df.iloc[i]['Basename']
        filenum=df.iloc[i]['Filenumber']
        ptnum=df.iloc[i]['Point']
        match=Integquantlog[(Integquantlog['Filenumber']==filenum) & (Integquantlog['Point']==ptnum) & (Integquantlog['Basename']==basename)] # find integ data for this filenumber
        # match=match[match['Area']==1] (select area 1)
        basis=0.0 #
        for j, elem in enumerate(elemlist): # handle the single peak elements            
            temp=match[match['Element']==elem] # finds entry for this element 
            if len(temp)==1:
                # thisthresh=thresholds.get(elem) # sig level for this element
                if temp.iloc[0]['Significance']>sigthreshold: # if above set threshold then calculate elem's value and add to basis
                    df=df.set_value(i, elem, temp.iloc[0]['Correctedcounts']) # copy adjusted counts of this element
                    cntsname=elem+'cnts'                    
                    df=df.set_value(i, cntsname, temp.iloc[0]['Subtractedcounts']) # 
                    basis+=temp.iloc[0]['Correctedcounts'] # add this element's value to AES basis
        df=df.set_value(i,'SEMbasis',basis) # basis is sum of all chosen corrcnts 
        # Now compute at.% for each listed element (incl errors)
        for j, elem in enumerate(elemlist):
            colname='%'+elem
            ratio=df.iloc[i][elem]/df.iloc[i]['SEMbasis'] # initialized to zero in cases where peak is below significance threshold
            df.set_value(i, colname, 100*ratio) # this ratio is the atomic percentage 
            temp=match[match['Element']==elem] # again find peak entry and get finds entry for this peak
            # TODO maybe check threshold again (although element's value will be zero)
            if len(temp)==1: 
                thiserr=temp.iloc[0]['Errcorrcnts'] # already includes countings stats and k-factor error %
                # this error calculated during SEMimport is not a percentage (fractional error) but is an absolute error 
                atpercerr=thiserr/df.iloc[i]['SEMbasis'] # fractional error calculation (use entire basis as demon
                errname='err%'+elem # error column
                abserr=atpercerr*ratio # gets absolute error in atomic percent calculation 
                df.set_value(i, errname, 100*abserr) # Writes fractional error in at%
                # TODO don't we want absolute error here
        # Also calculate for elements w/ multiple peaks (if present)              
    # re-organize data based on mycols template
    df=organizecolumns(df,mycols) # 
    return df

def recalccorrcounts(df, SEMquantparams):
    '''For each elemental peak in integquantlog, recalculate corrected counts with updated k-factor/ k-fact error and mass
    result stored in corrcnts column and used for subsequent compositional determinations
    can change AESquantresults and recalc at any time ''' 
    if 'Correctedcounts' not in df:
        df['Correctedcounts']=0.0 # new column for adjusted amplitude Correctedcounts (if not already present)
    if 'Errcorrcnts' not in df:
        df['Errcorrcnts']=0.0 # new column for error
    # loop for each element, mask df, get appropriate k-factor & mass
    df=df.reset_index(drop=True) # go ahead and reset index
    elemlist=np.ndarray.tolist(df.Element.unique()) # list of unique elements from df
    for i,elem in enumerate(elemlist):
        match=SEMquantparams[(SEMquantparams['element']==elem)]
        if not match.empty: # skip if this elem is not found
            match=match.reset_index(drop=True)
            kfactor=match.iloc[0]['kfactor'] # kfactor and mass for this element/peak
            mass=match.iloc[0]['mass'] 
            errkfact=match.iloc[0]['errkfact'] # % error in this element's k-factor (should be 2 sigma)
            elemmask=(df['Element']==elem) # mask for this element in loop 
            for j in range(0,len(df)): # loop and set adjamplitude to amp*kfact/mass
                if elemmask[j]==True: # row has this element
                    # check if this element has been corrected for pathological overlap/interference (recalcadjbatch)
                    if str(df.iloc[j]['Adjcounts'])!='nan': # 
                        integcnts=df.iloc[j]['Adjcounts'] # use adjusted value for corrcnts calculation
                    else: # no adjustment so just use subtracted data
                        integcnts=df.iloc[j]['Subtractedcounts']
                    backcnts=df.iloc[j]['Backcounts']
                    newcorrcnts=integcnts*kfactor/mass # standard value
                    if integcnts<2*np.sqrt(backcnts): # 2sigma of background as lower limit
                        integcnts=2*np.sqrt(backcnts)                 
                        newcorrcnts=2*np.sqrt(backcnts)*kfactor/mass # set to lower limit
                        fname=df.iloc[j]['Filename']
                        print ('2sigma of background limiting value used for ', elem, fname)
                        df=df.set_value(j,'Significance',0) # set to zero as marker of limiting value
                    df=df.set_value(j,'Correctedcounts',newcorrcnts)
                    # find combined 2 sigma error (counts and k-factor) as percent error
                    comberr=np.sqrt(errkfact**2+(2/np.sqrt(integcnts))**2)
                    df=df.set_value(j,'Errcorrcnts',newcorrcnts*comberr)
    return df

def calccomp(df, Integquantlog, elemlist, sigthreshold=2):
    '''Calculate elemental composition of given files based on input element list (for python imported/processed data) 
    should elements be eliminated if amplitude is less than 2x that of noise background?
    threshold - ratio of element peak to noise peak (0 means no threshold applied
    load element-dependent significance level from SEMquantparams'''
    # thresholds=getelemthresholds(elemlist, SEMquantparams) # Get list of sigma levels for significance/inclusion 
    # thresholds for both single and multipeak
    elemlist, multipeaklist = parseelemlist(elemlist) # list of single peak elements and dict with multipeaks
    # check if any of the single peaks are secondary (i.e. quant on Fe2 not main Fe)
    elemlist, multipeaklist= parseelem2(elemlist, multipeaklist)
    # two element lists needed (elements with one peak and elements with compositions averaged from two peaks i.e. Fe2, Fe3)
    df=df.reset_index(drop=True)
    df['SEMbasis']=0.0 # resets to zero if already present from calcamplitude
    mycols=['Filenumber', 'Project', 'Basename', 'Point', 'Filename', 'FilePath', 'Sample', 'Comments','SEMbasis','Phase']
    for i, col in enumerate((mycols)):
        if col not in df:
            df[col]=''
    for i, elem in enumerate(elemlist):  # add columns for basis
        df[elem]=0.0 # add col for each element to spelist
        mycols.append(elem)
    for i,elem in enumerate(list(multipeaklist.keys())): # get elements (keys) from dict
        df[elem]=0.0
        mycols.append(elem)
    for i, elem in enumerate(elemlist):  # now add at.% columns (e.g. %S, %Mg)
        colname='%'+elem # at % columns named %S, %Mg, etc.
        errname='err%'+elem
        mycols.append(colname)  # add to column list template
        mycols.append(errname)
        df[colname]=0.0
        df[errname]=0.0
    for i,elem in enumerate(list(multipeaklist.keys())): # add multipeak elements
        colname='%'+elem # at % columns named %S, %Mg, etc.
        errname='err%'+elem 
        mycols.append(colname)  # add to column list template
        mycols.append(errname)        
        df[colname]=0.0
        df[errname]=0.0
    for i in range(0,len(df)): # loop through each spectrum (match basename, filenumber and point number for psmsa)
        basename=df.iloc[i]['Basename']
        filenum=df.iloc[i]['Filenumber']
        ptnum=df.iloc[i]['Point']
        match=Integquantlog[(Integquantlog['Filenumber']==filenum) & (Integquantlog['Point']==ptnum) & (Integquantlog['Basename']==basename)] # find integ data for this filenumber
        # match=match[match['Area']==1] (select area 1)
        basis=0.0 #
        for j, elem in enumerate(elemlist): # handle the single peak elements            
            temp=match[match['Element']==elem] # finds entry for this element 
            if len(temp)==1:
                # thisthresh=thresholds.get(elem) # sig level for this element
                if temp.iloc[0]['Significance']>sigthreshold: # if above set threshold then calculate elem's value and add to basis
                    df=df.set_value(i, elem, temp.iloc[0]['Correctedcounts']) # copy adjusted counts of this element
                    basis+=temp.iloc[0]['Correctedcounts'] # add this element's value to AES basis
        # now handle the multipeak elements (get average value from both peaks)
        for key, value in multipeaklist.items(): # key is element (aka colname in df), value is list of peaks in Smdifpeakslog
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # this is number of lines that are average (i.e. 2 for Fe&Fe2)            
            avgval=0.0 # working value for averaged adjamplitude
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                temp=match[match['Element']==peak] # finds integquantlog entry for this peak (match already trimmed to filenum and area)
                if len(temp)==1:
                    if temp.iloc[0]['Significance']>sigthreshold:
                        avgval+=temp.iloc[0]['Correctedcounts']
                    else:
                        numlines=numlines-1 # if peak is zeroed out and not added, this reduces # peaks in average
            if numlines>0: # avoid divbyzero if peak is too small
                avgval=avgval/numlines # this is now average basis for given element
            df=df.set_value(i, key, avgval) # copy adjusted amplitude of this element
            # add value from this element to SEMbasis
            basis+=avgval
        # end of multipeak elements loop
        df=df.set_value(i, 'SEMbasis', basis) # write total basis value to df
        # Skip if basis is zero (no lines for this spectrum in integlog)
        if basis==0.0:
            print('Basis is zero for spectrum', df.iloc[i]['Filename'], '.. skipped.')
        # Now compute at.% for each listed element (incl errors)
        for j, elem in enumerate(elemlist):
            colname='%'+elem
            ratio=df.iloc[i][elem]/df.iloc[i]['SEMbasis'] # initialized to zero in cases where peak is below significance threshold
            df.set_value(i, colname, 100*ratio) # this ratio is the atomic percentage 
            temp=match[match['Element']==elem] # again find peak entry and get finds entry for this peak
            # TODO maybe check threshold again (although element's value will be zero)
            if len(temp)==1: 
                thiserr=temp.iloc[0]['Errcorrcnts'] # already includes countings stats and k-factor error %
                # this error calculated during SEMimport is not a percentage (fractional error) but is an absolute error 
                atpercerr=thiserr/df.iloc[i]['SEMbasis'] # fractional error calculation (use entire basis as demon
                errname='err%'+elem # error column
                abserr=atpercerr*ratio # gets absolute error in atomic percent calculation 
                df.set_value(i, errname, 100*abserr) # Writes fractional error in at%
                # TODO don't we want absolute error here
        # Also calculate for elements w/ multiple peaks (if present)
        for key, value in multipeaklist.items(): 
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # this is number of lines that are average (i.e. 2 for Fe&Fe2)
            colname='%'+key
            ratio=df.iloc[i][key]/df.iloc[i]['SEMbasis']
            df.set_value(i, colname, 100*ratio)
            # TODO need to propagate errors through Fe & Fe2
            errlist=[] # list of errors in % (usually max of two)
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                temp=match[match['Element']==peak] # finds entry for this peak
                if len(temp)==1: # getting more info from integquantlog
                    if temp.iloc[0]['Correctedcounts']>0: # skip negative values
                        err=temp.iloc[0]['Errcorrcnts']/temp.iloc[0]['Correctedcounts'] # fractional error for this line
                        errlist.append(err) # add this to list 
            # combine errors in quadrature
            totalerr=0.0
            for j, err in enumerate(errlist):
                totalerr+=err**2
            totalerr=np.sqrt(totalerr) # combined fractional errors from peaks1, 2 of same element
            # now get  actual error
            thisval=df.iloc[i][key] # this is averaged corrected counts value computed above (possibly zero if below thresholds )
            thiserr=thisval*totalerr # absolute error (in Fe corrected counts) as actual value based on average of multiple peaks
            atpercerr=thiserr/df.iloc[i]['SEMbasis'] # fractional error in at % calculation
            errname='err%'+ key  # error column
            abserr=atpercerr*ratio # convert fractional error to absolute error (of at % calculation)            
            df.set_value(i, errname, 100*abserr) # Writes absolute error in at% 
        # end of loop calculation for each spectrum 
                
    # re-organize data based on mycols template
    df=organizecolumns(df,mycols) # 
    return df

def calccomposition(df, SEMquantparams, elemlist, sigma=2, threshold='Ld'):
    '''Calculate basis and at % from existing SEM cnts file (Mg_cnts and Mg_level are cnts and Lc,Ld,Lq)
    less relevant and used if python import process has not been done
	use threshold var as Lc, Ld, Lq, 
    '''
    newelemlist, multipeaklist = parseelemlist(elemlist) # list of single peak elements and dict with multipeaks
    # check if any of the single peaks are secondary (i.e. quant on Fe2 not main Fe)
    newelemlist, multipeaklist= parseelem2(newelemlist, multipeaklist)
    
    df=df.reset_index(drop=True)
    df['EDXbasis']=0.0 # resets to zero if already present from calcamplitude
    mycols=['Filename', 'Sample', 'Comments','Livetime','Stored','Deadfraction', 'EDXbasis']
    # set up columns for adjusted basis and at % for each chosen element
    for i, elem in enumerate(elemlist):  # add columns for each peak to hold adjusted counts
        df[elem]=0.0 # add col for each element to spelist
        mycols.append(elem)
    for i, elem in enumerate(newelemlist):  # now add at.% columns for single peaks
        colname='%'+elem # at % columns named %S, %Mg, etc.
        errname='err%'+elem
        mycols.append(colname)  # add to column list template
        mycols.append(errname)          
        df[colname]=0.0
        df[errname]=0.0
    # For at.% the column should be %Fe, not %Fe2 so use keys
    for i,elem in enumerate(list(multipeaklist.keys())): # add columns for multipeak elements
        colname='%'+elem # at % columns named %S, %Mg, etc.
        errname='err%'+elem 
        mycols.append(colname)  # add to column list template
        mycols.append(errname) 
        df[colname]=0.0
        df[errname]=0.0
    # First compute adjusted counts for each element (stored in 'elemname') using original elemlist
    for i,elem in enumerate(elemlist):
        match=SEMquantparams[(SEMquantparams['element']==elem)]
        match=match.reset_index(drop=True)
        if len(match)!=1:
            print('Error finding kfactor info for ', elem)
        kfactor=match.iloc[0]['kfactor'] # kfactor and mass for this element/peak
        mass=match.iloc[0]['mass'] 
        errkfact=float(match.iloc[0]['errkfact']) # 1 sig error in underlying kfactor (not yet well established for SEM-EDX)
        colname=elem+'_cnts' # find counts column for this element
        threshname=elem+'_level'
        for j in range(0,len(df)): # loop and set adjamplitude to amp*kfact/mass
            # get counts value for this element
            thisval=df.iloc[j][colname]
            thisthresh=df.iloc[j][threshname]
            val=thresholdSEM(thisval,thisthresh,threshold) # knock out values below desired EDX threshold 
            df=df.set_value(j,elem,val*kfactor/mass) # will be zero if level is below desired threshold
    
    # Now compute EDXbasis based on adjusted count intensities
    # use newelemlist and multipeaklist in case we're using averages of multiple lines (e.g.FeK and FeL)
    for i in range(0,len(df)):    
        basis=0.0
        # colname is now elem (with adjusted counts for this element)
        # if not above threshold level, value in elem column will be zero (thresholdSEM function above)
        # first add normal single peak elements to basis for this spectrum        
        for j, elem in enumerate(newelemlist): 
            basis+=df.iloc[i][elem] # add each to basis (will be zero if below entered threshold)
        # now handle the multipeak elements (use average value from multiple peaks)
        for key, value in multipeaklist.items(): # key is element (aka colname in df), value is list of peaks in Smdifpeakslog
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # this is number of lines that are average (i.e. 2 for Fe&Fe2)
            avgval=0.0 # working value for averaged adjamplitude
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                avgval+=df.iloc[i][peak] # column w/ adjcount amplitude should match this name
            if numlines>0: # avoid divbyzero if peak is too small
                avgval=avgval/numlines # this is now average basis for given element
            basis+=avgval # add averaged value from this (multipeaked) element to EDXbasis
        df=df.set_value(i, 'EDXbasis', basis) # write total basis value to df 

        # now that EDXbasis is known, compute at.% for each listed element (also including error calculation)
        for j, elem in enumerate(newelemlist):
            colname='%'+elem
            # errors in quadrature and including errkfact and 2/sqrt(N) for each peak
            cntcolname=elem+'_cnts'
            # need to grab associated error for each element again
            match=SEMquantparams[(SEMquantparams['element']==elem)]
            match=match.reset_index(drop=True)
            if len(match)!=1:
                print('Error finding kfactor info for ', elem)
            errkfact=float(match.iloc[0]['errkfact']) # 1 sig error in underlying kfactor
            # just use abs(counts)
            thiserr=sigma/np.sqrt(abs(df.iloc[i][cntcolname]))
            totalerr=np.sqrt(thiserr**2+errkfact**2) # Combine kfactor and counting statistical errors in quadrature
            print ('Count stat and total errors for', elem, df.iloc[i]['Sample'],' are ', thiserr, ' and ',totalerr)            
            # TODO for new python EDXquant, fix problem with possible negative counts
            if df.iloc[i]['EDXbasis']>0: # shouldn't happen but avoid divbyzero
                atperc=df.iloc[i][elem]/df.iloc[i]['EDXbasis']
                thiserr=atperc*totalerr # i.e. error in Mg at. % is Mg at. % times total combined error percentage
                # this is absolute error  (not percentage)                
                df=df.set_value(i, colname, atperc)
                errcolname='err%'+elem
                df=df.set_value(i, errcolname, thiserr) # write to at.% error column 
        # also calculate for elements w/ multiple peaks using average (if present); also includes combined error calculation
        for key, value in multipeaklist.items(): # key is element name (Fe), value(s) are peak names (Fe, Fe2)
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # number of lines that are averaged (i.e. 2 for Fe&Fe2)            
            avgval=0.0 # working value for averaged adjamplitude
            errmultipeaks=0 # variable for combined error calculation (e.g. Fe2 & Fe3 & kfactor)
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                avgval+=df.iloc[i][peak] # column w/ adjcount amplitude should match this name
                cntscolname=peak+'_cnts' # underlying counts and peak is Fe2 or comparable
                thiserr=sigma/np.sqrt(abs(df.iloc[i][cntscolname])) # % error for this peak (i.e. Fe2)
                errmultipeaks+=thiserr*2 # start quadrature combination
            if numlines>0: # avoid divbyzero if peak is too small
                avgval=avgval/numlines # this is now average element basis for given element from multiple peaks
            # need k-factor error and combined error from multiple peaks
            # look up error in k-factor for this element
            match=SEMquantparams[(SEMquantparams['element']==key)]
            match=match.reset_index(drop=True)
            if len(match)!=1:
                print('Error finding kfactor info for ', elem)
                
            errkfact=float(match.iloc[0]['errkfact']) # 1 sig error in underlying kfactor
            combinederr=np.sqrt(errmultipeaks+errkfact**2) # sqrt of sum of squares gives combined error percentage
            if df.iloc[i]['EDXbasis']>0: # shouldn't happen but avoid divbyzero
                atperc=avgval/df.iloc[i]['EDXbasis']
                colname='%'+key # key should be element name
                df=df.set_value(i, colname, atperc)
                # now calculated absolute error from error percentage 
                errcolname='err%'+key
                print('error column name is ', errcolname)
                thiserr=combinederr*atperc # at% error is at.% times com
                print('Abs error, combined % error, and at % are ', thiserr, combinederr,atperc)
                df=df.set_value(i, errcolname, thiserr) 
    # end of EDXbasis and at % calculation loop
       
    # organize data based on mycols template
    df=organizecolumns(df,mycols)
    df=df[df['SEMbasis']>0] # drop if basis is zero (often no entries in integlog)
    return df

def describecomps(compdf, Elements):
    ''' Quick output of subset of at.% data '''
    mycols=[]
    pd.set_option('display.float_format', lambda x: '%.1f' % x)
    for i, elem in enumerate(Elements):
        mycols.append('%'+elem)
    # only keep element subset in actual df
    mycols=[col for col in mycols if col in compdf]
    print('\n')
    print(compdf[mycols].describe())
    return

def printcomps(compdf, Elements, **kwargs):
    ''' Quick output of subset of at.% data '''
    mycols=['Filename']
    pd.set_option('display.float_format', lambda x: '%.1f' % x)
    for i, elem in enumerate(Elements):
        mycols.append('%'+elem)
    # only keep element subset in actual df
    mycols=[col for col in mycols if col in compdf]
    compdf['Filename']=compdf['Filename'].str.replace('.emsa','')
    if 'string' in kwargs:
        mystr=kwargs.get('string','')
        compdf=compdf[compdf['Filename'].str.contains(mystr)]
    compdf=compdf[mycols]
    # Renormalize to shown elements       
    print('\n')
    print(compdf[mycols])
    return compdf

def dropelement(EDXcomp, elem):
    ''' Remove element basis, comp and errors if never present (3 columns removed) '''
    dropcols=[elem, '%'+elem, 'err%'+elem]
    mycols=EDXcomp.columns.tolist()
    mycols=[col for col in mycols if col not in dropcols]
    EDXcomp=EDXcomp[mycols]
    return EDXcomp

def getelements(EDXcomp):
    ''' Retrieve actual element list in EDXcomp file '''
    mycols=EDXcomp.columns.tolist()
    elemlist=[col[1:] for col in mycols if col.startswith('%')]
    return elemlist

def changeelemsGUI(EDXquantparams, Elements):
    ''' Quick method of interactively selecting elements for plotting 
    has some hard-coded presets that can be changed using preset dictionaries below
    only elements with info in quant params csv files are selectable
    Note.. only tkinter variables exist after root.destroy
    '''
    # Enter default dictionaries for preset buttons (two currently available)
    preset1={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1,'FeL':1}
    preset2={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1,'O':1}
    # All available elements are those with entries in edxquantparams.csv
    elems=np.ndarray.tolist(EDXquantparams.element.unique()) 
    # Populate on/off with current element list
    elemdict={}
    for i, elem in enumerate(Elements):
        elemdict.update({elem:1})
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

def crosschecknames(TEMimagelog,EDXcomp, difflog):
    ''' Output of inconsistent names between EDX, TEM images and diff data sheets (diff is optional)
    typically for manual correction'''
    TEMnames=np.ndarray.tolist(TEMimagelog.Sample.unique())
    TEMnames=[col for col in TEMnames if str(col)!='nan']
    TEMnames=[col.lower() for col in TEMnames]
    EDXnames=np.ndarray.tolist(EDXcomp.Sample.unique())
    EDXnames=[col for col in EDXnames if str(col)!='nan']
    EDXnames=[col.lower() for col in EDXnames]
    # diffraction indexing is optional (can pass empty frame)
    if 'Sample' not in difflog:
        difflog['Sample']=''
    diffnames=np.ndarray.tolist(difflog.Sample.unique())
    diffnames=[col for col in diffnames if str(col)!='nan']
    diffnames=[col.lower() for col in diffnames]
    # don't want any mismatched diff names as identified phase info could be missed
    diffonly=[col for col in diffnames if col not in TEMnames+EDXnames]
    EDXonly=[col for col in EDXnames if col not in TEMnames+diffnames]
    if len(diffonly)>0:
        print('unique diff data samples:',', '.join(diffonly))
    if len(EDXonly)>0:
        print('unique EDX data samples:',', '.join(EDXonly))
    # output csv of filenames 
    return

def findphase(EDXcomp, difflog):
    ''' Find samples with known phase from diffraction and copy to EDX comp phase col'''
    temp=difflog.dropna(subset=['Phase'])
    mycols=EDXcomp.columns.tolist()
    EDXcomp=pd.merge(EDXcomp, temp, on=['Sample'], how='left', suffixes=('','_2'))
    EDXcomp['Phase']=EDXcomp['Phase'].replace('', np.nan)
    mask=EDXcomp[(pd.isnull(EDXcomp['Phase'])) & (pd.notnull(EDXcomp['Phase_2']))]
    EDXcomp.loc[mask.index, mask.columns]=mask
    EDXcomp['Phase']=EDXcomp['Phase_2']
    EDXcomp=EDXcomp[mycols]
    return EDXcomp

