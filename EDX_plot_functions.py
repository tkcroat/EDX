"""
Created on Aug 24th 2016
Plotting functions for SEM-EDX (or TEM-EDX ) spectra

@author: tkc
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re, math, sys
from matplotlib.backends.backend_pdf import PdfPages
import scipy
import scipy.stats # load in this sequence to get linregress working
from statsmodels.formula.api import ols # ordinary least squares
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX')
from EDX_import_functions import getelemmaps # called by reportmaps
import ternary
import tkinter as tk
# from sympy import Interval, Union # for overlapping plot range removal
#%%

font = {'size'   : 22}
plt.rc('font', **font)

MPL_STYLE = {
    "text.color":"k",
    "axes.labelcolor":"black",
    "axes.edgecolor":"0.4",
    "axes.facecolor":"white",   
    "xtick.color": "k",
    "ytick.color": "k",
    "figure.facecolor":"white",
    "figure.edgecolor":"white",
    "text.usetex":False,
    "axes.labelsize":"large"
}
plt.rcParams.update(MPL_STYLE)

def elemcompareplot(Integlog, elem1, elem2, thresh=0.1, errbars='xy'):
    ''' Pass two elements and make scatter plot of corrected counts
    useful for setting relative k-factors
    
    uses inner merge to select only subset with values from each df
	use either sample or filenumber'''
    #TODO have calccomp copy error in basis to allow errbars if basis=True
    el1=Integlog[Integlog['Element']==elem1]
    el2=Integlog[Integlog['Element']==elem2]

    fig, axes = plt.subplots(nrows=1, ncols=1) # axes is array
    # Merge dfs with comp1 and comp2 using inner join    
    compdata=pd.merge(el1, el2, how='inner', on=['Basename','Filenumber','Point'
            ,'Filename','Filepath','Sample','Comments'], suffixes=('','b'))
    
    compdata.plot.scatter(x='Correctedcounts', y='Correctedcountsb', ax=axes) # single plot axes has no [#,#]
    
    # linear regression: fitting, plot and add labels
    xdata=compdata['Correctedcounts'].as_matrix() # this data column as np array
    ydata=compdata['Correctedcountsb'].as_matrix()
    # 
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xdata, ydata) # imported from scipy.stats
    # set x range for linear plot 
    text1=str(round(slope,3))+' *x +' + str(round(intercept,3))
    text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,4))
    xmax=max(xdata)
    x=np.linspace(0,xmax,100) # setting range for 
    axes.text(0.025,0.9, text1, fontsize=12, transform=axes.transAxes)
    axes.text(0.025,0.8, text2, fontsize=12, transform=axes.transAxes)
    plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
    # Now test, plot and return outliers
    theseoutliers=returnoutliers(xdata.tolist(), ydata.tolist()) # index # 
    compdata= pd.concat([compdata, theseoutliers], axis=1, join='inner') # same length so just join by index
    # Add residual and pval to each compositional comparison line
    theseoutliers=compdata[compdata['Pval']<thresh] # now filter by threshold for outliers (with all cols)
    if not theseoutliers.empty:
        # error bars from errcorrcnts 
        if errbars=='xy':
            theseoutliers.plot.scatter(x='Correctedcounts', y='Correctedcountsb', 
                xerr='Errcorrcnts', yerr='Errcorrcntsb', ax=axes, color='r')
        elif errbars=='x': # plottable x error column exists
            theseoutliers.plot.scatter(x='Correctedcounts', y='Correctedcountsb', 
                    xerr='Errcorrcnts', ax=axes, color='r')
        elif errbars=='y': # plottable y error column exists
            theseoutliers.plot.scatter(x='Correctedcounts', y='Correctedcountsb', 
                    yerr='Errcorrcntsb', ax=axes, color='r')
        else: # no plottable errors for outliers
            theseoutliers.plot.scatter(x='Correctedcounts', y='Correctedcountsb', 
                    ax=axes, color='r')
    return compdata, theseoutliers

def reportmaps(NSScsvparams, Elements, PDFreport='NSSmaps_report.pdf'):
    '''Batch plot/PDF report of all selected elements for all NSS extracted csv image from NSS param log '''
    plt.ioff() # turn off interactive mode (for plot to pdf)
    with PdfPages(PDFreport) as pdf:
        for index, row in NSScsvparams.iterrows():
            thisrow=NSScsvparams.loc[[index]]
            elementmaps=getelemmaps(thisrow, Elements)
            fig=plotmaps(elementmaps, thisrow, savename='', imode=False) # no separate save and return figure
            pdf.savefig(fig)
    plt.ion() # turn interactive plotting back on
    return
        
def plotmaps(elementmaps, thismaprow, savename='', imode=True):
    ''' For plot arbitrary number of element maps (passed as list of numpy arrays) into single figure; designed for x-ray image maps 
    extracted from spectral images; name passed in thismap
    imode to false if called by PDF report loop'''
    if imode:
        plt.ion() # ensure interactive mode is on (on for single plots, off for PDF reports)
    # determine which of the elementmaps will be plotted 
    
    nummaps=len(elementmaps) # all passed element maps plotted (selection/filtering occurs with getelemmaps)
    # Determine shape of figure
    if nummaps<=3:
        numrows=1
    else:
        numrows=2
    numcols=math.ceil(nummaps/numrows)
    
    # fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False)
    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, squeeze=False)
    plottitle=thismaprow.iloc[0]['Basename']+': '+thismaprow.iloc[0]['SIname']+'.si'
    fig.suptitle(plottitle)
    
    for i, [elem, thismap] in enumerate(elementmaps): # thismap is element string followed by 512 x 512 array
        thisrow=i//numcols
        thiscol=i%numcols
        axindex=thisrow, thiscol # tuple to index axes 
        axes[axindex].set_aspect('equal')
        
        if elem=='Grey':
            axes[axindex].imshow(thismap, cmap='gray')  # make SE image grayscale
            axes[axindex].set_title('SE')
        else:
            axes[axindex].imshow(thismap, cmap='hot') # plots element map to correct subplot
            axes[axindex].set_title(elem)
    fig.tight_layout()
    # Hide empty subplots
    for i in range(0,numrows*numcols):
        if i>len(elementmaps)-1:
            thisrow=i//numcols
            thiscol=i%numcols
            axindex=thisrow, thiscol # tuple to index axes 
            axes[axindex].set_visible(False)            
    if savename!='':
        fig.savefig(savename) # optional saving of figure
    if imode==False:
        return fig
    else:
        return # just plot for interactive mode
    
def organizecomp(df):
    '''Get rid of common duplicated columns  '''
    removelist=['Projectb','Filenameb','FilePathb','Sampleb','Commentsb']
    singleelemlist=['Ca','Mg','Si','S']
    for i, val in enumerate(removelist):
        if val in df:
            df=df.drop(val,axis=1)
    for i, val in enumerate(singleelemlist): # don't drop basis as these will differ if comparing smdif and integ compositions
        if val+'amplb' in df:
            df=df.drop(val+'amplb',axis=1)
    return df
            
def returnoutliers(xdata, ydata):
    '''pass xcol and ycol as lists, makes plot and return outliers with pvals below specified threshold'''    
    # convert pandas series to lists
    regression= ols("data ~ x", data=dict(data=ydata, x=xdata)).fit()
    outliers=regression.outlier_test()    
    # df with cols as student_resid, unadj_p and bonf (bonferroni)
    colnames=['Resid','Pval','Bonf']
    outliers.columns=colnames # rename columns
    return outliers
    
def scattercompplot(comp1, comp2, elemlist, joinlist=['Filename'], thresh=0.1, basis=False, errbars='xy'):
    '''Pass two versions of composition calculation (using different lines or whatever) and compare 
    major elements using scatter graphs .. single point for each sample
    
    uses inner merge to select only subset with values from each df
	use either sample or filenumber'''
    #TODO have calccomp copy error in basis to allow errbars if basis=True
    elemlist=[re.match('\D+',i).group(0) for i in elemlist] 
    # strip number from peaks like Fe2 if present; columns will be element names (Fe) not peak names (Fe2)
    if basis==False: # use atomic % (which is the default), not basis for each element
        elemlist=['%'+s for s in elemlist]
    numregions=len(elemlist)
    
    # set nrows and ncols for figure of proper size
    cols=divmod(numregions,2)[0]+ divmod(numregions,2)[1]
    if numregions>1:
        rows=2
    else:
        rows=1        
    fig, axes = plt.subplots(nrows=rows, ncols=cols) # axes is array
    # merge dfs with comp1 and comp2 using inner join    
    compdata=pd.merge(comp1, comp2, how='inner', on=joinlist, suffixes=('','b'))
    mycols=compdata.dtypes.index # list of same columns
    mycols=mycols.tolist()
    #mycols=mycols.append('Element')
    outliers=pd.DataFrame(columns=mycols) # empty dataframe for outlying points
    fulldata=pd.DataFrame(columns=mycols)
    newcols=['Resid','Pval','Bonf', 'Element'] # single column for residuals but separate value needed for each row per element
    mycols.extend(newcols)
    for i, cname in enumerate(newcols):
        outliers[cname]=''
        fulldata[cname]=''
    for i, elem in enumerate(elemlist):
        # new version of base compositional data for each loop (otherwise reindexing problems)
        compdata=pd.merge(comp1, comp2, how='inner', on=joinlist, suffixes=('','b')) 
         # determine which subplot to use
        if (i+1)%2==1:
            rownum=0
        else:
            rownum=1
        colnum=int((i+1)/2.1)       
        xcol=elem
        ycol=elem+'b' # same element from second dataset
        if numregions==1: # deal with single subplot separately
            compdata.plot.scatter(x=xcol, y=ycol, ax=axes) # single plot axes has no [#,#]
        else:
            compdata.plot.scatter(x=xcol, y=ycol, ax=axes[rownum,colnum])
            # linear regression: fitting, plot and add labels
        xdata=compdata[elem].as_matrix() # this data column as np array
        colname=elem+'b'
        ydata=compdata[colname].as_matrix()
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xdata, ydata) # imported from scipy.stats
        # set x range for linear plot 
        text1=str(round(slope,3))+' *x +' + str(round(intercept,3))
        text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,4))
        if basis==False: # compositional plot so generally 0 to 1 for x and y axes
            xmax=max(max(xdata),max(ydata))*1.1 # set to slightly larger than max of dataset
        else: # if plotting elemental basis set to appropriate xmax (and y range will follow)
            xmax=max(xdata)
        x=np.linspace(0,xmax,100) # setting range for 
        if numregions==1: # deal with single subplot separately
            axes.text(0.025,0.9, text1, fontsize=12, transform=axes.transAxes)
            axes.text(0.025,0.8, text2, fontsize=12, transform=axes.transAxes)
            plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        else: # typical multiarea plot
            axes[rownum,colnum].text(0.025,0.9, text1, fontsize=12, transform=axes[rownum,colnum].transAxes)
            axes[rownum,colnum].text(0.025,0.8, text2, fontsize=12, transform=axes[rownum,colnum].transAxes)
            plt.axes(axes[rownum,colnum]) # set correct axes as active
            plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        # Now test, plot and return outliers
        theseoutliers=returnoutliers(xdata.tolist(), ydata.tolist()) # index # 
        # Add residual and pval to each compositional comparison line
        compdata= pd.concat([compdata, theseoutliers], axis=1, join='inner') # same length so just join by index
        compdata['Element']=elem # set after concat operation
        theseoutliers=compdata[compdata['Pval']<thresh] # now filter by threshold for outliers (with all cols)
        if not theseoutliers.empty and numregions==1:
            if errbars=='xy' and 'err'+elem in theseoutliers and 'err'+elem+'b' in theseoutliers: # plottable x and y error column
                theseoutliers.plot.scatter(x=xcol, y=ycol, xerr='err'+elem, yerr='err'+elem+'b', ax=axes, color='r')
            elif errbars=='x' and 'err'+elem in theseoutliers: # plottable x error column exists
                theseoutliers.plot.scatter(x=xcol, y=ycol, xerr='err'+elem, ax=axes, color='r')
            elif errbars=='y' and 'err'+elem in theseoutliers: # plottable y error column exists
                theseoutliers.plot.scatter(x=xcol, y=ycol, yerr='err'+elem, ax=axes, color='r')
            else: # no plottable errors for outliers
                theseoutliers.plot.scatter(x=xcol, y=ycol, ax=axes, color='r')
        if not theseoutliers.empty and numregions>1:
            if errbars=='xy' and 'err'+elem in theseoutliers and 'err'+elem+'b' in theseoutliers: # plottable x and y error column
                theseoutliers.plot.scatter(x=xcol, y=ycol, xerr='err'+elem, yerr='err'+elem+'b', ax=axes[rownum,colnum], color='r')
            elif errbars=='x' and 'err'+elem in theseoutliers: # plottable x error column exists
                theseoutliers.plot.scatter(x=xcol, y=ycol, xerr='err'+elem, ax=axes[rownum,colnum], color='r')
            elif errbars=='y' and 'err'+elem in theseoutliers: # plottable y error column exists
                theseoutliers.plot.scatter(x=xcol, y=ycol, yerr='err'+elem, ax=axes[rownum,colnum], color='r')
            else: # no plottable errors for outliers
                theseoutliers.plot.scatter(x=xcol, y=ycol, ax=axes[rownum,colnum], color='r')
        outliers=outliers.append(theseoutliers) # outliers from all elements
        fulldata=fulldata.append(compdata)
        # possibly could use ignore index but probably either is fine
    fulldata=fulldata[mycols] # put back in original column order
    outliers=outliers[mycols]
    fulldata=organizecomp(fulldata) # duplicated elements for same filenum have element-specific fit results
    outliers=organizecomp(outliers)
    return fulldata, outliers
    
def getplotboundaries(EDXfile, plotelems, EDXquantparams, colname='Energy'):
    ''' Gets typical boundary of plots of given line from EDXquantparams from plotelems but remove duplicates
    or remove if no data in range
    do not return range for duplicates (ie. Fe2 is in Fe plot range, Ca in C plot range so return 0s 
    defaults to checking for vals in energy column but can check peaks or other cols'''
    plotranges=[] # returns list of length 2 lists for valid elements
    for i, elem in enumerate(plotelems):
        ''' can use this to combine nearby elements
        if elem=='Ca' and 'C' in plotelems: # skip Ca if C is present
            continue
        '''
        thiselemdata=EDXquantparams[(EDXquantparams['element']==elem)]
        if len(thiselemdata)==1:
            thisrange=thiselemdata.iloc[0]['plotrange']
            try:
                match= re.finditer(r'\d+.\d+', thisrange)
                if match:
                    thisrange=[m.group(0) for m in match] # parse range into lower upper
                    # energy vals in keV for SEM-EDX so ensure boundaries are floats
                    thisrange=[float(i) for i in thisrange] 
                    EDXslice=EDXfile[(EDXfile[colname]>thisrange[0]) & (EDXfile[colname]<thisrange[1])]
                    if not EDXslice.empty:
                        plotranges.append(thisrange)
            except:
                pass # combining multiple lines in single window isn't an error
                # print('Duplicate plot range for ', elem)
        else:
            print ('Problem finding plot range for element ', elem)
    return plotranges

def getelemenergy(plotelems, plotrange, SEMquantparams):
    ''' Pass plotted data range and list of desired elements for plotting, return paired list of elemental and ideal energies 
    only returns element if in range of the given plot'''
    # elems=[str(s) for s in plotelems.split(' ')]
    elemlines=[] # energies list
    if '-' in plotrange: # determine range for plot if passed as hyphenated string
        plotrange=(float(plotrange.split('-')[0]),float(plotrange.split('-')[1])) # convert to 2 element tuple
    SEMquantparams=SEMquantparams[SEMquantparams['element'].isin(plotelems)] # select rows in element list
    SEMquantparams=SEMquantparams[(SEMquantparams['energy']>plotrange[0]) & (SEMquantparams['energy']<plotrange[1])]
    for index,row in SEMquantparams.iterrows() :
        elemlines.append([SEMquantparams.loc[index]['energy'],SEMquantparams.loc[index]['element']])
    return elemlines # list with [energy, element name] pairs for plotting

def setplotrange(plotrange, SEMfile):
    ''' Set range of plot based on element, numerical range or default to max range of spectrum  
    values are in eV (and often approx same as index for standard SEM-EDX '''
    # TODO ensure index # is approx same as energy in eV 
    plotranges=[]
    if '-' in plotrange: # determine range for plot
        plotranges.append([float(plotrange.split('-')[0]),float(plotrange.split('-')[1])])
    elif plotrange=='C':
        plotranges.append([2.36,3.16])
    elif plotrange=='Ca':
        plotranges.append([2.36,3.36]) 
    elif plotrange=='O':
        plotranges.append([4.70,5.40])
    elif plotrange=='Fe':
        plotranges.append([5.60,7.47])
    elif plotrange=='Mg':
        plotranges.append([1.145,1.225])
    elif plotrange=='Al':
        plotranges.append([1.350,1.430])
    elif plotrange=='Si':
        plotranges.append([1.570,1.650])
    else: # defaults to full data range (units)
        lower=SEMfile.Energy.min()        
        upper=SEMfile.Energy.max()
        plotranges.append([lower, upper])
    return plotranges

def getbackpoints(Backfitparamslog, SEMfileName):
    '''Get full list of index #s used for background fits  from backfitlog '''
    fullbackptslist=[]
    Backfit=Backfitparamslog[Backfitparamslog['Filename']==SEMfileName]
    for index, row in Backfit.iterrows():
        tempstr=Backfit.loc[index]['Backfitpts'] # range stored as string during fitting
        if type(tempstr)==list:
            fullbackptslist.extend(tempstr) 
        elif type(tempstr)==str:
            tempstr=tempstr.split('[')[1] # remove brackets from list
            tempstr=tempstr.split(']')[0]
            thesepts=[int(s) for s in tempstr.split(',')] # convert string to list of break index values
            fullbackptslist.extend(thesepts)        
        else:
            print('Problem opening background fit points for ', SEMfileName)
        # convert each string to list (one for each background region )
    return fullbackptslist
    
def reportcounts(EDXfiles, EDXquantparams, **kwargs):
    ''' Plot of list of files in EDXfiles, elements for plotting in plotelems, info on background fitted regions from backfitdf (optional) 
    background fits themselves stored with SEM files
    optional pass of backfitlog (w/ points defining region boundary for background fitting useful for troubleshooting fits)
    4/3/17 generalized version
    kwargs - 
        PDFname -- optional change of save name
        plotbackfits - bool ; add backfit col or not
        plotbackpts - optional scatter plot of background pts used in baseline fit 
             pass backfitdf dataframe
        plotelems - separate subplots around each chosen line energy
        labelelems - labeling of line energies in normal range plot
        xrange - optional range ... if not passed, plotted as separate elemental ranges
        yrange
    '''
    plt.ioff() # turn off interactive plot mode
    PDFname=kwargs.get('PDFname','countsback_elemreport.pdf') # use default or pass different name 
    with PdfPages(PDFname) as pdf:
        for index,row in EDXfiles.iterrows(): # iterrows avoids reindexing problems
            EDXfileName=EDXfiles.loc[index]['Filename']
            # find and open the EDX file
            csvname=EDXfileName.split('.')[0]+'.csv' # make sure it's csv not psmsa
            try: # works if in current path
                EDXfile=pd.read_csv(csvname) # reads entire spectra into df (all areas)
            except: # use full path
                path=EDXfiles.loc[index]['FilePath']
                fullname=path+"\\"+csvname
                EDXfile=pd.read_csv(fullname)
            if 'xrange' in kwargs:
                xrange=kwargs.get('xrange','0.250-10')
                plotranges=setplotrange(xrange, EDXfile) # list w/ single list of lower, upper
                plotelems=kwargs.get('plotelems',[]) # plotelems still needed for labeling
            elif 'plotelems' in kwargs:                
                plotelems=kwargs.get('plotelems',[]) # include default elements?
                plotranges=getplotboundaries(EDXfile, plotelems, EDXquantparams, colname='Energy') # returns plot ranges for all regions with data from plotelems
            else: # just use default energy range and no labelled elements
                plotranges=setplotrange('0.25-10', EDXfile)
            # optional plot of boundaries of backfit range from backfitlog 
            if 'plotbackpts' in kwargs:
                backfitdf=kwargs.get('plotbackpts','')
                thisfilebackpts=backfitdf[backfitdf['Filename']==EDXfileName]
                plotbackpts=True # don't add backpts without backfit col
            else:
                plotbackpts=False
            plotbackfits=kwargs.get('plotbackfits', True) # also plot subset of points used for backfit
            if len(plotranges)>0: # Skips if no data in selected element ranges (shouldn't happen for SEM-EDX)
                if plotbackpts==True: # this gets all the lower1, lower2, upper1, upper2 index point boundaries
                    # Just get all index #s from backfitpts for all regions
                    indexptslist=[]
                    for ind, row in thisfilebackpts.iterrows():
                        liststring=thisfilebackpts.loc[ind]['Backfitpts'] # can be string (if close/opened) or list 
                        # remove brackets and convert string to list of integer index #s
                        if isinstance(liststring,str):
                            liststring=liststring.split('[')[1] 
                            liststring=liststring.split(']')[0]
                            ptslist=[int(s) for s in liststring.split(',')]
                        else: # better be a list
                            ptslist=liststring
                        indexptslist.extend(ptslist)
                        # convert comma-separated string to actual 
                    indexptslist.sort()    
                # Determine # rows and columns from len(plotranges)
                numrows=min(len(plotranges),2) # 1 or 2 columns
                numcols=math.ceil(len(plotranges)/2)
                # Single plot for each SEM-EDX spectrum
                try:
                    if len(plotranges)<7:                    
                        fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
                    else:
                        fig, axes = plt.subplots(nrows=numcols, ncols=numrows, figsize=(16,9), squeeze=False) # switch to 2 row style for >7 subplots
                    # make plot title (include spec# and pt # if non-zero)
                    mytitle=EDXfiles.loc[index]['Basename']
                    if EDXfiles.loc[index]['Filenumber']!=1:
                        mytitle+='sp'+ str(EDXfiles.loc[index]['Filenumber'])
                    if EDXfiles.loc[index]['Point']!=1:
                        mytitle+='pt'+ str(EDXfiles.loc[index]['Point'])
                    plt.suptitle(mytitle)
                    
                    # now loop over the elemental plot ranges
                    for j, bounds in enumerate(plotranges):
                        [lower, upper]=bounds
                        if len(plotranges)<7:
                            thisrow=j%numrows
                            thiscol=j//numrows
                        else:
                            thiscol=j%numrows
                            thisrow=j//numrows
                        EDXslice=EDXfile[(EDXfile['Energy']>=lower) & (EDXfile['Energy']<=upper)] # already known that this isn't empty
                        EDXslice.plot(x='Energy', y='Counts', ax=axes[thisrow,thiscol]) # plot counts
                        if plotbackpts==True:
                            # Now add scatter plot points at fit region boundaries
                            backpts=EDXslice[EDXslice.index.isin(indexptslist)] # gets background fitted pts but only from this data slice
                            if not backpts.empty: # show fitted pts from counts
                                backpts.plot.scatter(x='Energy', y='Counts', ax=axes[thisrow,thiscol])
                        if plotbackfits==True:
                            EDXslice.plot(x='Energy', y='Backfit', ax=axes[thisrow,thiscol]) 
                        # Section for labeling plotelements (if range is passed, elems are still labeled)
                        if 'plotelems' in kwargs:
                            plotelems=kwargs.get('plotelems',[])
                            elemlines=getelemenergy(plotelems, bounds, EDXquantparams) # can pass plot range as lower,upper tuple
                            # list of tuples with energy,elemname
                            for k, elemtuple in enumerate(elemlines):
                                # elemtuple[0] is energy and [1] is element symbol
                                # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                try:
                                    axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                    yval=(EDXslice['Counts'].max()-EDXslice['Counts'].min())*0.9+EDXslice['Counts'].min() # setting y range
                                    axes[thisrow,thiscol].text(elemtuple[0],yval, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                                except:
                                    print('Problem labeling elements')
                    pdf.savefig(fig)
                    plt.close(fig) # closes recently displayed figure (since it's saved in PDF report)
                    print(EDXfileName,' plotted')
                except:
                    print('Problem plotting file ', EDXfileName)
            elif len(plotranges)==0:  # no data in selected plotelems range
                print('No data in plot range for all areas of ', EDXfileName)
    plt.ion() 
    return
        
def plotderivthresh(SEMfile, Backfitlog, Fitregionsdf, plotrange): # single file already selected by number
    '''Plot single SEM-EDX spectrum and its savgol derivative and background regions
    also available in report mode in reportcounts
    '''
    # TODO fix handing... plotrange is now list of lists (in keV)
    myplotrange=setplotrange(plotrange, SEMfile) # passes ev range or element or defaults to max range (returns tuple)
    
    fig, axes = plt.subplots(nrows=2, ncols=1) # axes is array
    
    # Slice to desired energy plotrange
    SEMslice=SEMfile[(SEMfile['Energy']>myplotrange[0]) & (SEMfile['Energy']<myplotrange[1])]

    Fitregions=unpackfitregs(Fitregionsdf) # turn df into list of lists
    
    # get max threshold for deriv knockout applicable for this plotted range
    threshold=findthreshold(Fitregions, plotrange)    
    # select only the subset in the background fitting region
    fullbackptslist=getbackpoints(Backfitlog, SEMfile) 
    Backfitpts=SEMslice.ix[[x for x in fullbackptslist]]
    Backfitpts=Backfitpts.dropna(subset=['Counts']) # above ix puts na in counts rather than dropping
    
    SEMslice.plot(x='Energy', y='Counts', ax=axes[0]) # single plot axes has no [#,#]
    SEMslice.plot(x='Energy', y='Backfit', color='r', ax=axes[0])
    Backfitpts.plot.scatter(x='Energy', y='Counts', color='b', ax=axes[0]) # single plot axes has no [#,#]
    SEMslice.plot(x='Energy', y='Savgol', ax=axes[1]) # single plot axes has no [#,#]
    axes[1].axhline(y=threshold, color='r') 
    axes[1].axhline(y=-threshold, color='r')
    axes[1].set_ylim([-threshold-10,threshold+10]) # setting y axes
    return

def reportsubdatamajor(paramlog, Integquantlog, PDFname='Subcounts_major_report'):
    ''' 2x3 SEM-EDX plot of subtracted data near O/FeL, Mg, Si, S, Ca and Fe
    pass pre-sliced params and integquantlog (with integration centers)  9/16/16
   '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            SEMfileName=paramlog.loc[index]['Filename']
            SEMfile=pd.read_csv(SEMfileName.replace('.psmsa','.csv')) # reads entire spectra into df (all areas)
            myplotrange=(SEMfile['Energy'].min(),SEMfile['Energy'].max()) # same range for all areas in spe            
            fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
            mytitle=str(SEMfileName)
            plt.suptitle(mytitle)
            # retrieve integration centers for all peaks from integquantlog
            match=Integquantlog[Integquantlog['Filename']==SEMfileName]
            # O/FeL regions
            if myplotrange[0] < 0.4 and myplotrange[1] > 0.9:
                SEMslice=SEMfile[(SEMfile['Energy']>0.4) & (SEMfile['Energy']<0.9)]
                if not SEMslice.empty:                    
                    SEMslice.plot(x='Energy', y='Subdata', ax=axes[0,0])
                ELmatch=match[match['Element']=='Fe2'] # find integ center for Fe2
                if len(ELmatch)==1:
                    center=ELmatch.iloc[0]['Energy']/100
                    axes[0,0].axvline(x=center-0.07, color='r') # assumes 15 channel integ window
                    axes[0,0].axvline(x=center+0.07, color='r')
                    axes[0,0].axvline(x=0.704, color='b') # ideal FeL
                    # Label with corrected counts and full 2 sigma error
                    if ELmatch.iloc[0]['Significance']>2:
                        corrcnts=ELmatch.iloc[0]['Correctedcounts']
                        errcorrcnts=ELmatch.iloc[0]['Errcorrcnts']
                        titlestring='FeL corrected counts ='+ '%.1f' %corrcnts +' +/- ' + '%.1f' %errcorrcnts
                        axes[0,0].set_title(titlestring, fontsize=10)
                    else:
                        titlestring='FeL region'
                        axes[0,0].set_title(titlestring, fontsize=10)
            # Mg region
            # TODO need scaling for Mg next to large Al peak
            if myplotrange[0] < 1.15 and  myplotrange[1] > 1.4:
                SEMslice=SEMfile[(SEMfile['Energy']>1.15) & (SEMfile['Energy']<1.4)]    
                if not SEMslice.empty:                    
                    SEMslice.plot(x='Energy', y='Subdata', ax=axes[1,0]) # Mg region
                ELmatch=match[match['Element']=='Mg']
                if len(ELmatch)==1:
                    center=ELmatch.iloc[0]['Energy']/100
                    axes[1,0].axvline(x=center-0.07, color='r') # assumes 15 channel integ window
                    axes[1,0].axvline(x=center+0.07, color='r')
                    axes[1,0].axvline(x=1.254, color='b') # ideal Mg
                    # Label with corrected counts and full 2 sigma error (if peak is significant)
                    if ELmatch.iloc[0]['Significance']>2:
                        corrcnts=ELmatch.iloc[0]['Correctedcounts']
                        errcorrcnts=ELmatch.iloc[0]['Errcorrcnts']
                        titlestring='Mg corrected counts ='+ '%.1f' %corrcnts +' +/- ' + '%.1f' %errcorrcnts
                        axes[1,0].set_title(titlestring, fontsize=10)
                    else:
                        titlestring='Mg region'
                        axes[1,0].set_title(titlestring, fontsize=10)
    			# Si region
            if myplotrange[0] < 1.64 and  myplotrange[1] > 1.94:   
                SEMslice=SEMfile[(SEMfile['Energy']>1.64) & (SEMfile['Energy']<1.94)]
                if not SEMslice.empty:                    
                    SEMslice.plot(x='Energy', y='Subdata', ax=axes[0,1]) # Si2 region                  
                ELmatch=match[match['Element']=='Si']
                if len(ELmatch)==1:
                    center=ELmatch.iloc[0]['Energy']/100
                    axes[0,1].axvline(x=center-0.07, color='r') # assumes 15 channel integ window
                    axes[0,1].axvline(x=center+0.07, color='r')  
                    axes[0,1].axvline(x=1.74, color='b') # ideal Si
                    # Label with corrected counts and full 2 sigma error (if peak is significant)
                    if ELmatch.iloc[0]['Significance']>2:
                        corrcnts=ELmatch.iloc[0]['Correctedcounts']
                        errcorrcnts=ELmatch.iloc[0]['Errcorrcnts']
                        titlestring=' Si corrected counts ='+ '%.1f' %corrcnts +' +/- ' + '%.1f' %errcorrcnts
                        axes[0,1].set_title(titlestring, fontsize=10)
                    else:
                        titlestring='Si region'
                        axes[0,1].set_title(titlestring, fontsize=10)
            # S region
            if myplotrange[0] < 2.2 and  myplotrange[1] > 2.5:        
                SEMslice=SEMfile[(SEMfile['Energy']>2.2) & (SEMfile['Energy']<2.5)]
                if not SEMslice.empty:                    
                    SEMslice.plot(x='Energy', y='Subdata', ax=axes[1,1]) # S region
                ELmatch=match[match['Element']=='S']
                if len(ELmatch)==1:
                    center=ELmatch.iloc[0]['Energy']/100
                    axes[1,1].axvline(x=center-0.07, color='r') # assumes 15 channel integ window
                    axes[1,1].axvline(x=center+0.07, color='r') 
                    axes[1,1].axvline(x=2.307, color='b') 
                    # Label with corrected counts and full 2 sigma error (if peak is significant)
                    if ELmatch.iloc[0]['Significance']>2:
                        corrcnts=ELmatch.iloc[0]['Correctedcounts']
                        errcorrcnts=ELmatch.iloc[0]['Errcorrcnts']
                        titlestring='S corrected counts ='+ '%.1f' %corrcnts +' +/- ' + '%.1f' %errcorrcnts
                        axes[1,1].set_title(titlestring, fontsize=10)
                    else:
                        titlestring='S region'
                        axes[1,1].set_title(titlestring, fontsize=10)
            # Ca region
            if myplotrange[0] < 3.6 and  myplotrange[1] > 3.9:         
                SEMslice=SEMfile[(SEMfile['Energy']>3.6) & (SEMfile['Energy']<3.9)]
                if not SEMslice.empty:                    
                    SEMslice.plot(x='Energy', y='Subdata', ax=axes[0,2]) # C/Ca region               
                ELmatch=match[match['Element']=='Ca']
                if len(ELmatch)==1:
                    center=ELmatch.iloc[0]['Energy']/100
                    axes[0,2].axvline(x=center-0.07, color='r') # assumes 15 channel integ window
                    axes[0,2].axvline(x=center+0.07, color='r') 
                    axes[0,2].axvline(x=3.69, color='b') 
                    # Label with corrected counts and full 2 sigma error (if peak is significant)
                    if ELmatch.iloc[0]['Significance']>2:
                        corrcnts=ELmatch.iloc[0]['Correctedcounts']
                        errcorrcnts=ELmatch.iloc[0]['Errcorrcnts']
                        titlestring='Ca corrected counts ='+ '%.1f' %corrcnts +' +/- ' + '%.1f' %errcorrcnts
                        axes[0,2].set_title(titlestring, fontsize=10)
                    else:
                        titlestring='Ca region'
                        axes[0,2].set_title(titlestring, fontsize=10)
            # Fe region
            if myplotrange[0] < 6.3 and  myplotrange[1] > 6.6:
                SEMslice=SEMfile[(SEMfile['Energy']>6.3) & (SEMfile['Energy']<6.6)]
                if not SEMslice.empty:                    
                    SEMslice.plot(x='Energy', y='Subdata', ax=axes[1,2]) # Fe region 
                ELmatch=match[match['Element']=='Fe']
                if len(ELmatch)==1:
                    center=ELmatch.iloc[0]['Energy']/100
                    axes[1,2].axvline(x=center-0.07, color='r') # assumes 15 channel integ window
                    axes[1,2].axvline(x=center+0.07, color='r')
                    axes[1,2].axvline(x=6.398, color='b')
                    # Label with corrected counts and full 2 sigma error (if peak is significant)
                    if ELmatch.iloc[0]['Significance']>2:
                        corrcnts=ELmatch.iloc[0]['Correctedcounts']
                        errcorrcnts=ELmatch.iloc[0]['Errcorrcnts']
                        titlestring='Fe corrected counts ='+ '%.1f' %corrcnts +' +/- ' + '%.1f' %errcorrcnts
                        axes[1,2].set_title(titlestring, fontsize=10)   
                    else:
                        titlestring='Fe region'
                        axes[1,2].set_title(titlestring, fontsize=10)                         
            pdf.savefig(fig)
            plt.close('all') # close all open figures
    plt.ion()
    return

def plotternary(df, ternelems, **kwargs):
    ''' Take compositional data, compute as 3-tuples and plot on ternary diagram 
    kwargs: symbolsize (default 40)
            colorgroup -- does groupby for multiples on this column - plotted in similar color
            title -- 
            symboltype
            plotunknown -- include unidentified phases in groupby ternary plot
            '''
    # Calculate ternary quantities for all in duplicates dataframe (allows each axis to be sum of elements)
    comp1=df.copy() # avoids problems with phase col alteration in source data
    hyphensplit = re.compile('(\+[a-zA-Z]+)').split
    ternlist=[part for img in ternelems for part in hyphensplit(img) if part]
    ternlist=[x.replace('+','') for x in ternlist]
    try:
        comp1['Tbasis']=comp1[ternlist].sum(axis=1) # sum a list of columns 
    except:
        print('Failed sum.. missing data for given element?')
    # calculate T1, T2, T3 values for ternary plot (auto normalized to unity)
    for i, val in enumerate(ternelems):
        num=str(i+1)
        if '+' not in val:
            comp1['T'+num]=comp1[val]/comp1['Tbasis']
        else:
            elems=[str(s) for s in val.split('+')]
            comp1['T'+num]=comp1[elems].sum(axis=1)/comp1['Tbasis']
    # create ternary plot
    figure,tax = ternary.figure(scale=1.0)
    fontsize=kwargs.get("fontsize",20)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=0.1, linewidth=0.5) # denser lines with higher multiple
    title=kwargs.get('title',"Ternary Composition Plot")
    tax.set_title(title, fontsize=fontsize)
    tax.left_axis_label(ternelems[2], fontsize=fontsize)
    tax.right_axis_label(ternelems[1], fontsize=fontsize)
    tax.bottom_axis_label(ternelems[0], fontsize=fontsize)
    tax.ticks(axis='lbr', linewidth=1, multiple=0.1) # set ticks
    # would be nice to be able to change size of tick labels
    # tax.ticks.label.set_fontsize(16) ... doesn't work
    tax.clear_matplotlib_ticks() # remove default Matplotlib axes (non-existent x and y axes)
    symbsize=kwargs.get('symbolsize',40)
    marktype=kwargs.get('symboltype','s') # symbol type
    if 'colorgroup' not in kwargs:
        plotpts=[]
        # Create list with 3 points as tuples (plus optional color)
        for index, row in comp1.iterrows():
            plotpts.append((comp1.loc[index]['T1'],comp1.loc[index]['T2'],comp1.loc[index]['T3']))
        tax.scatter(plotpts, marker='s', s=symbsize, color='b')  # s is point size
    # optional color groupby plot
    # http://stackoverflow.com/questions/26139423/plot-different-color-for-different-categorical-levels-using-matplotlib        
    else:
        # optional plotting of compositions w/o known phase from diffraction
        if kwargs.get('plotunknown',True):
            # replace blanks w/ unknown
            comp1['Phase']=comp1['Phase'].replace('','unknown')
            comp1['Phase']=comp1['Phase'].replace(np.nan,'unknown')
        groupcol=kwargs.get('colorgroup','')

        # pd groupby on passed column (done in uniform manner)
        groups=comp1.groupby(groupcol)
        # Define colorlist
        colorlist=['b','r','g','c','m','y','k', 'olive','pink','purple','plum']
        # Could also use 26 unique RGBA color tuple list but probably not needed
        colornum=0
        for key, group in groups:
            plotpts=[]
            for index, row in group.iterrows():
                plotpts.append((group.loc[index]['T1'],group.loc[index]['T2'],group.loc[index]['T3']))
            tax.scatter(plotpts, marker=marktype, s=symbsize, color=colorlist[colornum], label=key)  # s is point size
            colornum+=1
        tax.legend()
    return # return with ternary compositional data

def reportcountspeakfits(SEMfiles, Backfitparamslog, plotrange, plotelems, SEMquantparams):
    ''' Report of SEM counts + background and deriv of all files in passed log
    Fitregionsdf is loaded SEM background regions csv file
    plotrange as string in eV (aka index #); 
    Threshold sets deriv value above which points are eliminated from backfit (use same threshold as for prior background fits)
    
    '''            
    threshold=15   
    plt.ioff() # turn off interactive mode
    with PdfPages('countsback_report.pdf') as pdf:
        for index,row in SEMfiles.iterrows(): # iterrows avoids reindexing problems
            SEMfileName=SEMfiles.loc[index]['Filename']
            csvname=SEMfileName.split('.')[0]+'.csv' # make sure it's csv not psmsa
            try: # works if in current path
                SEMfile=pd.read_csv(csvname) # reads entire spectra into df (all areas)
            except: # use full path
                path=SEMfiles.loc[index]['FilePath']
                fullname=path+"\\"+csvname
                SEMfile=pd.read_csv(fullname)                 
            myplotrange=setplotrange(plotrange, SEMfile) # could be done globally but might hit data out of range problem
            fig, axes = plt.subplots(nrows=2, ncols=1) # axes is array
            mytitle=SEMfiles.loc[index]['Basename']+'sp#'+ str(SEMfiles.loc[index]['Filenumber'])
            plt.suptitle(mytitle)
            SEMslice=SEMfile[(SEMfile['Energy']>myplotrange[0]) & (SEMfile['Energy']<myplotrange[1])]
            
            # Get subset of points used for background fit for scatter plot
            fullbackptslist=getbackpoints(Backfitparamslog, SEMfileName)            
            Backfitpts=SEMslice.ix[[x for x in fullbackptslist]]
            Backfitpts=Backfitpts.dropna(subset=['Counts']) # above ix puts na in counts rather than dropping
    
            SEMslice.plot(x='Energy', y='Counts', ax=axes[0]) # single plot axes has no [#,#]
            SEMslice.plot(x='Energy', y='Backfit', color='r', ax=axes[0])
            
            # find elemental lines in this keV range
            elemlines=getelemenergy(plotelems, plotrange, SEMquantparams)
            for i, elemtuple in enumerate(elemlines):
                # elemtuple[0] is energy and [1] is element symbol
                axes[0].axvline(x=elemtuple[0], color='b') # O line
                axes[0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val
            Backfitpts.plot.scatter(x='Energy', y='Counts', color='b', ax=axes[0]) # single plot axes has no [#,#]
            SEMslice.plot(x='Energy', y='Savgol', ax=axes[1]) # single plot axes has no [#,#]
            axes[1].axhline(y=threshold, color='r') 
            axes[1].axhline(y=-threshold, color='r')
            axes[1].set_ylim([-threshold-10,threshold+10]) # setting y axes
            pdf.savefig(fig)
            plt.close('all') # close all open figures    
    plt.ion() # turn interactive back on
    return

def EDXplot_gui(EDXlog, Elements, Backfitlog, EDXquantparams):
    ''' tk interface for args/kwargs of single AES interactive plot
	filtering by filenumber, filename incorporated in tk interface 
    all args/dataframes must be passed through to plot functions 
    '''
    # first print out existing info in various lines
    root = tk.Tk()
    root.title('Interactive EDX plot interface')
    # Set up all the tk variables 
    filterstr=tk.StringVar() # comma separated or range of filenumbers for plot or string for sample name
    filenums=tk.StringVar()  # For file number filtering (single or multiple comma separated)
    filtercol=tk.StringVar()  # to which column is optional filter applied
    filtercol.set('Filename') # set default to filename (not sample)
    xrangestr=tk.StringVar()  # energy range in eV 
    xrangestr.set('0.3-10')

    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    backfitbool.set(1) # default true for EDX background fits 
    backptbool=tk.BooleanVar() 
    elemstr=tk.StringVar()
    mytext=', '.join(Elements) # elements for labelling 
    elemstr.set(mytext)
    legendsize=tk.IntVar()
    legendsize.set(8)
    plotelemsbool=tk.BooleanVar()  # optional labelling of elements
    plotelemsbool.set(True) # default to true
    choice=tk.StringVar()  # plot or abortw

    # Optional filtering of chosen spectra by filenumber or string
    tk.Label(root, text='Filename string filter').grid(row=0, column=0)
    tk.Label(root, text='Filenumber(s) filter').grid(row=1, column=0)
    tk.Checkbutton(root, variable=filtercol, text='Filter on sample column').grid(row=2, column=0)
    tk.Entry(root, textvariable=filterstr).grid(row=0, column=1)
    tk.Entry(root, textvariable=filenums).grid(row=1, column=1)
    
    backfit=tk.Checkbutton(root, variable=backfitbool, text='Plot background fits?')
    backfit.grid(row=0, column=2) # can't do immediate grid or nonetype is returned
    backpt=tk.Checkbutton(root, variable=backptbool, text='Plot background points?')
    backpt.grid(row=1, column=2)   
    tk.Checkbutton(root, variable=plotelemsbool, text='Label element peaks?').grid(row=2, column=2)
     
    rownum=3
    # Choose x plotting range
    tk.Label(root, text='Set x data range').grid(row=rownum, column=0)
    tk.Entry(root, textvariable=xrangestr).grid(row=rownum, column=1)
    rownum+=1

    tk.Label(root, text='Elements:').grid(row=rownum, column=0)
    elementry=tk.Entry(root, textvariable=elemstr)
    elementry.grid(row=rownum, column=1)

    rownum+=1
    tk.Label(root, text='Legend size').grid(row=rownum, column=0)
    tk.Entry(root, textvariable=legendsize).grid(row=rownum, column=1)
    rownum+=1
    # TODO option to reselect labeled elemental peaks 

    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    a=tk.Button(root, text='Plot')
    a.bind('<Button-1>', plot)
    a.grid(row=rownum, column=1)
    
    a=tk.Button(root, text='Abort')
    a.bind('<Button-1>', abort)
    a.grid(row=rownum, column=2)

    root.mainloop()
        
    mychoice=choice.get()
    
    if mychoice=='plot':
        kwargs={} # kwargs for EDXplot1 
        if filenums.get()!='': # filter by entered filenumber(s)
            filenums=[int(i) for i in filenums.get().split(',')]
            EDXlog=EDXlog[EDXlog['Filenumber'].isin(filenums)]
        if filterstr.get()!='':            
            if filtercol.get()=='Filename':
                EDXlog=EDXlog[EDXlog['Filename'].str.contains(filterstr.get())]
            elif filtercol.get()=='Sample':
                EDXlog=EDXlog[EDXlog['Sample'].str.contains(filterstr.get())]
        # Can handle multiple spectra in single plot
        if len(EDXlog)==0:
            print('No files remain after', filterstr.get(),' filter applied.')
            # return         
        # Handle plot x range choice (used for all but peaks report)
        # modified from Auger but probably overly complicated for EDX scan 
        plotrange=[]
        tempstr=xrangestr.get()
        if '-' in tempstr: # should be of form 0.3-10
            plotrange.extend([float(f) for f in tempstr.split('-')[0:2]]) # parse into strings if necessary (list must be passed)
        else: 
            print('Problem parsing ev range... using default')
            plotrange=[[0.3,10]]
         # Pass kwargs to single EDX plotting window
        if backfitbool.get():
            kwargs.update({'plotbackfit':True}) # plotting of background
        if backptbool.get():
            kwargs.update({'plotbackpts':True}) # plot pts used to fit background
            # these are retrieved from backfitlog
        if plotelemsbool.get(): # optional labeling of elements
            kwargs.update({'plotelems':Elements}) # pass elements list        
        kwargs.update({'legend':int(legendsize.get())}) # always pass legend size
        EDXplot1(EDXlog, plotrange, Backfitlog, EDXquantparams, **kwargs)
    return kwargs

def EDXplot1(EDXlog, plotrange, Backfitlog, EDXquantparams, **kwargs):
    ''' tk launched single plot window (modeled on AESplot1 and its interface)
    
    Plot a single frame over desired range with passed filenumbers/ filenames
    plotelems kwarg is list of elements who's ideal position will be shown 
    as blue vert line on plot
    filenumbers - Single # or list of numbers
    kwargs:
        plotelems - set of elements to label 
        plotbackfit -- optional background fitting (data direct from EDXfile) 
        plotbackpts -- optional plot of pts used for background fits
        legend - font size for added legend
        # TODO add gaussian fit option to Counts plots (using params from Integquantlog )
        '''
    fig, axes = plt.subplots(nrows=1, ncols=1, squeeze=False) # axes is array
    plt.ticklabel_format(axes='y', style='sci', scilimits=(-3,3))
    # Use color and linestyles to differentiate filenumber(s) and areanumber(s)
    mylegend=[]
    colorlist=['b','r','g','c','m','y','k', 'olive','pink','purple']
    '''
    TESTING    
    index=0  row=EDXlog.loc[index]
    '''
    filecount=0  
    for index, row in EDXlog.iterrows():
        thisfilenum=row.Filenumber
        try:
            filename=row.Filename.split('.')[0]+'.csv' # psmsa or emsa
            EDXfile=pd.read_csv(filename)
            if EDXfile.empty: # file not found
                continue
            # Also check in sub directory
            filecount+=1
            EDXslice=EDXfile[(EDXfile['Energy']>plotrange[0]) & (EDXfile['Energy']<plotrange[1])]

            # Default is color for areas and linestyle for 
            plkwargs={'color':colorlist[filecount%10]}

            mylegend.append(filename.split('.')[0])

            EDXslice.plot(x='Energy', y='Counts', ax=axes[0,0], **plkwargs)
            # Optional additional plotting of backgrounds (only for direct counts plotting)
            if 'plotbackfit' in kwargs:
                # plot background itself always in axes 0, 0
                EDXslice.plot(x='Energy', y='Backfit', ax=axes[0,0])
            if 'plotbackpts' in kwargs:
                # Scatter plot points over which background was fitted
                indexptslist=getbackfitpts(Backfitlog, row.Filename)
                backpts=EDXslice[EDXslice.index.isin(indexptslist)] # gets background fitted pts but only from this data slice
                if not backpts.empty: # show fitted pts from counts
                    backpts.plot.scatter(x='Energy', y='Counts', ax=axes[0,0])
        except:
            print('Error plotting file', str(thisfilenum))
    # Alter legend with filenumber+ areanumber
    if 'legend' in kwargs:
        fs=kwargs.get('legend',0)
        axes[0,0].legend(mylegend, loc='best', fontsize=fs)
        try:
            axes[1,0].legend(mylegend, loc='best', fontsize=fs)
        except:
            pass
    for axis in ['top', 'bottom','left','right']:
        axes[0,0].spines[axis].set_linewidth(2.5) # change axes thickness
        try: # sometimes plotting both counts and deriv
            axes[1,0].spines[axis].set_linewidth(2.5) 
        except:
            pass
    # Add elemental peaks as vert lines
    plotelems=kwargs.get('plotelems',[])
    if len(plotelems)>0:
        elemlines = findelemlines(plotelems, plotrange, EDXquantparams)
        ypos=axes[0,0].get_ylim()[1]*.95
        for i, [val,elem] in enumerate(elemlines):
            axes[0,0].axvline(x=val, color='b')
            xpos=val+0.05
            # Add element text label  
            axes[0,0].text(xpos,ypos, elem, fontsize=12, rotation=90)
    return
    

def findelemlines(plotelems, xrange, EDXquantparams):
    ''' Pass list of element peaks and plot range; return energies if in range
    ''' 
    # plotelems=plotelems.replace('Fe','Fe Fe1 Fe2') # choose all Fe lines
    # elems=[str(s) for s in plotelems.split(' ')]
    
    elemlines=[] # energies list    
    peaks=EDXquantparams[EDXquantparams['element'].isin(plotelems)] # select rows in element list
    peaks=peaks[(peaks['energy']>xrange[0]) & (peaks['energy']<xrange[1])]
    for index, row in peaks.iterrows():
        elemlines.append([row.energy, row.element])
    return elemlines

def getbackfitpts(Backfitlog, filename):
    ''' Find points over which background fit was performed (for scatter plots)
    called by various plot and reporting functions '''
    
    thisfilebackpts=Backfitlog[(Backfitlog['Filename']==filename)]
    if len(thisfilebackpts)==0:
        print('Backfitpts not found for ', filename)
        return []
    backptslist=[]
    for index, row in thisfilebackpts.iterrows():
        if isinstance(row.Backfitpts, str): # should be a string
            pts=row.Backfitpts.replace('[','').replace(']','')
            pts=pts.split(',')
            pts=[int(i) for i in pts]
            backptslist.extend(pts)
        elif isinstance(row.Backfitpts, list): # each is a list
            backptslist.extend(row.Backfitpts)
    return backptslist

def reportSEMpeaks(paramlog, plotelems, SEMquantparams, addgauss=True, PDFname='peak_report.pdf'):
    ''' Two col plot of low and high energy SEM-EDX from logfile
    plotelems is list of elements to label on the plot
    '''
    # TODO put line on center of integration
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows():
            SEMfileName=paramlog.iloc[index]['Filename']
            if '.psmsa' in SEMfileName: # always open csv not psmsa
                SEMfileName=SEMfileName.split('.')[0]+'.csv'
            SEMfile=pd.read_csv(SEMfileName) # reads entire spectra into df (all areas)
            myplotrange=(SEMfile['Energy'].min(),SEMfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.iloc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            basename=Params.Basename # retrieve filenumber
            ptnum=Params.Point # retrieve filenumber

            fig, axes = plt.subplots(nrows=2, ncols=1) # 2 by 3 axes array
            mytitle=str(basename)+ 'sp'+ str(filenumber)+'pt'+str(ptnum)
            plt.suptitle(mytitle)
            # low energy region
            if myplotrange[0] < 0.1 and  myplotrange[1] > 4.0:
                SEMslice=SEMfile[(SEMfile['Energy']>0.1) & (SEMfile['Energy']<4.0)]
                thisrange='0.1-4.0'
                # find elemental lines in this keV range
                elemlines=getelemenergy(plotelems, thisrange, SEMquantparams)
                for i, elemtuple in enumerate(elemlines):
                    # elemtuple[0] is energy and [1] is element symbol
                    axes[0].axvline(x=elemtuple[0], color='b') # O line
                    axes[0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val
                if not SEMslice.empty:                    
                    SEMslice.plot(x='Energy', y='Subdata', ax=axes[0]) # S region
                    if addgauss==True and 'Gauss' in SEMslice.dtypes.index: # ensure Gaussian fit col exists
                        SEMslice.plot(x='Energy', y='Gauss', ax=axes[0]) # S region
            # high energy region
            if myplotrange[0] < 4.0 and  myplotrange[1] > 8.0:         
                SEMslice=SEMfile[(SEMfile['Energy']>4.0) & (SEMfile['Energy']<8.0)]
                thisrange='4.0-8.0'
                # find elemental lines in this keV range
                elemlines=getelemenergy(plotelems, thisrange, SEMquantparams)
                for i, elemtuple in enumerate(elemlines):
                    # elemtuple[0] is energy and [1] is element symbol
                    axes[0].axvline(x=elemtuple[0], color='b') # O line
                    axes[0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val
                axes[1].axvline(x=6.4, color='b') # O line
                axes[1].text(6.4,-10, 'Fe',rotation=90)
                if not SEMslice.empty:                    
                    SEMslice.plot(x='Energy', y='Subdata', ax=axes[1]) # C/Ca region               
                    if addgauss==True and 'Gauss' in SEMslice.dtypes.index: 
                        SEMslice.plot(x='Energy', y='Gauss', ax=axes[1]) # C/Ca region                    
            pdf.savefig(fig)
            plt.close('all') # close all open figures
    return

# LEGACY FUNCTIONS (most likely)
def unpackfitregs(df, SEMfile):
    ''' df is fitregionsdf; Loaded data frame has ev range, list of background regions and fit type
    unpack from dataframe into list of each 
    1) fitrange is total ev range (i.e.0-100), 2 )fitpts are index #s (or energy in eV) of regirons commonly without peaks
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
            mandminlist=findminpoints(df.iloc[i]['Minpoint'], SEMfile)
        Fitregions.append([df.iloc[i]['Fitrange'],indexrange, mandminlist, df.iloc[i]['Fittype'], df.iloc[i]['Threshold']])
    return Fitregions
    
def findthreshold(Fitregions, plotrange):
    ''' Find best derivative threshold for this plot 
    Fitregions is unpacked list of lists for each commonly fitted region, plotrange is normally text str with numerical range in keV '''
    if '-' in plotrange: # determine range for plot converted to eV (from keV)
        try:
            myplotrange=(int(float(plotrange.split('-')[0])*100),int(float(plotrange.split('-')[1])*100))
            for i in range(0,len(Fitregions)): # proceeding from highest threshold (low eV) to lowest
                fitrange=Fitregions[i][0]
                upper=int(fitrange.split('-')[1])
                if upper>myplotrange[0]: #
                    threshold=int(Fitregions[i][3])
                    return threshold
        except:
            threshold=20 # just use max value
            return threshold
    else: # just set to 0-1000 which will get max value
        threshold=20
        return threshold
        
def findminpoints(rangestrings, SEMfile):
    '''Pass index energy range (or ranges) and find/return minimum(s) in SEM file within these ranges 
    rangestring is index # range as string '''
    stringlists=[str(s) for s in rangestrings.split(',')]
    addlist=[]
    for i, val in enumerate(stringlists):
        thisminrange=rangefromstring(val)
        dfslice=SEMfile[min(thisminrange):max(thisminrange)] # slice based on this index range
        thismin=dfslice['Counts'].idxmin()
        addlist.append(int(thismin))
    return addlist # list of index #s of minimums within the ranges specified 
      
def rangefromstring(x):
    ''' Pass plot range string and return  '''
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

''' Legacy and problem with sympy Union/Interval?
def getboundaries(EDXfile, plotrange, EDXquantparams):
    # Slightly modified version of getplotboundaries... same basic idea
    
    plotranges=[] # returns list of length 2 lists for valid elements
    # optional filtering colname (normally used for setting range for Peaks)

    for i, thisrange in enumerate(plotrange):
        # Already split as list with xmin, xmax 
        if isinstance(thisrange, list):
            plotranges.append([float(thisrange[0]),float(thisrange[1])])
        if '-' in thisrange: # direct ev range specified in length 1 list (not list of elements)
            plotranges.append([float(thisrange.split('-')[0]),float(thisrange.split('-')[1])])
        else: # if not hyphenated it's a named Auger peak 
            thiselemdata=EDXquantparams[EDXquantparams['element']==thisrange]
            if len(thiselemdata)==1:
                thisrange=thiselemdata.iloc[0]['plotrange']
                try:
                    match= re.finditer(r'\d+', thisrange)
                    if match:
                        thisrange=[m.group(0) for m in match] # Parse range into lower upper
                        thisrange=[int(i) for i in thisrange]
                        EDXslice=EDXfile[(EDXfile['Energy']>thisrange[0]) & (EDXfile['Energy']<thisrange[1])]
                        # now drop na on passed column
                        if not EDXslice.empty:
                            plotranges.append(thisrange)
                except:
                    print ('Problem finding plot range for element ', thisrange)
                    pass # combining multiple lines in single window isn't an error
            # knock out duplicated plot ranges (i.e. just plot Fe and Fe2 together)
    # now knock out any duplicate or overlapping plot ranges, such as C and Ca, Fe and Fe2
    # this is controllable via plotrange specified in AESquantparams
    if len(plotranges)>1:
        plotranges=combineranges(plotranges)
    return plotranges

def combineranges(mylist):
    # Remove exact duplicates and combine ranges if overlapping
    # uses sympy functions Interval (real intervals) and Union
    intervals=[Interval(begin,end) for [begin, end] in mylist]
    u=Union(*intervals)
    # Convert all the FiniteSets (distinct elements in unions) into nested list with min, max values for subplots
    for i in range(0, len(u.args)):
        newrange=[]
        for i in range(0, len(u.args)):
            pair=[]
            for i in iter(u.args[i].boundary):
                pair.append(int(i))
            newrange.append(pair)
    return newrange
'''