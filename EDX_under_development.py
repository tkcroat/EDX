# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 14:08:57 2016

@author: tkc
"""

# SEM under development
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os, re
from matplotlib.backends.backend_pdf import PdfPages
import scipy
import scipy.stats # load in this sequence to get linregress working
import datetime
from io import StringIO
from math import factorial # used by Savgol matrix
from scipy import optimize
from PIL import Image, ImageDraw, ImageFont
from decimal import Decimal 
import os
if 'C:\\Users\\tkc\\Documents\\Python_Scripts' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts')
import tkinter as tk
import datetime

EDXfile=pd.read_csv('C:\Temp\SiC\MD2d_11Jun10\MD2d18Jun10(1)_pt1.csv')


# Determine overlap between pathological peaks  (origin gaussian fits)
y0=1.16
xc1=9.45
w1=0.152
A1=62.7
xc2=9.257
w2=0.18
A2=51.8

gauss1=(A1/w1*np.sqrt(math.pi()/2)*exp
        
        
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def remove_freq(train):
    ''' Turn removed and added points into frequency histograms
    '''
    # make single 1D array with all vals concatenated
    vals=np.empty(0)
    # testing  row=train.iloc[0]
    for index, row in train.iterrows():
        newvals=row.Xrem.replace('[','').replace(']','')
        if newvals=='':
            continue
        new=newvals.split(',')
        new=[int(i) for i in new]
        vals=np.append(vals, new)
    remove=stats.itemfreq(vals)
    remove.sort(axis=1)
    return remove

def plotternary(df, ternelems, **kwargs):
    ''' Take compositional data, compute as 3-tuples and plot on ternary diagram 
    kwargs: symbolsize (default 40)
            colorgroup -- does groupby for multiples on this column - plotted in similar color
            symboltype -            
            title -- 
            
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
    fontsize=20
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=0.1, linewidth=0.5) # denser lines with higher multiple
    title=kwargs.get('title',"Ternary Composition Plot")
    tax.set_title(title, fontsize=fontsize)
    tax.left_axis_label(ternelems[2], fontsize=fontsize)
    tax.right_axis_label(ternelems[1], fontsize=fontsize)
    tax.bottom_axis_label(ternelems[0], fontsize=fontsize)
    tax.ticks(axis='lbr', linewidth=1, multiple=0.1) # set ticks
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
        if kwargs.get('plotundetermined',False):
            # replace blanks w/ unknown
            comp1['Phase']=comp1['Phase'].replace('','undetermined')
            comp1['Phase']=comp1['Phase'].replace(np.nan,'undetermined')
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



from sklearn.lda import LDA
from sklearn import decomposition
from sklearn import datasets
from sklearn.preprocessing import StandardScaler

fig, axes = plt.subplots(nrows=1, ncols=1, squeeze=False)

coords1=[] # Global coords list for matching correspondence points between im1 and im2
coords2=[]
cid = fig.canvas.mpl_connect('button_press_event', onclick)
cid2 = fig.canvas.mpl_connect('button_release_event', onrelease)

coords1=[(1.2, 1600)]
coords2=[(1.8, 2400)]
fig.canvas.mpl_disconnect(cid)
fig.canvas.mpl_disconnect(cid2)

len(coords1)


EDXcols=pdutils.pickcols_tk(EDXcomp)
EDXcols=pickcols_tk(EDXcomp)

iris = datasets.load_iris()
X=iris.data
Y=iris.target # categories as number (0- setosa 1- versicolor 2-virginica)

# Attempt PCA (unlabelled) and LDA (labelled) on same dataset
X=EDXcols
pca = decomposition.PCA(n_components=3)
pca.fit(X)
Y = pca.transform(X)

X_std = StandardScaler().fit_transform(X) # 
np.mean(X_std[:0])
np.mean(X_std[0])
Ser=X_std.loc[]
Ser=X_std.iloc[[3]]
       
def findbadpts(coords1, coords2, backpts):
    ''' Quick way of finding/returning troublesome points for backfitdf 
    backpts returned from single plot '''
    for i in range(0,len(coords1)):
        print(i)
    myranges=
    pts=pts[]
    test=backpts[backpts['Energy'] 

def findbadpts():
    ''' need an interactive mpl method of removing points from std backfit regions ... create custom 
    backfit version for later use   '''
    # TODO go through a stack of spectra and select/keep track of keV of bad points ... return and 
    # make custom backfit regions 
    # mpl_connect ... check match_images_workflow
    
    def onclick(event):
        ''' Grabs pixel coords from click on image... defind locally to work with global coords'''
        global ix, iy
        ix, iy = event.xdata, event.ydata
        print ('x = %d, y = %d'%(ix, iy))
        global coords1
        coords1.append((ix, iy))
        return coords1
    def onrelease(event):
        ''' Grabs pixel coords from click on image... defind locally to work with global coords'''
        global ix, iy
        ix, iy = event.xdata, event.ydata
        print ('x = %d, y = %d'%(ix, iy))
        global coords2
        coords2.append((ix, iy))
        return coords2

Ser=EDXfiles.iloc[1]
Ser.Filename in Backfitlog.Filename.unique()
plotcounts(Ser, EDXquantparams, **kwargs)

len(Backfitlog.Filename.unique())

backpt=plotcounts(Ser, EDXquantparams, **kwargs)

def plotcounts(Ser, EDXquantparams, **kwargs):
    ''' Single plot of chosen EDX file, elements for plotting in plotelems, info on background fitted regions from backfitdf (optional) 
    background fits themselves stored with SEM files
    optional pass of backfitlog (w/ points defining region boundary for background fitting useful for troubleshooting fits)
    kwargs - 
        backfitdf -  optional scatter plot of background pts used in baseline fit 
        plotelems - separate subplots around each chosen line energy
        labelelems - labeling of line energies in normal range plot
        xrange - optional range ... if not passed, plotted as separate elemental ranges
        yrange
    '''
    EDXfileName=Ser['Filename']
    # find and open the EDX file
    csvname=EDXfileName.split('.')[0]+'.csv' # make sure it's csv not psmsa
    try: # works if in current path
        EDXfile=pd.read_csv(csvname) # reads entire spectra into df (all areas)
    except: # use full path
        path=Ser['FilePath']
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
        plotelems=[]
    # optional plot of boundaries of backfit range from backfitlog 
    if 'backfitdf' in kwargs:
        backfitdf=kwargs.get('backfitdf','')
        thisfilebackpts=backfitdf[backfitdf['Filename']==EDXfileName]
        plotbackpts=True
    else:
        plotbackpts=False
        thisfilebackpts=pd.DataFrame()
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
            mytitle=Ser['Basename']
            if Ser['Filenumber']!=1:
                mytitle+='sp'+ str(Ser['Filenumber'])
            if Ser['Point']!=1:
                mytitle+='pt'+ str(Ser['Point'])
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
                    EDXslice.plot(x='Energy', y='Backfit', ax=axes[thisrow,thiscol]) 
                # Section for labeling plotelements (if range is passed, elems are still labeled)
                # 
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
            print(EDXfileName,' plotted')
        except:
            print('Problem plotting file ', EDXfileName)
    elif len(plotranges)==0:  # no data in selected plotelems range
        print('No data in plot range for all areas of ', EDXfileName)
    return thisfilebackpts
    
clf=LDA()
clf.fit(X, y)

# Copy sample names between EDXcomp and EDXlog (phase typically only in EDXcomp)

synclogs('EDXparamlog.csv', 'EDXcomp.csv', colnames)
synclogs('EDXcomp.csv', 'EDXcomp_metalsbasis.csv', colnames)
    
def printcomps(df, Elements, **kwargs):
    ''' Quick output of subset of at.% data 
    string search filters w/ wildcard and compositional filters '''
    mycols=[]
    compdf=df.copy()
    pd.set_option('display.float_format', lambda x: '%.1f' % x)
    for i, elem in enumerate(Elements):
        mycols.append('%'+elem)
    mycols=[col for col in mycols if col in compdf]
    compdf['Total']=0
    for i, col in enumerate(mycols):
        compdf['Total']+=compdf[col]
    # Only keep element subset in actual df
    compdf['Filename']=compdf['Filename'].str.replace('.emsa','')
    # choose subset that match substring
    if 'string' in kwargs:
        mystr=kwargs.get('string','')
        compdf=compdf[compdf['Filename'].str.contains(mystr)]
    # Renormalize to shown elements
    for i, col in enumerate(mycols):
        compdf[col]=compdf[col]*100/compdf['Total']
    mycols.append('Filename')
    # Apply compositional filters if present
    for i, col in enumerate(mycols):
        if col in kwargs: # > or < and 
            limstr=kwargs.get(col,'')
            if limstr[0]=='>':
                try:
                    val=int(limstr[1:])
                    compdf=compdf[compdf[col]>val]
                    print('Comps with', col,'> ', str(val))
                except:
                    pass
            elif limstr[0]=='<':
                try:
                    val=int(limstr[1:])
                    compdf=compdf[compdf[col]<val]
                    print('Comps with', col,'< ', str(val))
                except:
                    pass
    compdf=compdf[mycols]
    print('\n')
    print(compdf[mycols].to_string(index=False))
    print(compdf[mycols].describe())
    return



EDXcomp=EDXcomp[~EDXcomp['Filename'].str.contains('wide', na=False, case=False)]

def lda(EDXcomp, Elements, n_components=4):
    ''' Perform PCA on compositional data; data is at.% for selected elements
    target is phase, however should pick up unidentified components'''
    np.random.seed(5) # seed generator
    elem=['%'+s for s in Elements]
    EDXcomp=EDXcomp[elem]
    EDXcomp=EDXcomp.dropna() # drop rows if any are nan
    X=EDXcomp[elem].as_matrix()
    y=EDXcomp['Phase'] # target classifiers
    lda = LDA(n_components)
    lda.fit(X)
    # Look at PCA parameters
    
def pca(EDXcomp, Elements, n_components=4):
    ''' Perform PCA on compositional data; data is at.% for selected elements
    target is phase, however should pick up unidentified components'''
    np.random.seed(5) # seed generator
    elem=['%'+s for s in Elements]
    EDXcomp=EDXcomp[elem]
    EDXcomp=EDXcomp.dropna() # drop rows if any are nan
    X=EDXcomp[elem].as_matrix()
    y=EDXcomp['Phase'] # target classifiers
    pca = decomposition.PCA(n_components)
    pca.fit(X)
    # Look at PCA parameters
    
    
    
def getcomp(Ser):
    ''' Pass single compositional row as series, remove trace elems and C/O; 
    create renormalized metals basis composition 
    '''
    dropcols=['%C','%O', '%PtM','%Ga','%Cu'] # drop excluded elements
    thresh=1.0 # threshold for inclusion in at. % 
    mycols=Ser.index.tolist()    
    atcols=[col for col in mycols if '%' in col]
    atcols=[col for col in atcols if col not in dropcols]
    # drop error columns
    atcols=[col for col in atcols if 'err' not in col]
    # drop secondary L lines
    atcols=[col for col in atcols if not col.endswith('L')]
    renorm=0.0
    incl=[]
    for i, col in enumerate(atcols):
        if (Ser[col]>thresh):
            incl.append(col)
            renorm+=Ser[col]
    compstr=''
    for i, col in enumerate(incl):
        compstr+=col[1:]
        val=100*Ser[col]/renorm
        compstr+='%.0f' % val 
        
        print(Ser[col])

# Principal components analysis for compositions


# K-means clustering

 
knn = KNeighborsClassifier()    



def alterbackfit():
    ''' Using known set of elements present (above some threshold) use EDXquantparams
    excluderange column to remove points from background (i.e. for Pt or Ga contamination
    presented ranges are typical for larger time constants (34000,56000) but need to 
    be enlarged for shorter time constants (12400 or less) '''
    

def scatterGUI(axes, fig):
    ''' Interactive selection of points from scatter plots?
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

# working on some way to synchronize EDXcomp identified phases with EDXlog ... 
# is this even necessary or wise? 
colname='Phase'

    
    
# Blob detection with skimage

blobs_log=blob_log(Omap, min_sigma=1,max_sigma=30,threshold=0.1)
blobs_log=pd.DataFrame(blobs_log, columns=['Y','X','size'])
blobs_log['size']=blobs_log['size']*np.sqrt(2) # convert sigma to feature radius




apffilename='JEOLSEMarr3by3at3000overlap20.apf'       

# Plotting all numpy histograms from elementmaps
def plothistograms(elementmaps, Elements):
    
elemlist=[str(el) for el in elements if el not in availableelems]
elemlist=[str(el) for el in elements]

plt.ion()

def makeratioimage

def plotoutliers(df,elem1,elem2):
    ''' Make df scatter plot, fit and find most significant outliers, and plot spectra via report '''
    fig, axes = plt.subplots(nrows=1, ncols=1)
    df.plot.scatter(x=elem1, y=elem2, s=50, color=colorlist[i], ax=axes)
    xcol=df[elem1]
    ycol=df[elem2]
    
    np.polyfit(xcol,ycol, full=True, )
    try:
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xcol, ycol)
    except: # deal with common problems with linregress
        print('Unspecified fitting error')
        return
    xmin=df[elem1].min()
    xmax=df[elem1].max()
    ymin=intercept+slope*xmin
    ymax=intercept+slope*xmax
    axes.plot([xmin, xmax], [ymin, ymax], c='r')
    
    linefunct=lambda x, a, b: a*x+b
    p, cov = curve_fit(linefunct, xcol, ycol)
    
    # residuals
    difference = linefunct(x, *p) - ycol
    plot(x, difference,color='r')
	

def findfitregion(df, fitregion, mandminpts, threshold, fitrange, SEMfileName):
    '''Passing single list of allowable index #s for background fits (no duplicates) 
    remove those with high from list of allowable indices any that show high smoothed-derivatives (i.e. not good for background fitting 
    fitrange and SEMfilename -- error handling only'''
    # loop through Fitregions
    fullptslist=[]
    for i, [fitrange, fitpts, mandminlist, fittype, threshold] in enumerate(Fitregions):
        Backfitdf=df.ix[[x for x in fitpts]] # filter out those not in allowable background ranges
        # these are loaded from SEM_backfit_regions.csv
        Backfitdf=Backfitdf.dropna(subset=['Counts']) # drops above (set to na by ix)
        # now additionally filter out those with derivative above threshold value
        Backfitslice=Backfitdf[(Backfitdf['Savgol']<threshold) & (Backfitdf['Savgol']>-threshold)]
        if Backfitslice.empty==True:
            print('Threshold too low for ', fitrange, ' in ', SEMfileName)
            try:        
                while len(Backfitslice)<4: # reslice until at least 3 points appear for fitting this region
                    threshold=threshold+1 # incrementally raise threshold
                    Backfitslice=Backfitdf[(Backfitdf['Savgol']<threshold) & (Backfitdf['Savgol']>-threshold)]
                print ('Threshold reduced to ', str(threshold))
            except KeyboardInterrupt: # probably not necessary
                print('Threshold while loop interrupted!')
        # TODO add these index #s to fullptslist
        
        # Add in the mandatory minimum fitting points (no threshold applied).. .list of ints
        for i, val in enumerate(mandminlist): 
            if val not in fullptslist:# test if this index # is included in above
                fullptslist.append(val)
    return fullptslist
                
def appendcomment(dfsubset, comment):
    
# Alt version of findfitregion that reduces threshold ... seems problematic
    
def findfitregion(df, fitregion, threshold, fitrange, SEMfileName):
    '''Passing single list of allowable index #s for background fits (no duplicates) 
    remove those with high from list of allowable indices any that show high smoothed-derivatives (i.e. not good for background fitting 
    fitrange and SEMfilename -- error handling only'''
    Backfitdf=df.ix[[x for x in fitregion]] # filter out those not in allowable background ranges
    # these are loaded from SEM_backfit_regions.csv
    Backfitdf=Backfitdf.dropna(subset=['Counts']) # drops above (set to na by ix)
    # now additionally filter out those with derivative above threshold value
    Backfitslice=Backfitdf[(Backfitdf['Savgol']<threshold) & (Backfitdf['Savgol']>-threshold)]
    if Backfitslice.empty==True:
        print('Threshold too low for ', fitrange, ' in ', SEMfileName)
        try:        
            while len(Backfitslice)<4: # reslice until at least 3 points appear for fitting this region
                threshold=threshold+1 # incrementally raise threshold
                Backfitslice=Backfitdf[(Backfitdf['Savgol']<threshold) & (Backfitdf['Savgol']>-threshold)]
            print ('Threshold reduced to ', str(threshold))
        except KeyboardInterrupt: # probably not necessary
            print('Threshold while loop interrupted!')
    return Backfitslice
    
# This is custom version for non-standard energy structure in spectrum (findelemregions is global version)
def findfitregions(SEMfile, Elements, EDXquantparams, logmatch):
    ''' Takes element strings and element list and returns tuple for each elem symbol containing all params 
    needed for finding and quantifying each SEM-EDX peak from given spectrum
    tuple for integ peak is symbol, ideal peak index #, and integ kfactor
    don't apply energy shifts here... apply later when doing integrate''' 
    
    Elemdatamod=[] # returns list of length5 tuples for all elements
    Energyvals = SEMfile.Energy # 
    for i, elem in enumerate(Elements):
        # find row in EDXquantparams for this element
        thiselemdata=EDXquantparams[(EDXquantparams['element']==elem)]
        thiselemdata=thiselemdata.squeeze() # series with this elements params
        
        # integ peak position value is relative to negpeak in smooth-diff (i.e. -5 is 5 eV below ideal negpeak)
        idealev=thiselemdata.energy # ideal energy value of SEM-EDX peak
        
        # TODO find shift most appropriate for element's energy value 
        # find index # for ideal peak position with lambda funct.
        # convert each energy value into index #
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-idealev)) # tuple with index and value of closest energy
        idealindex=temptuple[0] # first of tuple is index #
        peakinrange=temptuple[1]-idealev # should be ~0 if desired peak is in data range
        if abs(peakinrange)<1: # Must skip entire desired element here if it's out of range of the data in this particular spectrum
            kfact=thiselemdata.kfactor # typical sensitivity k-factor associated with element for integration
            errkfact=thiselemdata.errkfact 
            mass=thiselemdata.mass
            # full peak width in keV from EDXquantparams (usually 0.15keV or 15 channels at 0.1eV/chan)
            width=int((thiselemdata.fullwidth*10-1)/2) # integration width in channels for direct integration for this element
            # total # of channels in EDXquantparams but include n-1/2 channels on either side of peak center (usually width is 8 channels)
            
            #Elemdata is a list (of length number of elements) containing length5 tuples
            elemtuple=(elem, idealindex, width, kfact, errkfact, mass) # add tuple with info for this element
            Elemdatamod.append(elemtuple) # now contains proper limits on fitting regions 
        else:
            SEMfileName=logmatch.Filename # logmatch is series
            print('Warning: No quant for ', elem,' for ',SEMfileName, 'data not collected in this energy range.')
    return Elemdatamod

# Need function to compare compositions derived from duplicate spectra 

groupname='Timeconst'
elem1='Fecnts'
elem2='Fe2cnts'
df.plot.scatter(x=elem1, y=elem2, s=50, color=colorlist[i], ax=axes)
from scipy.optimize import curve_fit

def scatterplot(df,elem1, elem2, groupname=''):
    '''Make a scatter plot of one or more datasets of elem1 (x) vs elem2 (y) '''
    fig, axes = plt.subplots(nrows=1, ncols=1)

    grouplist=df[groupname].unique()
    grouplist.sort()
    colorlist=['b','g','r','c','m','y','k']
    for i, group in enumerate(grouplist):
        dfslice=df[df[groupname]==group]
        dfslice.plot.scatter(x=elem1, y=elem2, s=50, color=colorlist[i], ax=axes)
    return
axes.set_yscale('log') 
axes.set_xscale('log')

# Direct scatter plots from these datasets

def scattercompplot(comp1, comp2, elemlist, basis=False):
    '''Pass two versions of composition calculation (using different lines or whatever) and compare 
    major elements using scatter graphs .. single point for each sample
    uses inner merge to select only subset with values from each df 
    basis=False means use the atomic percent columns %Fe, %Mg, etc.; otherwise using Fe, Mg (adjusted counts)'''
    comp1=C2010Woldquant    
    comp2=C2010WEDXcomp
    elemlist=Elements
    
    elemlist=[re.match('\D+',i).group(0) for i in elemlist] 
    # strip number from peaks like Fe2 if present; columns will be element names (Fe) not peak names (Fe2)
    if basis==False: # use atomic % (which is the default), not basis for each element
        elemlist=['%'+s for s in elemlist]
    numareas=len(elemlist)
    
    # set nrows and ncols for figure of proper size
    cols=divmod(numareas,2)[0]+ divmod(numareas,2)[1]
    if numareas>1:
        rows=2
    else:
        rows=1        
    fig, axes = plt.subplots(nrows=rows, ncols=cols) # axes is array
    # merge dfs with comp1 and comp2 using inner join    
    df=pd.merge(comp1, comp2, how='inner', on=['Sample'], suffixes=('','b'))    
    for i,elem in enumerate(elemlist):
         # determine which subplot to use
        if (i+1)%2==1:
            rownum=0
        else:
            rownum=1
        colnum=int((i+1)/2.1)       
        xcol=elem
        ycol=elem+'b' # same element from second dataset
        if numareas==1: # deal with single subplot separately
            df.plot.scatter(x=xcol, y=ycol, ax=axes) # single plot axes has no [#,#]
        else:
            df.plot.scatter(x=xcol, y=ycol, ax=axes[rownum,colnum])
            # linear regression: fitting, plot and add labels
        data1=df[elem]
        colname=elem+'b'
        data2=df[colname]
        # slope,intercept=np.polyfit(data1, data2, 1)  numpy version
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(data1, data2)
        # set x range for linear plot 
        text1=str(round(slope,2))+' *x +' + str(round(intercept,2))
        text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,3))
        xmax=max(max(data1),max(data2))*1.1 # set to slightly larger than max of dataset
        x=np.linspace(0,xmax,100) # setting range for 
        if numareas==1: # deal with single subplot separately
            axes.text(0.025,0.9, text1, fontsize=12, transform=axes.transAxes)
            axes.text(0.025,0.8, text2, fontsize=12, transform=axes.transAxes)
            plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        else: # typical multiarea plot
            axes[rownum,colnum].text(0.025,0.9, text1, fontsize=12, transform=axes[rownum,colnum].transAxes)
            axes[rownum,colnum].text(0.025,0.8, text2, fontsize=12, transform=axes[rownum,colnum].transAxes)
            plt.axes(axes[rownum,colnum]) # set correct axes as active
            plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
    return