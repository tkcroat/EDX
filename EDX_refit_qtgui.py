# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 07:45:29 2017

@author: tkc
"""
import os, glob
from PyQt5 import QtGui, QtCore, QtWidgets
import numpy as np 
import pandas as pd

from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib as mpl
mpl.use('Qt5Agg') # choose option for interactive mpl 
from matplotlib.widgets import RectangleSelector

# TODO ... roll EDXcanvas into EDX_app class?
class EDX_app(QtWidgets.QMainWindow):
    ''' Main app/GUI for SEM-EDX background fit corrections 
    old version with EDXcanvas split from EDX_app '''
    def __init__(self, parent=None):
        super(EDX_app, self).__init__(parent)
        # load UI file
        # uic.loadUi('gui\\ui\\pyLATTICE_GUI.ui')

        self.setGeometry(50,50,1050,650) # main window sizing
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("EDX refitting")

        self.main_widget = QtWidgets.QWidget(self) # define main matplotlib widget                
        self.vbl = QtWidgets.QVBoxLayout(self.main_widget) # simple layout w widgets stacked

        self.main_widget.setFocus() # keyboard input to this widgets
        self.setCentralWidget(self.main_widget)

        self.EDXdata=EDXdataset() # instance w/ 6 pandas dataframes with params

        # sets EDX canvas (mpl plot and spinbox) as main widget
        sc = EDXCanvas(self, self.main_widget, self.EDXdata)
        self.vbl.addWidget(sc)
        self.setLayout(self.vbl) # finalizes layout 
        self.show()
        
    def file_Quit(self):
        ''' Done button within child (EDXCanvas) calls quit for main window 
        '''
        self.close()

class EDXCanvas(QtWidgets.QWidget):
    """Qwidget with mpl window and assorted buttons for EDX refitting procedure
    (not main window which is separate)
    EDXdataset contains the list of plottable files
    """
    # passed EDXdata as argument
    def __init__(self, parent, main_widget, EDXdata, **kwargs):     
        # super(EDXCanvas, self).__init(parent)
        # still needs call of super_class init
        QtWidgets.QWidget.__init__(self)
        self.parent=parent
        # six pandas dataframes passed as arguments
        self.EDXdata=EDXdata
        # Instantiate matplotlib plot (declares it as figurecanvas)
        self.canvas=MplCanvas() # instantiate matplotlib canvas
        # Standard mpl toolbar for zooming, panning
        self.mpl_toolbar = NavigationToolbar(self.canvas, self) # needs a parent frame
        
        # Basic vertical box layout 
        self.vbl = QtWidgets.QVBoxLayout(self)
        self.vbl.addWidget(self.mpl_toolbar)
        self.vbl.addWidget(self.canvas) # add main mpl plot 
        
        self.vbl.addWidget(QtWidgets.QLabel('Spectral filenumber index'))
        self.spectralspin = QtWidgets.QSpinBox()
        self.spectralspin.setRange(0, self.EDXdata.numfiles-1)
        self.vbl.addWidget(self.spectralspin)
        
        '''
        # set up for mpl pick events 
        self.canvas.mpl_connect('pick_event', self.on_pick()) 
        ''' 
        print("Running plot_init method")        
        self.plot_init() # for plotting first spectrum and background (works!)
        # QtCore.QObject.connect(self.spinBox_spacegroup, QtCore.SIGNAL(_fromUtf8("valueChanged(int)")), self.plot_number)
        self.spectralspin.valueChanged.connect(self.loadplot_spe) # connect to plot_spe function 
        
        # Save new background
        self.removeptbut= QtWidgets.QPushButton("Remove background pts")
        self.removeptbut.clicked.connect(self.remove_backpts)
        self.vbl.addWidget(self.removeptbut)
        self.refitbut= QtWidgets.QPushButton("Refit backgrounds")
        self.refitbut.clicked.connect(self.refit_backgrounds)
        self.vbl.addWidget(self.refitbut)
        self.savebut= QtWidgets.QPushButton("Save modified background")
        self.savebut.clicked.connect(self.save_changes)
        self.vbl.addWidget(self.savebut)
        
        '''# Listbox for holding evranges (pyqtWrapperType)
        # 9/21/17 now stored (and altered) in backlists attribute 
        self.evrangelist=evrangelist(EDXdata,self.spectralspin.value())
        self.vbl.addWidget(self.evrangelist)
        print('Evrangelist of type', type(evrangelist))
        # TODO handle double-click selection of a given fitrange
        self.evrangelist.currentItemChanged.connect(self.plotbackpts)
        '''

        # done button connected to filequit in AppWindow
        self.myb = QtWidgets.QPushButton("Done")
        # connect button to file_quit (old qt4 style)        
        self.myb.clicked.connect(self.parent.file_Quit)
        
        self.vbl.addWidget(self.myb)

        self.setLayout(self.vbl) # finalizes layout 

    def remove_backpts(self):
        ''' Find subset of index #s in removal ev range '''
        # Gets currently visible subset of points
        xmin, xmax=self.canvas.axes.get_xlim()
        ymin, ymax=self.canvas.axes.get_ylim()
        # find index #s of EDX data within current xlim and ylim
        EDXslice=self.EDXfile[ (self.EDXfile['Energy']>xmin) & (self.EDXfile['Energy']<xmax) ]
        EDXslice=EDXslice[ (EDXslice['Counts']>ymin) & (EDXslice['Counts']<ymax) ]
        indrange=EDXslice.index.unique().tolist()
        currvals=self.backlists
        # now remove said points from backlists (list of lists)
        self.backlists = [[x for x in l if x not in indrange] for l in self.backlists]
        for i, thislist in enumerate(self.backlists):
            if thislist!=currvals[i]:
                self.refitflags[i]=True # mark change to guide refitting
        # Replot after removing offensive points 
        self.replot_spe()

    def alter_back(self, fitparams, minindex, maxindex):
        ''' Change EDXfile with new fit params over single subregion '''
        for index in range(minindex,maxindex+1):
            xval=self.EDXfile.loc[index]['Energy']
            self.EDXfile=self.EDXfile.set_value(index, 'Backfit', fitparams[0] * xval**2+ fitparams[1] * xval + fitparams[2])
        return self.EDXfile
    
    def refit_backgrounds(self):
        ''' Refit altered background points regions (modified fitbackgrounds function)
        backlists holds index #s (after interactive modification) for points to include in fit
        Background fit for each direct peak(opens source spectrum as EDXfile, 
        fits peak backgrounds above and below using Elemdata, saves background to source csv (overwrites existing fits), 
        also saves linear fit params to logdataframe with position/amplitude/etc;
        
        '''
        for i, indnums in enumerate(self.backlists):
            if self.refitflags[i]: # skip if not altered
                print('Refitting region', str(i),'from indices', str(indnums[0]),' to ', str(indnums[1]))
                thisfitreg=self.EDXfile[self.EDXfile.index.isin(indnums)]
                # Force low energy counts to zero
                for index,row in thisfitreg.iterrows():
                    if index < 5:
                        thisfitreg=thisfitreg.set_value(index,'Counts',0)
                # A, B, C fit params of parabolic fit
                fitparams = self.fitparabola(thisfitreg)
                # Update current value in fitparamslists
                self.fitparamslists[i]=[fitparams[0], fitparams[1],fitparams[2], np.nan, 'parabola']
                # recreate background over specific refit region
                self.EDXfile=self.alter_back(fitparams, indnums[0], indnums[-1])
        self.replot_spe()
        # TODO add returnoverlaps and crossfadeback for smooth transitions between backfit regions?

    def fitparabola(df):
        '''Pass chunk of EDX dataframe and perform polynomial/parabola fit
        return chunk with altered backfit column
        '''
        xcol=df['Energy']
        ycol=df['Counts']
        try:
            fitparams=np.polyfit(xcol, ycol, 2) # A, B and C params of 2nd order poly fit
        except: # deal with common problems with linregress
            print('Fitting error from ', "{0:.2f}".format(df.Energy.min()),'to ',"{0:.2f}".format(df.Energy.max()))
            print('df is empty =', df.empty) # indicate problem with empty frame probably due to thresholds
            fitparams=('n/a','n/a','n/a') # return all n/a
            return ['n/a','n/a','n/a']
        return fitparams

    def save_changes(self):
        ''' on button, save changes to backpts regions, and altered background itself for single EDXfile  '''
        # Save modified selected backfit points, new params, and EDXfile itself
        # save modified EDX file
        self.EDXfile.to_csv(self.filename, index=False)
        # Save modified backpoints and fit params
        #TODO have these already been updated??
        matches=self.EDXdata.Backfitlog[self.EDXdata.Backfitlog['Filename']==self.filename]
        if len(matches)!=len(self.backlists):
            print("Couldn't find matching backfitlog entries.")
            return
        for i in range(0, len(self.backlists)+1):
            if self.refitflags[i]:
                thisind=matches.index[i] # get index of ith row
                tempstr=",".join(self.backlists[i])
                self.EDXdata.Backfitlog=self.EDXdata.Backfitlog.set_value(thisind, 'Backfitpts', tempstr)
                self.EDXdata.Backfitlog=self.EDXdata.Backfitlog.set_value(thisind, 'A', self.fitparamslists[i][0])
                self.EDXdata.Backfitlog=self.EDXdata.Backfitlog.set_value(thisind, 'B', self.fitparamslists[i][1])
                self.EDXdata.Backfitlog=self.EDXdata.Backfitlog.set_value(thisind, 'C', self.fitparamslists[i][2])
                self.EDXdata.Backfitlog=self.EDXdata.Backfitlog.set_value(thisind, 'D', self.fitparamslists[i][3])
                self.EDXdata.Backfitlog=self.EDXdata.Backfitlog.set_value(thisind, 'Fittype', 'parabola')
        self.EDXdata.save() # use class method 
        print("Backfitparamslog saved")
    
    def openspefile(self, filename):
        ''' Open and returns named Augerfile and subset of background points from same '''
        print("Running open spe file")
        filename=filename.replace('.emsa','.csv')
        print('Filename at openspe is ', filename)
        if os.path.exists(filename):
            self.EDXfile=pd.read_csv(filename, encoding='cp437')
            self.filename=filename
            print(filename,' loaded.')
        elif os.path.exists('sub//'+filename):
            self.EDXfile=pd.read_csv('sub//'+filename, encoding='cp437')
            self.filename='sub//'+filename
            print(filename,' loaded.')
        else: # If not found skip to next file
            print(filename, ' not found.')
            self.EDXfile=pd.DataFrame() # Return empty frame
            self.filename=''
        
        matches=self.EDXdata.Backfitlog[self.EDXdata.Backfitlog['Filename']==filename]
        self.backlists=[] # list of backpoint index #s for each backfit region
        self.fitparamlists=[]
        for index, row in matches.iterrows():
            templist=row.Backfitpts.replace('[','').replace(']','').split(',')
            templist=[int(i) for i in templist]
            self.backlists.append(templist)
            self.fitparamslists.append([row.A, row.B, row.C, row.D, row.Fittype])

        # Change these flags to indicate regions for refitting
        self.refitflags=[False]*len(self.backlists)
        print('Loaded ', filename,' of type ',type(self.EDXfile))
        return self.EDXfile, self.backlists, self.fitparamslists, self.refitflags

    def loadplot_spe(self):
        ''' Load and plot spectral file associated with given index number 
        why was current previous needed '''
        mypath=os.getcwd()
        index = self.spectralspin.value()
        filename=mypath+'\\'+str(self.EDXdata.EDXlog.iloc[index]['Filename'])
        filename=filename.replace('.emsa','.csv')
        try:            
            self.EDXfile, self.backlists, self.refitflags=self.openspefile(filename)
            # backlists is list of index #s of background pts for each region 
            self.canvas.axes.cla() # clear existing plots
            self.EDXfile.plot(x='Energy', y='Counts', ax=self.canvas.axes)
            self.EDXfile.plot(x='Energy', y='Backfit', ax=self.canvas.axes)
            self.canvas.draw() # update plot
            # Also plot pts used for background fits (compress to single list)
            fullbackpts=[item for sublist in self.backlists for item in sublist]
            backsubset=self.EDXfile[self.EDXfile.index.isin(fullbackpts)]
            backsubset.plot.scatter(x='Energy', y='Counts', ax=self.canvas.axes)  
            self.canvas.draw() # update plot
        except:
            print('failed load')
        # TODO deal with multiple areas of same filename

    def replot_spe(self):
        ''' Replot spe and backfitpts after changes made '''
        try:            
            # backlists is list of index #s of background pts for each region 
            self.canvas.axes.cla() # clear existing plots
            self.EDXfile.plot(x='Energy', y='Counts', ax=self.canvas.axes)
            self.EDXfile.plot(x='Energy', y='Backfit', ax=self.canvas.axes)
            self.canvas.draw() # update plot
            # Also plot pts used for background fits (compress to single list)
            fullbackpts=[item for sublist in self.backlists for item in sublist]
            backsubset=self.EDXfile[self.EDXfile.index.isin(fullbackpts)]
            backsubset.plot.scatter(x='Energy', y='Counts', ax=self.canvas.axes)  
            self.canvas.draw() # update plot
        except:
            print('Failed replotting')
        # TODO deal with multiple areas of same filename
        
    def plot_init(self):
        ''' Plot first EDX file in stack along w/ background pts 
        .. seems to be working correctly '''
        mypath=os.getcwd()
        print("File path is", mypath)
        filename=mypath+'\\'+str(self.EDXdata.EDXlog.iloc[0]['Filename'])
        filename=filename.replace('.emsa','.csv')
        #filename='c:\\temp\\EDX\\area1hFeS1.csv'
        self.EDXfile=self.openspefile(filename)
        try:
            pass
            # Somehow this is causing problems?
            # self.EDXfile.plot(x='Energy', y='Counts', ax=self.canvas.axes)
            # self.EDXfile.plot(x='Energy', y='Backfit', ax=self.canvas.axes)
        except:
            print('Problem plotting initial spe file')
        
    # legacy method ... Remove after working 
    def plotbackpts_old(self, current, previous):
        ''' Plots background points when ev range clicked in list
        current list item (and previous one) sent as args by currentItemChanged '''
        index=self.spectralspin.value()
        # print(index)
        thisfile=self.EDXdata.EDXlog.iloc[index]['Filename']
        # print(current.text())
        thisrange=current.text() # energy range passed from evrangelist
        # Get sequence of background fits for this file
        match=self.EDXdata.Backfitlog[(self.EDXdata.Backfitlog['Filename']==thisfile) & (self.EDXdata.Backfitlog['Fitrange']==thisrange)]
        backpts=match.iloc[0]['Backfitpts'].replace('[','').replace(']','')
        backpts=backpts.replace('[','').replace(']','')
        backpts=backpts.split(',')
        backpts=[s.strip() for s in backpts]
        backpts=[int(i) for i in backpts] # final set of index values
        # select this subset of pts by index # from plotted spectrum
        self.EDXfile=self.openspefile(thisfile)
        # Plot this point subset
        backsubset=self.EDXfile[self.EDXfile.index.isin(backpts)]
        backsubset.plot.scatter(x='Energy', y='Counts', ax=self.canvas.axes)
        # now call interactive 
        # LassoMananger(self.canvas.axes, backsubset)
        self.canvas.draw() # update plot
        
    def file_quit(self):
        ''' Done button within child (EDXCanvas) calls quit for main window 
        '''
        self.close()
        
class MplCanvas(FigureCanvas):
    """ Standard matplotlib canvas / also Qwidget  """

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        # set parent?? 

class EDXfile():
    ''' Loads all info related to active EDX file '''
    def __init__(self, *args, **kwargs):
        # open files 
        self.EDXlog, self.Backfitlog, self.Integlog, self.Peakfitlog, self.EDXquantparams, self.Interferences=self.open_main_files()
        self.filelist=np.ndarray.tolist(self.EDXlog.Filenumber.unique())
        self.numfiles=len(self.EDXlog)
        self.EDXfile=
        print(str(self.numfiles),' loaded from EDXdataset.')
    
    def save(self):
    		''' Save backfitlog files  '''
    		fname=self.path+'\\Backfitparamslog.csv'
    		self.Backfitlog.to_csv(fname, index=False)
    		'''
    		fname=self.path+'\\EDXparamlog.csv'
    		self.EDXlog.to_csv(fname, index=False)
    		fname=self.path+'\\Peakfitlog.csv'
    		self.Peakfitlog.to_csv(fname, index=False)
    		fname=self.path+'\\Integquantlog.csv'
    		self.Integlog.to_csv(fname, index=False)
    		'''

class EDXdataset():
    ''' loads all dataframes with EDX parameters from current project folder '''
    def __init__(self, *args, **kwargs):
        # open files 
        self.EDXlog, self.Backfitlog, self.Integlog, self.Peakfitlog, self.EDXquantparams, self.Interferences=self.open_main_files()
        self.filelist=np.ndarray.tolist(self.EDXlog.Filenumber.unique())
        self.numfiles=len(self.EDXlog)
        # Autoload first file
        self.EDXfile=EDXfile()
        print(str(self.numfiles),' loaded from EDXdataset.')
        
    def change_EDXfile(self):
        ''' When index from spinbox changes load/instantiate different EDXfile '''
        pass 
    
    def open_main_files(self):
        ''' Auto loads EDX param files from working directory including 
        EDXparalog- assorted params associated w/ each SEM-EDX or TEM-EDX emsa file 
        Backfitparamslog - ranges and parameters for EDX background fits
        Integquantlog - subtracted and corrected counts for chosen elements
        Peakfitlog - params of gaussian fits to each element (xc, width, peakarea, Y0, rsquared)'''
        if os.path.exists('EDXparamlog.csv'):
            EDXlog=pd.read_csv('EDXparamlog.csv', encoding='cp437')
            start=len(EDXlog)
            EDXlog['Comments']=EDXlog['Comments'].replace(np.nan,'')
            EDXlog=EDXlog[~EDXlog['Comments'].str.contains("exclude",na=False, case=False)]
            if start-len(EDXlog)!=0:
                print('Dropped',str(int(start-len(EDXlog))), 'excluded spectral files.')
            self.path=os.getcwd()
        else:
            files=glob.glob('*paramlog.csv')
            if len(files)==1:
                print('Loaded params file', files[0])
                EDXlog=pd.read_csv(files[0], encoding='cp437')
            else:
                print("Couldn't find EDX params file in existing folder.")
                EDXlog=pd.DataFrame() # load blanks to avoid error but cd probably needed
                Integlog=pd.DataFrame()
                Backfitlog=pd.DataFrame()
                Peakfitlog=pd.DataFrame()
        if os.path.exists('Peakfitlog.csv'):
            Peakfitlog=pd.read_csv('Peakfitlog.csv', encoding='cp437')
        else:
            Peakfitlog=pd.DataFrame()
        if os.path.exists('Backfitparamslog.csv'):
            Backfitlog=pd.read_csv('Backfitparamslog.csv', encoding='cp437')
        else:
            Backfitlog=pd.DataFrame()
        if os.path.exists('Integquantlog.csv'):
            Integlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
        else:
            Integlog=pd.DataFrame()
        # Print TEM or SEM to console based on beam kV
        if EDXlog['Beamkv'].max()>30:
            print(EDXlog['Beamkv'].max(),'keV TEM spectra loaded.')
            EDXquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\TEMquantparams.csv', encoding='utf-8')
            Interferences=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\TEM_interferences.csv', encoding='utf-8')
        else:
            print(EDXlog['Beamkv'].max(),'keV SEM spectra and quant params loaded.')
            EDXquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEMquantparams.csv', encoding='utf-8')
            Interferences=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEM_interferences.csv', encoding='utf-8')
        return EDXlog, Backfitlog, Integlog, Peakfitlog, EDXquantparams, Interferences
    
    def save(self):
    		''' Save backfitlog files  '''
    		fname=self.path+'\\Backfitparamslog.csv'
    		self.Backfitlog.to_csv(fname, index=False)
    		'''
    		fname=self.path+'\\EDXparamlog.csv'
    		self.EDXlog.to_csv(fname, index=False)
    		fname=self.path+'\\Peakfitlog.csv'
    		self.Peakfitlog.to_csv(fname, index=False)
    		fname=self.path+'\\Integquantlog.csv'
    		self.Integlog.to_csv(fname, index=False)
    		'''
		
class evrangelist(QtWidgets.QListWidget):
    ''' QListWidget holding list of ev ranges used during fitting 
    QListView is another more expansive option  '''
    def __init__(self, EDXdata, index, parent=None):
        ''' populate list with evranges for this EDXdataset (use first row upon init)'''
        QtWidgets.QListWidget.__init__(self, parent)

        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection) # allows multiple selections
        self.EDXdata=EDXdata
        self.index=index
        self.populate(EDXdata, index)
        # self.setranges(self.EDXdata, self.index)
        # Why is parent of None needed here?
    
    def _ranges(self, EDXdata, index):
        ''' Actual list of ev ranges.. index is current value of spectral spin'''
        thisfile=self.EDXdata.EDXlog.iloc[0]['Filename']
        match=self.EDXdata.Backfitlog[self.EDXdata.Backfitlog['Filename']==thisfile]
        evranges=np.ndarray.tolist(match.Fitrange.unique())
        # self._populate()
        print('ev ranges initialized')
        for i, val in enumerate(evranges):
            print(val)
        return evranges
        
    def populate(self, EDXdata, index):
        ''' Initialization or repopulate of fit ranges as items to this list '''
        self.clear()
        for evrange in self._ranges(EDXdata,index):
            item = QtWidgets.QListWidgetItem(evrange)
            self.addItem(item) # add item to QListWidget
    # Need to make item change in list find and plot background points