# -*- coding: utf-8 -*-
"""
Interactive EDX background refitter 
Created on Wed Oct 11 00:44:29 2017

@author: tkc
"""
import sys
import numpy as np
import tkinter as tk
import os
import tkinter.messagebox as tkmess
from tkinter import filedialog
import matplotlib as mpl # using path, figure, rcParams
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.widgets import Lasso
from matplotlib import path
# import defined data classes 
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX')
from EDX_data_classes import EDXfile, EDXdataset

PLOT_SIZE = (10,6) # 8, 5 or 
MPL_STYLE = {
    "text.color":"k",
    "axes.labelcolor":"black",
    "axes.edgecolor":"0.4",
    "axes.facecolor":"white",   
    "xtick.color": "lightblue",
    "ytick.color": "lightblue",
    "figure.facecolor":"white",
    "figure.edgecolor":"white",
    "text.usetex":False
}
mpl.rcParams.update(MPL_STYLE)

#-------------------Misc.---------------------#

def launch_refitter():
    ''' Launcher function for tk refitter GUI '''
    root = tk.Tk()
    root.wm_title("EDX refitter")
    screensize=[root.winfo_screenwidth(),root.winfo_screenheight()]
    w, h=[int(1.0*i) for i in screensize]
    h-=100 # Taskbar
    x,y=[int(0.00*i) for i in screensize]
    root.geometry('%dx%d+%d+%d' % (w,h,x,y))
    # choose EDX data directory (starting at current)
    currdir=filedialog.askdirectory(initialdir =os.getcwd(), 
        title='Select EDX data directory')
    plotter = GUIMain(root, currdir)
    root.mainloop()
    return

class GUIMain():
    ''' Main container for plotter, options (at right), and fileloader (bottom) 
    pass current working directory as default directory'''
    def __init__(self,root, currdir):
        self.root = root
        self.root.wm_title("EDX refitter ")
        self.top_frame = tk.Frame(self.root)
        self.top_frame.pack(side=tk.TOP) 
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM)
        self.plot_frame = tk.Frame(self.top_frame)
        self.plot_frame.pack(side=tk.LEFT)
        self.refit_frame = tk.Frame(self.top_frame) 
        self.refit_frame .pack(side=tk.LEFT)
        self.loader_frame = tk.Frame(self.bottom_frame)
        self.loader_frame.pack(side=tk.LEFT,fill=tk.BOTH)
        self.plotter  = GUIPlotter(self.plot_frame,self)
        self.refitter  = GUIRefitter(self.refit_frame,self)
        self.loader  = GUIprojectloader(self.loader_frame,self, currdir)

class NavSelectToolbar(NavigationToolbar2TkAgg):
    ''' Custom matplotlib toolbar w/ lasso pt remover and point picker
    parent is GUIplotter
    '''
    def __init__(self, canvas, root, parent):
        self.canvas = canvas
        self.root   = root
        self.parent = parent # plotter is parent
        self.ax= self.parent.ax # axes needed for interaction
        self.xys = None # for xy vals later associated with plot
        self.select = None # lasso selected points for removal

        # Generic mpl toolbar using tkagg (with standard buttons)
        NavigationToolbar2TkAgg.__init__(self, canvas,root)
        # Create lasso and link to multi_select_callback
        self.lasso_button= tk.Button(master=self, text='lasso', padx=2, pady=2, 
                                     command=self.startlasso)
        self.lasso_button.pack(side=tk.LEFT,fill="y")
        self.remove_pts_button= tk.Button(master=self, text='Remove pts', 
                                          padx=2, pady=2, command=self.removepts)
        self.remove_pts_button.pack(side=tk.LEFT,fill="y")
        self.picker_button= tk.Button(master=self, text='add point', padx=2, 
                                      pady=2, command=self.startpicker)
        self.picker_button.pack(side=tk.LEFT,fill="y")
        self.show_button= tk.Button(master=self, text='Show backfit segments', 
                                    padx=2, pady=2, command=self.showbackseg)
        self.show_button.pack(side=tk.LEFT,fill="y")
        print('toolbar loaded')
        # temp definition of pick_button (invoked in GUIplotter)

    def startlasso(self):
        ''' Activated by lasso menu bar button on click; disconnects prior IDs, prep for lasso button press
        '''
        print('startlasso called')
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpresslasso)
        print('end of startlasso')

    def onpresslasso(self, event):
        ''' Create lasso when button pressed on active canvas/axes '''
        # ensure that current dataset is active
        print('onpress lasso called')
        self.xys = self.parent.xy # passed from plotter (parent)
        print('Length of xys is', len(self.xys))
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callbacklasso)
        #  self.canvas.widgetlock(self.lasso)  # skip... gives ValueError: already locked
        print('end of onpress lasso')

    def callbacklasso(self, verts):
        print('callback called')
        p = path.Path(verts)
        # true/false array
        ind = p.contains_points(self.xys)
        self.selected=[i for i in range(0, len(self.xys)) if ind[i]==True]
        print('Selected points are:', self.selected)
        self.canvas.draw_idle()
        # self.canvas.widgetlock.release(self.lasso) # valueerror you don't own this lock
        del self.lasso
        self.canvas.mpl_disconnect(self.cid) # disconnect lasso tool
        print('finished with callback')
        
    def startpicker(self):
        ''' Activated by lasso menu bar button on click; disconnects prior IDs, prep for lasso button press
        '''
        print('startpicker called')
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpresspick)
        print('end of startlasso')

    def onpresspick(self, event):
        ''' Make picker connection for adding pointsGet closest point in spectrum and 
        add to background points for refitting (in plotter)'''
        # just need event.xdata and ydata
        print('onpresspick called')
        print('X/y is', event.xdata, event.ydata)
        self.parent.point_add_callback(event.xdata, event.ydata)
        self.canvas.mpl_disconnect(self.cid)

    def showbackseg(self):
        ''' adds separate background fitted segments to plot in different colors 
        '''
        print('showbackseg called')
        # Shows current fit values over different ranges in different colors 
        self.parent.showfitsegments()

    def removepts(self):
        ''' Remove points currently in index of active lman 
        '''
        print('remove pts called')
        if self.selected == None:
            print('No active lassoed points')
            return
        print('Chosen indices are', self.selected)
        # Call point removal method 
        self.parent.points_removed_callback(self.selected)
        print('end of removepts call')

class GUIPlotter():
    def __init__(self,root, parent):
        self.root = root
        self.parent = parent
        self.xy = None # used by lasso selector (init below in plot_backfitpts)
        self.backsubset = None 
        self.figure = mpl.figure.Figure(figsize=PLOT_SIZE, dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
        self.canvas = FigureCanvasTkAgg(self.figure,self.root)
        # Custom navselecttoolbar w/ interactive buttons
        self.toolbar = NavSelectToolbar(self.canvas,self.root,self)
        self.toolbar.update()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.EDXfile = None
        self.canvas.show()

    def associate_EDXfile(self, EDXfile):
        ''' Associate EDX file created/loaded in GUIrefitter with GUIplotter 
        called from GUIrefitter '''
        # Is this a separate instance from that in guiopts?
        # print('New EDXfile associated with plot window')
        self.EDXfile = EDXfile
        self.plot() # 

    def plot(self,**kwargs):
        if self.EDXfile is None:return
        self.ax.cla() # clear axes
        self.current_plot = self.ax.plot(self.EDXfile.EDXdf['Energy'], 
            self.EDXfile.EDXdf['Counts'], color='k', picker=True,**kwargs)
        # add existing background fit in red 
        self.ax.plot(self.EDXfile.EDXdf['Energy'], self.EDXfile.EDXdf['Backfit'], 
                     color='r', picker=True,**kwargs)
        # Now plot backfitpts as scatter
        self.plot_backfitpts()

        self.canvas.show()
    
    def plot_backfitpts(self):
        ''' Get subset of points used for background fits '''
        # Make big list of all energy vals (in eV) used across all fit regions
        bpts=[]
        for i, ptlist in enumerate(self.EDXfile.backfitpts):
            bpts.extend(ptlist)

        self.backsubset=self.EDXfile.EDXdf[self.EDXfile.EDXdf.index.isin(bpts)]
        # plot backfit subset from filtered EDXdf dataframe
        self.backsubset.plot.scatter(x='Energy', y='Counts', color='b', ax=self.ax)

        # Initialize xy vals for use with lasso or selector (list of x,y tuples)
        self.xy=[]
        # print('length of backsubset is', len(self.backsubset))
        for index, row in self.backsubset.iterrows():
            self.xy.append((row.Energy, row.Counts))
            
    def point_add_callback(self, xpoint, ypoint):
        ''' Linked to pick_button in NavSelectToolbar.. find nearest spectral datapoint and 
        add to backfit points list'''
        if self.EDXfile is None:return
        print('point_add_callback called')
        # Remove this single value from every list in backptslist
        self.parent.refitter.add_pts(xpoint, ypoint)
        
    def points_removed_callback(self, inds):
        ''' Linked to lasso_button in NavSelectToolbar.. can launch GUI multiviewer 
        inds are lasso selected points (backpts xys in same order as original) '''
        #TODO fix problem with lasso selector that won't die
        if self.EDXfile is None:
            return
        # removepts=self.xy[]
        print('lasso callback called')
        badxvals=[]
        for i, index in enumerate(inds):
            badxvals.append(self.xy[index][0])
        print('Bad xvals are:', ",".join([str(i) for i in badxvals]))
        # background points stored as index numbers (starting at 0)
        # conversion is xval (in keV) = indexnum*.01 +0.01)
        badslice=self.backsubset[self.backsubset['Energy'].isin(badxvals)]
        badinds=badslice.index.tolist() # Need to remove these from 
        # TODO Use energy vals or use index numbers
        print('Bad ind #s are', ','.join([str(i) for i in badinds]))
        self.parent.refitter.remove_badpts(badinds)

    def showfitsegments(self):
        ''' Plot separate fitted segments in different colors on plot over appropriate ranges '''
        #TODO fix problem with lasso selector that won't die
        if self.EDXfile is None:
            return
        print('Running showfitsegments')
        colorlist=['g','c','m','olive','pink','purple']
        for i, ft in enumerate(self.EDXfile.fitorders):
            A=self.EDXfile.backfitparams[i][0]
            B=self.EDXfile.backfitparams[i][1]
            C=self.EDXfile.backfitparams[i][2]
            D=self.EDXfile.backfitparams[i][3]
            xvals=np.arange(min(self.EDXfile.backfitpts[i])/100, max(self.EDXfile.backfitpts[i])/100, 0.1)
            if ft=='linear':
                print('plotting linear')
                self.ax.plot(xvals, A*xvals+B, color=colorlist[i])
            elif ft=='parabola':
                self.ax.plot(xvals, A*xvals**2+B*xvals+C, color=colorlist[i])            
            elif ft=='cubic':
                self.ax.plot(xvals, A*xvals**3+B*xvals**2+C*xvals+D, color=colorlist[i])
        self.canvas.show() # Now show these lines
        
class GUIprojectloader():
    ''' Picks directory and loads main Auger param files
    needs current path (which should be set to working data directory) '''
    def __init__(self,root, parent, currdir):
        self.root = root
        self.parent = parent # GUImain is parent
        self.top_frame = tk.Frame(self.root)
        self.top_frame.pack(side=tk.TOP,anchor=tk.W)
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM,fill=tk.BOTH,expand=1)
        tk.Label(self.top_frame,text="Directory:",padx=8,pady=2,
                 height=1).pack(side=tk.LEFT,anchor=tk.W)
        self.directory_entry = tk.Entry(self.top_frame,width=90,bg="lightblue",
                fg="black",highlightcolor="lightblue",insertbackground="black",
                highlightthickness=2)
        self.directory_entry.pack(side=tk.LEFT,fill=tk.BOTH,expand=1,anchor=tk.W)
        self.directory_entry.insert(0,currdir)
        tk.Button(self.top_frame,text="Browse",command=self.launch_dir_finder
                  ).pack(side=tk.LEFT,fill=tk.BOTH,expand=1,anchor=tk.W)
        self.load_button = tk.Button(self.bottom_frame,text="Load EDX project folder",
            width=60, command=self.load_EDXdataset)
        self.load_button.pack(side=tk.BOTTOM,expand=1,anchor=tk.CENTER)
        self.autoload_EDXdataset(currdir)

    def autoload_EDXdataset(self, currdir):
        ''' Autoload directory chosen in launcher 
        '''
        EDXdata = EDXdataset(currdir)
        # pass to GUIrefitter and set spectral selector spinbox values
        self.parent.refitter.associate_EDXdataset(EDXdata)        
        
    def load_EDXdataset(self):
        ''' Load standard AES files (paramlogwith data returned to a DataManager '''
        directory = self.directory_entry.get()
        EDXdata = EDXdataset(directory)
        # pass to GUIrefitter and set spectral selector spinbox values
        self.parent.refitter.associate_EDXdataset(EDXdata)        
        # TODO associate EDXdataset with GUIplotter (or just selected EDXfile)
    
    def launch_dir_finder(self):
        directory = filedialog.askdirectory()
        self.directory_entry.delete(0,tk.END)
        self.directory_entry.insert(0,directory)

class GUIRefitter():
    ''' Parent is GUImain, manages EDXfile displayed in GUIplotter
    handles addition/removal of points for background (re)fitting'''
    def __init__(self,root,parent):
        self.root = root
        self.parent = parent
        self.EDXdataset = None  # created in GUIprojectloader but associated here
        # Instance of EDXfile local to the refitter
        self.EDXfile = None
        self.specselect_frame = tk.Frame(self.root,pady=10)
        self.specselect_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        self.currfile_frame = tk.Frame(self.root,pady=10)
        self.currfile_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        self.misc_opts_frame = tk.Frame(self.root,pady=10)
        self.misc_opts_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Frame for background fit ev ranges and points selected
        self.backregs_frame = tk.Frame(self.root,pady=10)
        self.backregs_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Simple spinbox for file selection in specselect frame
        self.specspin=tk.Spinbox(self.specselect_frame, command=self.on_specspin_change)
        # TODO does this need config before EDXdataset is loaded??
        self.specspin.pack(side=tk.TOP) # throw into specselect sub-frame
        # bools list that become true if any fitranges or backfitpts are altered
        self.fitflags = None
        # for readback of manually changed fitrange values 
        self.tkbegins= None # list with starting evs  of fitranges
        self.tkends= None # list with ending evs  of fitranges
        self.tkbackpts= None # list with background points in each fitrange
        
        # Replot button should link w/ plot in GUIplotter
        self.replot_button = tk.Button(
            self.misc_opts_frame,text="Replot",command=self.parent.plotter.plot(),
            padx=2, pady=6)
        self.replot_button.pack(side=tk.TOP,fill=tk.X,expand=1)

        self.refit_button = tk.Button(
            self.misc_opts_frame,text="Redo backfit", command=self.on_redo_backfit,
            padx=2, pady=6)
        self.refit_button.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.refit2_button = tk.Button(
            self.misc_opts_frame,text="Redo backfit all", command=self.on_redo_backfit_all,
            padx=2, pady=6)
        self.refit2_button.pack(side=tk.TOP,fill=tk.X,expand=1)

        self.readback_button = tk.Button(
            self.misc_opts_frame,text="Readback fitranges", command=self.read_backregs,
            padx=2, pady=6)
        self.readback_button.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.train_button = self._custom_button(
            self.misc_opts_frame,"Save backfit training", self.save_train)
        self.train_button.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.save_button = self._custom_button(
            self.misc_opts_frame,"Save EDXfile changes", self.on_save)
        self.save_button.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.quit_button = self._custom_button(
            self.misc_opts_frame,"Quit", self.on_quitapp)
        self.quit_button.pack(side=tk.TOP,fill=tk.X,expand=1)

    def associate_EDXdataset(self, EDXdataset):
        ''' associate loaded EDXdataset with GUIrefitter
        (passed as arg from GUIprojectloader) 
        called by GUIprojectloader '''
        self.EDXdataset= EDXdataset
        print('EDXdataset associated w/ GUIopts has ', len(EDXdataset.EDXlog),' files.')
        # Set specspin range (using zero-based indexing)
        self.specspin.config(from_=0, to=len(EDXdataset.EDXlog)-1)
        # clear any existing widgets in backreg frame
        for child in self.backregs_frame.winfo_children():
            child.destroy()
        # load first EDXfile into GUIplotter?
        self.load_EDXfile(0) # zero-based indexing so row zero
        # pass EDXfile laoded/created to GUIplotter
        self.parent.plotter.associate_EDXfile(self.EDXfile)
        # load background regions info from EDXfile into backregs_frame
        self.display_backregs() 
            
    def load_EDXfile(self, rowindex):
        ''' Load an EDXfile out of EDXdataset using dataframe row (.iloc) '''
        # Make instance of EDXfile class using parent (not EDXdataset itself) and rowindex 
        self.EDXfile=EDXfile(self.EDXdataset, rowindex)
        # Update displayed filename
        self.display_filename()
        # create fitflags of correct length
        self.fitflags=[False]*len(self.EDXfile.fitranges)
        ''' testing fit types problem
        for i, val in enumerate(self.EDXfile.fitorders):
            print('Region', i, 'fitorder is', val)
        '''

    def display_filename(self):
        ''' Displays csv name of currently-active emsa/csv file 
        called after every new load '''
        # clear filename display
        self.currfile_frame.grid_forget()
        tempstr='EDX Filename: '+self.EDXfile.filename
        tk.Label(self.currfile_frame, text=tempstr).pack()
        
    def on_specspin_change(self):
        ''' Load and plot chosen file, update backfit ranges, points, etc. 
        '''
        # clear old entries from any prior file 
        for child in self.backregs_frame.winfo_children():
            child.destroy()
        for child in self.currfile_frame.winfo_children():
            child.destroy()
        # EDXproject file must be loaded or no effect
        self.load_EDXfile(int(self.specspin.get()))
        # Update displayed fitregions, backfitpts
        self.display_backregs()

        # pass to GUIplotter
        self.parent.plotter.associate_EDXfile(self.EDXfile)

    def on_redo_backfit(self):
        ''' Update fitranges, backfitpts from display, then call refit method in EDXfile 
        linked to button '''
        # Note .. read back of fitranges, backpts done separately with button
        print('EDX background refitting initiated')
        self.EDXfile.process_refit(self.fitflags)
        print('EDX background refitting finished from GUIrefitter')
        # Pass updated EDXfile to plotter
        self.parent.plotter.associate_EDXfile(self.EDXfile)
        # reset fitflags to False
        self.fitflags=[False]*len(self.EDXfile.fitranges)

    def on_redo_backfit_all(self):
        ''' Refit of all regions ignoring fit flags
        Update fitranges, backfitpts from display, then call refit method in EDXfile 
        linked to button '''
        # Note .. read back of fitranges, backpts done separately with button
        print('EDX background refitting initiated')
        # set all to true to force complete refit (screwed up for some reason)
        self.EDXfile.process_refit([True]*len(self.fitflags))
        print('EDX background refitting finished from GUIrefitter')
        # Pass updated EDXfile to plotter
        self.parent.plotter.associate_EDXfile(self.EDXfile)
        # Reset fitflags to False
        self.fitflags=[False]*len(self.EDXfile.fitranges)
    
    def save_train(self):
        ''' call save training points method of currently-active EDXfile  
        training data about points added or removed from backfitpts 
        used to later improve fitting process
        '''
        # save any modified fitranges or backfitpts to backfitparamslog
        print('GUIrefitter save_train called')
        self.EDXfile.save_train()
        
    def on_save(self):
        ''' call save method of currently-active EDXfile  
        changeflags?? '''
        # save any modified fitranges or backfitpts to backfitparamslog
        print('GUIrefitter on_save called')
        self.EDXfile.save_backfits()
        # save EDXfile itself (with modified background column)
        self.EDXfile.save_csvfile()
    
    def on_quitapp(self):
        msg = "Quitting:\nUnsaved progress will be lost.\nDo you wish to Continue?"
        if tkmess.askokcancel("EDX refitter",msg):
            self.parent.root.destroy()
    
    def _custom_button(self,root,text,command,**kwargs):
        ''' use for lasso and point picker '''
        button = tk.Button(root, text=text,
            command=command,padx=2, pady=2,height=1, width=15,**kwargs)
        button.pack(side=tk.TOP,fill=None,expand=1)
        return button
    
    def remove_badpts(self, badinds):
        ''' Bad pt index #s in list returned from plotter after lasso-ing
        remove from EDXfile.Backfitpts, regenerate
        indices are effectively same as eV 
        Does badinds have index #s or actual energy vals in eV
        '''
        print('GUIrefitter remove_badpts called')
        # add points removed to existing list
        self.EDXfile.removedpts.extend(badinds)
        print(len(self.EDXfile.removedpts), ' points removed.')
        for fitnum, vals in enumerate(self.EDXfile.backfitpts):
            # See if badpts lie in this fit range
            common=[i for i in badinds if i in vals]
            if len(common)>0:
                # Troubleshoot remove pts error
                try:
                    self.fitflags[fitnum]=True # reset existing flag if change is made
                except:
                    print('Problem resetting flag', str(fitnum))
                # print('Removed ', len(common), 'faulty background fit points.')
                print(len(common),' values to remove for', fitnum)
            newvals=[i for i in vals if i not in badinds]
            # Check if rightmost or leftmost points in range have been removed
            [lowlim, hilim]=self.EDXfile.backptrange[fitnum]
            if lowlim in badinds or hilim in badinds:
                self.fix_badrange(fitnum, badinds)
                # alters background points and background ranges (but extends regions)
            else:
                self.EDXfile.backfitpts[fitnum]=newvals # list of lists
        self.display_backregs() # Update fitrange, backpts tkinter variables display
        # Update guiplot
        self.parent.plotter.associate_EDXfile(self.EDXfile)
    
    def fix_badrange(self, fitnum, badinds):
        ''' After bad point lasso removal, ensure that all regions still have 
        valid edge points 
        badinds are index numbers (should be same as vals stored in backfitparamslog)'''
        print('Fixing bad range after endpoint removed')
        # get all current backpoints from all regions
        allbackpts=self.EDXfile.get_allbackpts()
        # Make sure to remove bad points from list 
        allbackpts=[i for i in allbackpts if i not in badinds]
        # current boundaries for this fit region
        [lowlim, hilim]=self.EDXfile.backptrange[fitnum]
        # Current backpoints in this fit range
        currbackpts=self.EDXfile.backfitpts[fitnum]
        # Remove bad points
        currbackpts=[i for i in currbackpts if i not in badinds]
        if lowlim in badinds:
            # Find next smallest value
            print('Removing lower limit', lowlim)
            try:
                # Largest of negative differences is next lowest 
                newmin=lowlim+max([i-lowlim for i in allbackpts if i-lowlim<0])
            except:
                newmin=0
            # Reset backpts and associated range
            oldrange=self.EDXfile.backptrange[fitnum]
            self.EDXfile.backptrange[fitnum]=[newmin, oldrange[1]]
            # add 
            print('adding ', newmin, ' to region', fitnum,' backpts list')
            currbackpts.append(newmin)
        if hilim in badinds:
            print('Removing upper limit', lowlim)
            # Find next largest value 
            try:
                newmax=hilim+min([i-hilim for i in allbackpts if i-hilim>0])
            except:
                print('Problem removing largest background points value')
                #TODO fix for this problem 
            # Reset backpts and associated range
            oldrange=self.EDXfile.backptrange[fitnum]
            self.EDXfile.backptrange[fitnum]=[oldrange[0], newmax]
            # add 
            print('Adding ', newmax,' to region', fitnum,' backpts list')
            currbackpts.append(newmax)
        # Write changes back to this regions backfitpts
        self.EDXfile.backfitpts[fitnum]=currbackpts
            
    def add_pts(self, xval, yval):
        ''' Add closest single point (in energy) to background fit ranges'''
        # get index/eV of closest data point in energy col of EDX dataframe 
        print('GUIrefitter add_pts called')
        newval=self.EDXfile.EDXdf.Energy[(self.EDXfile.EDXdf.Energy-xval).abs().argsort()[:1]].index[0]
        # print('newval is', str(newval))
        # add points removed to existing list
        self.EDXfile.addedpts.append(newval)
        # Add newval to each backfitpts if within its fitrange
        for i, [fmin, fmax] in enumerate(self.EDXfile.fitranges):
            # evranges in format '0-100'
            if fmin  < newval < fmax:
                print('Newval', str(newval), 'added by addpts')
                self.fitflags[i]=True
                vals=self.EDXfile.backfitpts[i]
                vals.append(newval)
                vals.sort()
                self.EDXfile.backfitpts[i]=vals
            # Handle values greater than upper limit
            elif newval> fmax:
                self.fitflags[i]=True
                vals=self.EDXfile.backfitpts[i]
                vals.append(newval)
                vals.sort()
                self.EDXfile.backfitpts[i]=vals
                # alter total fit range
                self.EDXfile.fitranges[i]=[fmin, newval]
        self.display_backregs() # Update display
        self.parent.plotter.associate_EDXfile(self.EDXfile) # update guiplot

    def display_backregs(self):
        ''' Display fitranges, associated backpts for loaded EDXfile 
        '''
        # Clear any existing widgets in backreg frame
        for child in self.backregs_frame.winfo_children():
            child.destroy()
        # Write header row into backregs 
        rowframe=tk.Frame(self.backregs_frame)
        tk.Label(rowframe, text='Min').pack(side=tk.LEFT)
        tk.Label(rowframe, text='Max').pack(side=tk.LEFT)
        tk.Label(rowframe, text='Ptmin').pack(side=tk.LEFT)
        tk.Label(rowframe, text='Ptmax').pack(side=tk.LEFT)
        tk.Label(rowframe, text='#pts').pack(side=tk.LEFT)
        tk.Label(rowframe, text='Order').pack(side=tk.LEFT)
        rowframe.pack(fill=tk.X, expand=1)
        # Now display values associated w/ each 
        self.tkbegins=[] # list of tk string vars for fitrange beginnings 
        self.tkends=[] # fitrange ends
        # self.tkbackpts=[] # list of tk string vars for backpts
        self.tkptbegins=[]
        self.tkptends=[]
        self.tkfitorders=[] # fit types (linear (true) or parabola (false/default))
        # Unfortunately tk/mpl combo requires use of pack (not grid) 
        for i, [fmin, fmax] in enumerate(self.EDXfile.fitranges):
            # ev ranges are stored as "0-100" eV so needs parsing
            self.tkbegins.append(tk.IntVar())
            self.tkbegins[i].set(fmin)
            self.tkends.append(tk.IntVar())
            self.tkends[i].set(fmax)
            # self.tkbackpts.append(tk.StringVar())
            # bool var to keep track of linear (true) or parabola (false/default)
            self.tkfitorders.append(tk.IntVar())
            self.tkfitorders[i].set(self.EDXfile.fitorders[i])
            # backpoints are list of ints
            #print('backfitspts are of type ', type(self.EDXfile.backfitpts[i]))
            templist=self.EDXfile.backfitpts[i]
            # Something is rewriting EDXfile.backfitpts[i] to int
            # Set beginning of points included range
            self.tkptbegins.append(tk.StringVar())
            self.tkptbegins[i].set(str(min(templist)))
            # set end of points included range 
            self.tkptends.append(tk.StringVar())
            self.tkptends[i].set(str(max(templist)))
            
            # templist=[str(i) for i in templist]
            # tempstr=', '.join(templist)
            # self.tkbackpts[i].set(tempstr)
            # Add new row (via frame) .. .textvariable can be ints, right?
            rowframe=tk.Frame(self.backregs_frame)
            tk.Entry(rowframe, textvariable=self.tkbegins[i], width=5).pack(side=tk.LEFT)
            tk.Entry(rowframe, textvariable=self.tkends[i], width=5).pack(side=tk.LEFT)
            tk.Entry(rowframe, textvariable=self.tkptbegins[i], width=5).pack(side=tk.LEFT)
            tk.Entry(rowframe, textvariable=self.tkptends[i], width=5).pack(side=tk.LEFT)
            numpts=str(len(self.EDXfile.backfitpts[i]))
            tk.Label(rowframe, text=numpts, width=5).pack(side=tk.LEFT)
            # tk.Entry(rowframe, textvariable=self.tkbackpts[i]).pack(side=tk.LEFT)
            tk.Entry(rowframe, textvariable=self.tkfitorders[i], width=5).pack(side=tk.LEFT)
            rowframe.pack(fill=tk.X, expand=1)
            
    def read_backregs(self):
        ''' Readback manually altered fitranges, fitorders, and backpoints lists
        allows on-the-fly fit tweaking '''
        print('read_backregs started')
        # Set of bools keeping track of any altered params
        self.fitflags=[False]*len(self.EDXfile.fitranges)
        # Get old values for compare w/ readback 
        for i, [oldmin, oldmax] in enumerate(self.EDXfile.fitranges):
            # Also need to compare int lists of background points (readback vs current attributes)
            oldvals=self.EDXfile.backfitpts[i]
            oldbpmin=min(oldvals)
            oldbpmax=max(oldvals)
            # if backpoints ranges are changed, all backpts needed for every fitrange
            if int(self.tkbegins[i].get())!=oldmin or int(self.tkends[i].get())!=oldmax:
                # reset range and set of backpoints 
                self.EDXfile.fitranges[i]=[int(self.tkbegins[i].get()), 
                                      int(self.tkends[i].get())]
                self.fitflags[i]=True # keeping track of altered fitregions
                print('Background fit region', str(i),' changed')
            # Also need to check for backpoints changes w/o fitrange changes
            if int(self.tkptbegins[i].get())!=oldbpmin or int(self.tkptends[i].get())!=oldbpmax:
                # Need to reset backpoints in this fit region
                allback=self.EDXfile.get_allbackpts() # list of all backpts (ints)
                newvals=[i for i in allback if i >= oldbpmin and i <= oldbpmax]
                self.EDXfile.backfitpts[i]=newvals
                # reset range and set of backpoints 
                self.fitflags[i]=True # keeping track of altered fitregions
                print('Background fit points for region', str(i),' changed')
            # Check if fitorder has been changed
            if int(self.tkfitorders[i].get())!=self.EDXfile.fitorders[i]:
                self.EDXfile.fitorders[i]=int(self.tkfitorders[i].get())
                self.fitflags[i]=True
                print('Region',str(i),' changed to order', 
                      self.EDXfile.fitorders[i], ' polynomial')
    
    def runmenucommand(self, kwargs):
        ''' Method call from menu launched popup '''
        print('Running command', kwargs.get('command',''))
            
    def populate_specselector(self, spelist):
        ''' On project load, regenerate list of tk bools from spelist, update specselect frame view '''
        self.spec_tklist=[]
        for i, name in self.spelist:
            self.spec_tklist.append(tk.BooleanVar())
            self.spec_tklist[i].set(0) # Default unselected
            # Fill spectra selector frame w/ associated checkbuttons
            tk.Checkbutton(self.specselect_frame, text=name, variable=self.spec_tklist[i]).pack(side=tk.TOP)

