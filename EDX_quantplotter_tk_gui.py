# -*- coding: utf-8 -*-
"""
Interactive EDX background refitter 
Created on Wed Oct 11 00:44:29 2017

@author: tkc
"""
import os
import pandas as pd
import glob
import numpy as np
import tkinter as tk
import tkinter.messagebox as tkmess
from tkinter import filedialog
import matplotlib as mpl # using path, figure, rcParams
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import sys
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

#-----------Misc.-------------#

def launch_plotter(currdir):
    ''' Launcher function for tk quantplotter EDX GUI '''
    root = tk.Tk()
    root.wm_title("EDX refitter")
    GUIMain(root, currdir)
    root.mainloop()
    return

class GUIMain():
    ''' Main container for plotter, options (at right), and fileloader (bottom) 
    pass current working directory as default directory'''
    def __init__(self,root, currdir):
        self.root = root
        self.root.wm_title("EDX quant plotter")
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
        self.opts  = GUIOpts(self.refit_frame,self)
        self.loader  = GUIprojectloader(self.loader_frame,self, currdir)

class GUIPlotter():
    def __init__(self,root, parent):
        self.root = root
        self.parent = parent
        self.figure = mpl.figure.Figure(figsize=PLOT_SIZE, dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
        self.canvas = FigureCanvasTkAgg(self.figure,self.root)
        # just use standard toolbar
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,self.root)
        self.toolbar.update()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.EDXfile = None
        self.spectype= None # SEM or TEM spectra (from plotted EDXfile)
        self.canvas.show()

    def associate_EDXfile(self, EDXfile):
        ''' Associate EDX file created/loaded in GUIoptions with GUIplotter 
        called from GUIoptions '''
        # Is this a separate instance from that in guiopts?
        print('New EDXfile associated with plot window')
        self.EDXfile = EDXfile
        self.spectype=EDXfile.spectype
        # SEM-EDX or TEM-EDX ... different xautoscale
        self.plot() # 
    
    def auto_rescale(self):
        ''' Reset scaling depending on spectype '''
        if self.spectype=='SEM':
            self.ax.set_xlim([0.1,10]) # auto setting of scale
        elif self.spectype=='TEM':
            self.ax.set_xlim([0.1,20]) # auto setting of scale
            
    def plot(self,**kwargs):
        if self.EDXfile is None:return
        self.ax.cla() # clear axes
        self.ax.plot(self.EDXfile.energy, self.EDXfile.EDXdf['Counts'], color='k', picker=True,**kwargs)
        # Add existing background fit in red 
        self.ax.plot(self.EDXfile.energy, self.EDXfile.EDXdf['Backfit'], color='r', picker=True,**kwargs)
        self.auto_rescale()
        self.canvas.show()

    def plot_elems(self, elemparams):
        if self.EDXfile is None:return
        self.ax.cla() # clear
        self.plot() # add plot and backfit again
        # Add vertical lines at known element energies (empty elemparams passed to toggle off)
        for i, [elem, energyval] in enumerate(elemparams):
            self.ax.axvline(x=energyval, color='b', linestyle='dashed', label=elem)
        if len(elemparams)!=0:
            xmax=int(max(map(lambda x:x[1],elemparams))*1.2) # max of 2nd element in list of lists 
            self.ax.set_xlim([0.1,xmax])
        # Add text labels at approx position
        self.canvas.show()

    def label_quant(self, elems, vals):
        ''' Add quant text label with active elems and at. % values '''
        if self.EDXfile is None:return
        # Add vertical lines at known element energies
        fullstr=''
        for i, (elem,val) in enumerate(zip(elems, vals)):
            tempstr=r'$%s_{%.0f}$' %(elem, float(val))
            fullstr+=tempstr            
        # transform=ax.transAxes uses coords from 0 to 1 (instead of true x and y vals)
        self.ax.text(0.05,0.95, fullstr, fontsize=30, verticalalignment='top', transform=self.ax.transAxes)
        self.canvas.show()
        
class GUIprojectloader():
    ''' Picks directory and loads main Auger param files
    needs current path (which should be set to working data directory) '''
    def __init__(self,root,parent, currdir):
        self.root = root
        self.parent = parent # GUImain is parent
        self.top_frame = tk.Frame(self.root)
        self.top_frame.pack(side=tk.TOP,anchor=tk.W)
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM,fill=tk.BOTH,expand=1)
        tk.Label(self.top_frame,text="Directory:",padx=8,pady=2,height=1).pack(side=tk.LEFT,anchor=tk.W)
        self.directory_entry = tk.Entry(self.top_frame,width=90,bg="lightblue",
                                    fg="black",highlightcolor="lightblue",insertbackground="black",
                                    highlightthickness=2)
        self.directory_entry.pack(side=tk.LEFT,fill=tk.BOTH,expand=1,anchor=tk.W)
        self.directory_entry.insert(0,currdir)
        tk.Button(self.top_frame,text="Browse",command=self.launch_dir_finder
                  ).pack(side=tk.LEFT,fill=tk.BOTH,expand=1,anchor=tk.W)
        self.load_button = tk.Button(self.bottom_frame,text="Load EDX project folder",width=60,
                                       command=self.load_EDXdataset)
        self.load_button.pack(side=tk.BOTTOM,expand=1,anchor=tk.CENTER)
        
    def load_EDXdataset(self):
        ''' Load standard AES files (paramlogwith data returned to a DataManager '''
        directory = self.directory_entry.get()
        EDXdata = EDXdataset(directory)
        # pass to GUIoptions and set spectral selector spinbox values
        self.parent.opts.associate_EDXdataset(EDXdata)        
        # TODO associate EDXdataset with GUIplotter (or just selected EDXfile)
    
    def launch_dir_finder(self):
        directory = filedialog.askdirectory()
        self.directory_entry.delete(0,tk.END)
        self.directory_entry.insert(0,directory)

class GUIOpts():
    ''' Parent is GUImain, manages EDXfile displayed in GUIplotter
    handles addition/removal of points for background (re)fitting'''
    def __init__(self,root,parent):
        self.root = root
        self.parent = parent
        self.EDXdataset = None  # created in GUIprojectloader but associated here
        # instance of EDXfile local to the refitter
        self.EDXfile = None
        # spinbox 
        self.specselect_frame = tk.Frame(self.root,pady=10)
        self.specselect_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # for display of currently plotted file
        self.currfile_frame = tk.Frame(self.root,pady=10)
        self.currfile_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Element selector checkboxes
        self.elems_frame = tk.Frame(self.root,pady=10)
        self.elems_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Frame for display of counts/quant results
        self.quant_frame = tk.Frame(self.root,pady=10)
        self.quant_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        self.misc_frame = tk.Frame(self.root,pady=10)
        self.misc_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Simple spinbox for file selection in specselect frame
        self.specspin=tk.Spinbox(self.specselect_frame, command=self.on_specspin_change)
        # TODO does this need config before EDXdataset is loaded??
        self.specspin.pack(side=tk.TOP) # throw into specselect sub-frame
        self.tkelems=[] # bools list for elem display or quant
        self.activequant=[] # for at % results 
        self.showelems=False # toggle for showing elemental lines on plot

        # Element presets (top of misc frame)
        rowframe=tk.Frame(self.misc_frame)
        tk.Button(rowframe, text='Clear all', command=self.clearall).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Button(rowframe, text='Select all', command=self.selectall).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Button(rowframe, text='Si Fe Mg S Ca O', command=self.elempreset).pack(side=tk.LEFT,fill=tk.X,expand=1)
        rowframe.pack(fill=tk.X, expand=1)
        
        # permanent buttons in misc_frame
        self.label_button = tk.Button(
            self.misc_frame,text="Label elements",command=self.label_elems)
        self.label_button.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.quant_button = tk.Button(
            self.misc_frame,text="Update quant", command=self.do_quant)
        self.quant_button.pack(side=tk.TOP,fill=tk.X,expand=1)

        self.quant_button2 = tk.Button(
            self.misc_frame,text="Add quant label", command=self.label_quant)
        self.quant_button2.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.quit_button = tk.Button(
            self.misc_frame, text="Quit", command=self.on_quitapp)
        self.quit_button.pack(side=tk.TOP,fill=tk.X,expand=1)
    
    def elempreset(self):
        ''' Clear selected elements '''
        if not self.EDXfile: return
        presets=['S','Fe','Mg','Si','Ca','O']
        for i, elem in enumerate(self.EDXfile.quantelems):
            if elem in presets:
                self.tkelems[i].set(1)
            else:
                self.tkelems[i].set(0)
    
    def selectall(self):
        ''' Clear selected elements '''
        for i, tkbool in enumerate(self.tkelems):
            self.tkelems[i].set(1)

    def clearall(self):
        ''' Clear selected elements '''
        for i, tkbool in enumerate(self.tkelems):
            self.tkelems[i].set(0)
            
    def display_elems(self):
        ''' Display available quant elements in elems_frame (self.EDXfile.quantelems only 
        contains elements that already have quant results (not other random peaks) '''
        for child in self.elems_frame.winfo_children():
            child.destroy()
        # Need to call frame header Label and 
        # Write header row into backregs
        tk.Label(self.elems_frame, text='Available Elements').pack(side=tk.TOP,fill=tk.X,expand=1)

        # tkelems bool variables for active/inactive for each element
        self.tkelems=[]
        for i, elem in enumerate(self.EDXfile.quantelems):
            self.tkelems.append(tk.BooleanVar())
            self.tkelems[i].set(True)
        # Unfortunately tk/mpl combo requires use of pack (not grid) 
        for i in range(0, len(self.EDXfile.quantelems), 3):
            # associated checkbutton for each quantelem
            elemlistframe=tk.Frame(self.elems_frame)
            tk.Checkbutton(elemlistframe, variable=self.tkelems[i]).pack(side=tk.LEFT)
            tk.Label(elemlistframe, text=self.EDXfile.quantelems[i]).pack(side=tk.LEFT)
            try:
                tk.Checkbutton(elemlistframe, variable=self.tkelems[i+1]).pack(side=tk.LEFT)
                tk.Label(elemlistframe, text=self.EDXfile.quantelems[i+1]).pack(side=tk.LEFT)
            except: # out of list range problem
                pass
            try:
                tk.Checkbutton(elemlistframe, variable=self.tkelems[i+2]).pack(side=tk.LEFT)
                tk.Label(elemlistframe, text=self.EDXfile.quantelems[i+2]).pack(side=tk.LEFT)
            except:
                pass
            elemlistframe.pack(fill=tk.X, expand=1)
            
    def associate_EDXdataset(self, EDXdataset):
        ''' associate loaded EDXdataset with GUIoptions 
        (passed as arg from GUIprojectloader) 
        called by GUIprojectloader '''
        self.EDXdataset= EDXdataset
        print('EDXdataset associated w/ GUIopts has ', len(EDXdataset.EDXlog),' files.')
        # Set specspin range (using zero-based indexing)
        self.specspin.config(from_=0, to=len(EDXdataset.EDXlog)-1)
        # clear any existing widgets in backreg frame
        for child in self.elems_frame.winfo_children():
            child.destroy()
        # load first EDXfile into GUIplotter?
        self.load_EDXfile(0) # zero-based indexing so row zero
        # pass EDXfile laoded/created to GUIplotter
        self.parent.plotter.associate_EDXfile(self.EDXfile)
        # loads quant elements into elems frame
        self.display_elems()
            
    def load_EDXfile(self, rowindex):
        ''' Load an EDXfile out of EDXdataset using dataframe row (.iloc) '''
        # Make instance of EDXfile class using parent (not EDXdataset itself) and rowindex 
        self.EDXfile=EDXfile(self.EDXdataset, rowindex)
        # Update displayed filename
        self.display_filename()

    def label_quant(self):
        ''' Add a text label with current quant results to plotter 
        launched via button  '''
        elems=[i[0] for i in self.activequant]
        vals=[f[1] for f in self.activequant]
        vals=[int(i) if i>1 else "%0.1f" % i for i in vals]
        self.parent.plotter.label_quant(elems, vals)
        
    def display_filename(self):
        ''' Displays csv name and sample name of currently-active emsa/csv file 
        called after every new load '''
        # clear filename display
        for child in self.currfile_frame.winfo_children():
            child.destroy()
        tempstr='EDX Filename: '+self.EDXfile.filename
        tk.Label(self.currfile_frame, text=tempstr).pack()
        tempstr='Sample name: '+self.EDXfile.sample
        tk.Label(self.currfile_frame, text=tempstr).pack()
        
    def on_specspin_change(self):
        ''' Load and plot chosen file, update backfit ranges, points, etc. 
        '''
        # clear old entries from any prior file 
        for child in self.currfile_frame.winfo_children():
            child.destroy()
        for child in self.elems_frame.winfo_children():
            child.destroy()
        for child in self.quant_frame.winfo_children():
            child.destroy()
        # EDXproject file must be loaded or no effect
        self.load_EDXfile(int(self.specspin.get()))
        # Update displayed fitregions, backfitpts
        self.display_filename()
        persistlist=[] 
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                persistlist.append(True)
            else:
                persistlist.append(False)
        self.display_elems() # recreates tkelems bools list
        # reselect the same elements
        if len(persistlist)==len(self.tkelems):
            for i, thisbool in enumerate(persistlist):
                self.tkelems[i].set(thisbool)

        # pass to GUIplotter
        self.parent.plotter.associate_EDXfile(self.EDXfile)
        # Also rerun label_elements (lines stay on if set to on)
        if self.showelems: # needs double toggle to stay the same
            self.showelems=False
        else:
            self.showelems=True
        self.label_elems()
    
    def on_quitapp(self):
        msg = "Quitting:\nUnsaved progress will be lost.\nDo you wish to Continue?"
        if tkmess.askokcancel("EDX refitter",msg):
            self.parent.root.destroy()
    
    def label_elems(self):
        ''' Get currently selected elements from elems frame and associated energy
        pass to plotter.label_elems
        toggle style button
        '''
        if self.showelems:
            self.showelems=False
        else:
            self.showelems=True
        elemparams=[] # list of [elem, energyval]
        # active (checked) tk elements will automatically be updated on check, right?
        if self.showelems:
            for i, tkbool in enumerate(self.tkelems):
                if tkbool.get():
                    match=self.EDXdataset.EDXquantparams[self.EDXdataset.EDXquantparams['element']==self.EDXfile.quantelems[i]]
                    ev=match.iloc[0]['energy']
                    elemparams.append([self.EDXfile.quantelems[i],ev])
        # now pass info to plotter (can be empty list)
        self.parent.plotter.plot_elems(elemparams)
    
    def do_quant(self):
        ''' Generate at % and error for subset of selected elements
        then display quant ... linked to quant button '''
        # clear current values
        self.activequant=[] # list w/ [elem symb, atperc, erratperc, adjcnts, corrcnts, errcorrcnts, ]
        # note erratperc (3rd in list) is actual error not % err
        temperrperc=[] # error % calc by Errcorrcnts/corrcnts (integlog has actual error in corrcnts)
        corrsum=0.0
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                self.activequant.append([self.EDXfile.quantelems[i], 0.0, 0.0, self.EDXfile.adjcnts[i],self.EDXfile.corrcnts[i],self.EDXfile.errcorrcnts[i]])
                if self.EDXfile.corrcnts[i]>0:
                   temperrperc.append(self.EDXfile.errcorrcnts[i]/self.EDXfile.corrcnts[i])
                   corrsum+=self.EDXfile.corrcnts[i]
                else: # Exclude negative values but list entry needed
                    temperrperc.append(0.0) # percent error (intermediate calculation)
        for i, errperc in enumerate(temperrperc):
            self.activequant[i][1]=100*self.activequant[i][4]/corrsum
            self.activequant[i][2]=(100*self.activequant[i][4]/corrsum)*temperrperc[i]
        self.display_quant()    
        
    def display_quant(self):
        ''' Display fitranges, associated backpts for loaded EDXfile 
        '''
        # Clear any existing widgets in backreg frame
        for child in self.quant_frame.winfo_children():
            child.destroy()
        # sort active quant elements by at percent
        self.activequant.sort(key=lambda x: float(x[1]), reverse=True)
        # Write header row into backregs 
        rowframe=tk.Frame(self.quant_frame)
        tk.Label(rowframe, text='Element').pack(side=tk.LEFT)
        tk.Label(rowframe, text='At%').pack(side=tk.LEFT)
        tk.Label(rowframe, text='err at %').pack(side=tk.LEFT)
        tk.Label(rowframe, text='Adjcnts').pack(side=tk.LEFT)
        #tk.Label(rowframe, text='Corrcnts').pack(side=tk.LEFT)
        #tk.Label(rowframe, text='Err corrcnts').pack(side=tk.LEFT)
        rowframe.pack(fill=tk.X, expand=1)
        # For values display len(quantelems)!=len(activeelems)
        for i, [elem, atper, erratper, adjcnts, corrcnts, errcorrcnts] in enumerate(self.activequant):
            rowframe=tk.Frame(self.quant_frame)
            tempstr=elem+'   '
            tk.Label(rowframe, text=tempstr).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.1f" % atper).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.1f" % erratper).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.0f" % adjcnts).pack(side=tk.LEFT)
            #tk.Label(rowframe, text="%.1f" % corrcnts).pack(side=tk.LEFT)
            #tk.Label(rowframe, text="%.1f" % errcorrcnts).pack(side=tk.LEFT)
            rowframe.pack(fill=tk.X, expand=1)