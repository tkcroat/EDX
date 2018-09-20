# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 17:37:42 2017

@author: tkc
"""
import tkinter as tk

def EDXplot_tk(EDXlog, Elements, Backfitlog, AESquantparams, **kwargs):
    ''' tk interface for args/kwargs of AESplot1 function 
    all args/dataframes must be passed through to plot functions 
    '''
    # first print out existing info in various lines
    root = tk.Tk()
    filestr=tk.StringVar() # comma separated or range of filenumbers for plot
    filestr.set(kwargs.get('fileset','')) # get entry from prior run
    xrangestr=tk.StringVar()  # energy range in eV 
    xrangestr.set(kwargs.get('xrangestr',''))
    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    plotelemstr=tk.StringVar()
    mytext='Peaks to be labeled:'+', '.join(Elements)
    plotelemstr.set(mytext)
    newelems=tk.StringVar() # string for new element choices if made
    newelems.set('')
    choice=tk.StringVar()  # plot or abort

    a=tk.Label(root, text='Enter filenumbers for plotting').grid(row=0, column=0)
    b=tk.Entry(root, textvariable=filestr).grid(row=0, column=1)
    a=tk.Label(root, text='Enter xrange in eV (default 100-1900eV)').grid(row=2, column=0)
    b=tk.Entry(root, textvariable=xrangestr).grid(row=2, column=1)
        
    d=tk.Checkbutton(root, variable=backfitbool, text='Plot background fits?')
    d.grid(row=3, column=2)
    
    # option to reselect labeled elemental peaks 
    # TODO fix this nested tk interface ... chosen elements are not changing
    def changeelems(event):
        # pickelemsGUI imported from AESutils
        newelemlist=pickelemsGUI(AESquantparams) # get new elements/peaks list 
        newelems.set(', '.join(newelemlist))
        newtext='Peaks to be labeled: '+', '.join(newelemlist)
        plotelemstr.set(newtext)
        
    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    a=tk.Label(root, text=plotelemstr.get()).grid(row=5, column=0)
    d=tk.Button(root, text='Change labelled element peaks')
    d.bind('<Button-1>', changeelems)
    d.grid(row=6, column=0)

    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=6, column=2)

    d=tk.Button(root, text='Plot')
    d.bind('<Button-1>', plot)
    d.grid(row=6, column=1)

    root.mainloop()
        
    mychoice=choice.get()
    
    if mychoice=='plot':
        # Set up kwargs for plot 
        kwargs={}
        if newelems!='': # use newly assigned values
            kwargs.update({'plotelems':newelems.get()})            
        elif len(Elements)>0: # if not updated use passed list
            kwargs.update({'plotelems':Elements})
        if xrangestr.get()!='' and '-' in xrangestr.get(): # should be hyphenated range if altered
            xmin=int(xrangestr.get().split('-')[0])
            xmax=int(xrangestr.get().split('-')[1])
            kwargs.update({'xrange':(xmin, xmax)})       
            kwargs.update({'xrangestr':xrangestr.get()})
        if backfitbool.get():
            kwargs.update({'backfitdf':Backfitlog})
        if filestr.get()!='': # also put into kwargs to initialize next run
            kwargs.update({'fileset':filestr.get()})
            
        EDXplot1(filestr.get(), spelist, AESquantparams, **kwargs)

    return kwargs
