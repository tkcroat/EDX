
# EDX refit launcher 
root = tk.Tk()
root.wm_title("EDX refitter")
plotter = GUIMain(root)
root.mainloop()

# Working with OO data structures(used by GUIrefitter)
EDXset=EDXdataset('H:\\Research_data\\Thin_sections\\TEM_data\\LAP031117_3Nov14\\EMSA_files')
EDXset=EDXdataset('C:\\Temp\\EDX')
EDXset=EDXdataset(os.getcwd())
EDXset.EDXlog.iloc[0]['Filename']
len(EDXset.EDXlog)

# Load single dataframe file from above datasetEDXdf.energy
thisEDXfile=EDXfile(EDXset, 0)
thisEDXfile.filename
thisEDXfile.EDXdf
thisEDXfile.corrcnts
thisEDXfile.quantelems

launch_refitter(os.getcwd())

# test version with pop-up menus to get args 
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
        # Menubars
        self.menubar=tk.Menu(self.root)
        filemenu = tk.Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="Open", command=self.hello)
        filemenu.add_command(label="Save", command=self.hello)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.refitter.on_quitapp)
        self.menubar.add_cascade(label="File", menu=filemenu)
        self.root.config(menu=self.menubar)

        specmenu = tk.Menu(self.menubar, tearoff=0)
        specmenu.add_command(label="Popup", command=lambda: 
            self.args_menu({'command':'median_filter', 'filter size':1}))
        self.menubar.add_cascade(label="Spectral Commands", menu=specmenu)
        self.root.config(menu=self.menubar)

        # TODO pop up on right click not really what is needed 
        # new temporary window for entering kwargs (multi-level menu)

    def args_popup_menu(self, kwargs):
        ''' Menu launched toplevel window for args/kwargs entry and calls
        to 
        kwargs: command -  name of method in 
            param name & associated value (e.g.  kwargs={'filter size':1})
        '''
        def abort():
            ''' close popup '''
            t.destroy()
        def run():
            ''' close popup '''
            self.refitter.runmenucommand(kwargs)
            t.destroy()
        t = tk.Toplevel(self.root) # open new toplevel window
        tkvars=[] # display and alter params passed in kwargs
        for i, (key, val) in enumerate(kwargs.items()):
            tkvars.append(tk.StringVar())
            tkvars[i].set(val)
            if key!='command': # 
                prow=tk.Frame(t)
                tk.Label(prow, text=key).pack(side=tk.LEFT)
                tk.Entry(prow, textvariable=tkvars[i]).pack(side=tk.LEFT)
                prow.pack(side="top", fill="both", expand=True, padx=100, pady=100)
        prow=tk.Frame(t)
        tk.Button(prow, text='Abort', command=abort).pack(side=tk.LEFT)
        mystr='Run '+kwargs.get('command','')
        tk.Button(prow, text=mystr, command=run).pack(side=tk.LEFT)
        prow.pack(side="top", fill="both", expand=True, padx=100, pady=100)
        

            
            
