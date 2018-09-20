# -*- coding: utf-8 -*-
"""
Legacy functions for SEM-EDX import 
Created on Tue Apr  5 12:51:54 2016

@author: tkc
"""
# old method for entering sample names, project names, etc.
def fetch(entries): # used by get_project_info
   for entry in entries:
      field = entry[0]
      text  = entry[1].get()
      print('%s: "%s"' % (field, text)) 

def makeform(root, fields):
   entries = []
   for field in fields:
      row = Frame(root)
      lab = Label(row, width=15, text=field, anchor='w')
      ent = Entry(row)
      row.pack(side=TOP, fill=X, padx=5, pady=5)
      lab.pack(side=LEFT)
      ent.pack(side=RIGHT, expand=YES, fill=X)
      entries.append((field, ent))
   return entries
   
def get_project_info():
    fields='Project Name', 'Sample Category','Data Path' 
    root = Tk()
    ents=makeform(root, fields)
    root.bind('<Return>', (lambda event, e=ents: fetch(e)))   
    b1 = Button(root, text='Show',
        command=(lambda e=ents: fetch(e)))
    b1.pack(side=LEFT, padx=5, pady=5)
    b2 = Button(root, text='Quit', command=root.quit)
    b2.pack(side=LEFT, padx=5, pady=5)
    root.mainloop()
    return entries
entries=[]
projinfo=get_project_info()

def extractParams(emsafile): # pandas DataFrame
	'''Old method of extracting parameters from SEM-EDX files  ''' 
    emsaParams={'project':'','samplecategory':'','datapath':datapath,'spectra name':emsaFileName,'Date':'', 'Time':'', 'Title':'', 'Livetime':'', 'Realtime':'', 'Beamkv':'', 'Detect':'', 'Convert':'','Store':'','Timeconst':'','Deadfraction':''}    
    for i in range(0,100):
        tempstring=emsafile.ix[i,0]
        #if tempstring[:6]=="#TITLE": # just use spectra file name not embedded name
        #    emsaParam.update(title=emsafile.ix[i,1]
        if tempstring.startswith('#DATE'):
            date=emsafile.ix[i,1]
        if tempstring[:5]=="#TIME":
            time=emsafile.ix[i,1]
        if tempstring[:9]=="#LIVETIME":
            livetime=float(emsafile.ix[i,1])
        if tempstring[:9]=="#REALTIME":
            realtime=float(emsafile.ix[i,1])
        if tempstring[:9]=="##DETECTS":
            detects=int(emsafile.ix[i,1])
        if tempstring[:10]=="##CONVERTS":
            converts=int(emsafile.ix[i,1])
        else:
            converts=""        
        if tempstring[:8]=="##STORES":
            stores=int(emsafile.ix[i,1])
        else:
            stores=""        
        if tempstring[:12]=="##TIMECNSTNT":
            timeconst=int(emsafile.ix[i,1])
        else:
            timeconst=""        
        if tempstring[:7]=="#BEAMKV": # default 200kV for TEM but variable for SEM (5-20)
            beamkv=float(emsafile.ix[i,1])
        else:
            beamkv=""
        if tempstring[:9]=="#SPECTRUM":
            firstdatarow=i+1 # data begins on next row
    deadfraction=(realtime-livetime)/realtime
    
    # append params from emsa file i to emsaParamLog data frame
    newrow=[ None, date, time, title, livetime, realtime, beamkv, detects, converts, stores, timeconst, deadfraction ]
    emsaParamLog.loc[len(emsaParamLog)] = newrow # append to params dataframe
    if isinstance(lastrow,int):
        return lastrow  # should return lastrow to guide header removal from DataFrame
    else:
        return None

# Old compositional calculation before basis and errbasis were rolled into SEMimport		
def calccomp(df, Integquantlog, elemlist, SEMquantparams, sigthreshold=2):
    '''Calculate elemental composition of given files based on input element list 
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
    mycols=['Filenumber', 'Project', 'Basename', 'Point', 'Filename', 'FilePath', 'Sample', 'Comments','SEMbasis']
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
                    df=df.set_value(i, elem, temp.iloc[0]['Adjcnts']) # copy adjusted counts of this element
                    basis+=temp.iloc[0]['Subtractedcounts'] # add this element's value to AES basis
        # now handle the multipeak elements (get average value from both peaks)
        for key, value in multipeaklist.items(): # key is element (aka colname in df), value is list of peaks in Smdifpeakslog
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # this is number of lines that are average (i.e. 2 for Fe&Fe2)            
            avgval=0.0 # working value for averaged adjamplitude
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                temp=match[match['Element']==peak] # finds integquantlog entry for this peak (match already trimmed to filenum and area)
                if len(temp)==1:
                    if temp.iloc[0]['Significance']>sigthreshold:
                        avgval+=temp.iloc[0]['Subtractedcounts']
                    else:
                        numlines=numlines-1 # if peak is zeroed out and not added, this reduces # peaks in average
            if numlines>0: # avoid divbyzero if peak is too small
                avgval=avgval/numlines # this is now average basis for given element
            df=df.set_value(i, key, avgval) # copy adjusted amplitude of this element
            # add value from this element to SEMbasis
            basis+=avgval
        # end of multipeak elements loop
        df=df.set_value(i, 'SEMbasis', basis) # write total basis value to df
        # Now compute at.% for each listed element (incl errors)
        for j, elem in enumerate(elemlist):
            colname='%'+elem
            ratio=df.iloc[i][elem]/df.iloc[i]['SEMbasis'] # initialized to zero in cases where peak is below significance threshold
            df.set_value(i, colname, ratio)
            temp=match[match['Element']==elem] # again find peak entry and get finds entry for this peak
            # TODO maybe check threshold again (although element's value will be zero)
            if len(temp)==1: 
                thiserr=temp.iloc[0]['Erradjcnts']
                atpercerr=thiserr/df.iloc[i]['SEMbasis']
                errname='err%'+elem # error column
                df.set_value(i, errname, atpercerr) # Writes absolute error in at% 
        # Also calculate for elements w/ multiple peaks (if present)
        for key, value in multipeaklist.items(): 
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # this is number of lines that are average (i.e. 2 for Fe&Fe2)
            colname='%'+key
            ratio=df.iloc[i][key]/df.iloc[i]['SEMbasis']
            df.set_value(i, colname, ratio)
            # TODO need to propagate errors through Fe & Fe2
            errlist=[] # list of errors in % (usually max of two)
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                temp=match[match['Element']==peak] # finds entry for this peak
                if len(temp)==1:
                    if temp.iloc[0]['Adjcnts']>0: # skip negative values
                        err=temp.iloc[0]['Erradjcnts']/temp.iloc[0]['Adjcnts']
                        errlist.append(err) # add this to list 
            # combine errors in quadrature
            totalerr=0.0
            for j, err in enumerate(errlist):
                totalerr+=err**2
            totalerr=np.sqrt(totalerr) # percent error in at % 
            # now get  actual error
            thisval=df.iloc[i][key] # this is averaged value computed above (possibly zero if below thresholds )
            thiserr=thisval*totalerr # error (in Fe) as actual value based on average of multiple peaks
            atpercerr=thiserr/df.iloc[i]['SEMbasis']
            errname='err%'+ key  # error column
            df.set_value(i, errname, atpercerr) # Writes absolute error in at% 
        # end of loop calculation for each spectrum 
                
    # organize data based on mycols template
    df=organizecolumns(df,mycols)
    return df
	

# If log not created during session, just use params file instead of separate xls
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
        if match:
            basename=filename[0:match.start()] # pulls out base name
            specnum=int(filename[match.start()+1:match.end()-1]) # gets spectral number within parentheses
        else: # must have different naming convention (not auto-named)
            basename=filename.split('.')[0]
            specnum=1
        if '.psmsa' in filename: # emsa doesn't have point #
            ptnum=int(filename.split('_pt')[1].split('.')[0]) # gets pt number from 'name(10)_pt1.psmsa'
            Samplelogrow=Samplelogrow.set_value(0,'Point', ptnum)
        elif '.emsa' in filename: # no pt number for emsa ... set to 1
            ptnum=1
            Samplelogrow=Samplelogrow.set_value(0,'Point', ptnum)
        Samplelogrow=Samplelogrow.set_value(0,'Basename', basename)
        Samplelogrow=Samplelogrow.set_value(0,'Filenumber', specnum)

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
