# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:05:50 2017

@author: tkc
"""
import numpy as np
import glob
import os
import pandas as pd
from scipy import optimize

''' TESTING
EDXdf=EDXdataset('C:\\Temp\\SiC\\MD2d_11Jun10')
MyEDX=EDXfile(EDXdf,0)
''' 

class EDXfile():
    ''' Single spectral file loaded from row of EDXdataset
    same for EDX refitter and quantplotter '''
    def __init__(self, EDXdataset, rowindex, **kwargs):        
        ''' Single instance created from EDX dataframe row (pandas series) '''
        # retrieve desired row from EDXlog file
        self.EDXdataset=EDXdataset # needed to access within methods
        self.row=EDXdataset.EDXlog.iloc[rowindex] # this file's row in EDXlog
        self.filename=self.row.Filename # emsa filename
        self.sample=str(self.row.Sample) # sample name(could be nan)   
        self.beamkv=self.row.Beamkv
        self.livetime=self.row.Livetime
        self.deadfraction=self.row.Deadfraction
        self.EDXdf = None # entire dataframe (w/ energy, counts, backfit, subdata columns)
        self.energy=None # occasionally used but never altered
        self.fitranges = None # index numbers (same as eV) used for fitting
        self.fitorders= None  # lits of fittypes for each region
        self.backfitpts = None # subset of pts used for fits 
        # this is list of lists w/ points for each backfit region
        self.backptrange = None # list of ranges of backpoints used (min/max)
        self.backfitparams = None # A, B, C coeffs resulting from fitting
        # altered backfitpoints for training (vals set in GUIrefitter)
        self.origbackfitpts = None # single flattened list from backfitpts
        self.removedpts = [] # emply lists are best
        self.addedpts = []
        
        self.open_csvfile() # opens data file
        self.get_backfitregs() # Opens associated backfitlog from parent
        # attribs pulled from integration log (integlog of EDXdataset)
        self.quantelems = None
        self.adjcnts = None
        self.corrcnts = None
        self.errcorrcnts = None
        self.get_quantelems() # find elements w/ existing quant info
        self.Elemdata = None
        self.spectype=EDXdataset.spectype # SEM or TEM
            
    def get_quantelems(self):
        ''' Finds element quant already performed from integlog (within EDXdataset)
        '''
        match=self.EDXdataset.Integlog.loc[ (self.EDXdataset.Integlog['Filename']==self.filename)]
        # should contain row for each element included in quant
        self.quantelems=[]
        self.adjcnts=[]
        self.corrcnts=[]
        self.errcorrcnts=[]
        for index, row in match.iterrows():
            self.quantelems.append(row.Element)
            self.corrcnts.append(row.Correctedcounts)
            self.errcorrcnts.append(row.Errcorrcnts)
            if str(row.Adjcounts)=='nan':
                self.adjcnts.append(row.Subtractedcounts)
            else:
                self.adjcnts.append(row.Adjcounts)
        
    def get_backfitregs(self):
        ''' From backfitlog get list of lists of fitting ranges and backfitpts
        only called on init '''
        # df with ev ranges , backpoints, fittype, and parabola or lin params (ie. A,B,C)
        match=self.EDXdataset.Backfitlog.loc[ (self.EDXdataset.Backfitlog['Filename']==self.filename)]
        # list of int lists containing backfitpts fo each region
        self.fitranges=[]
        self.backfitpts=[]
        self.backptrange=[]
        self.fitorders=[]
        # TODO Error on reopen related to string or list conversion
        for index, row in match.iterrows():
            # Backfitpts can be comma-sep string or python list (first open or reopen)
            if isinstance(row.Backfitpts, list):
                self.backfitpts.append(row.Backfitpts)
            else:
                tempstr=row.Backfitpts.replace('[','').replace(']','')
                mylist=tempstr.split(',')
                mylist=[int(i) for i in mylist]
                self.backfitpts.append(mylist)
                self.backptrange.append([min(mylist), max(mylist)])
            # Convert fittype to fitorder (i.e. linear=1, parabola =2, cubic=3)
            if row.Fittype.strip()=='linear':
                self.fitorders.append(1)
            elif row.Fittype.strip()=='parabola':
                self.fitorders.append(2)
            if row.Fittype.strip()=='cubic':
                self.fitorders.append(3)
            # Convert fitranges from string 0-100 to min, max
            self.fitranges.append([int(i) for i in row.Fitrange.split('-')])            
            
        self.backfitparams=[] # Current fit params for each subregions
        for index, row in match.iterrows():
            self.backfitparams.append([row.A, row.B, row.C, row.D])
        # Go ahead and set original backfit points list after loading
        self.origbackfitpts=self.get_allbackpts()

    def get_allbackpts(self):
        ''' get all unique active background points from all fit regions
        often called by refitter
        '''
        allbackpts=[]
        for i, bplist in enumerate(self.backfitpts):
            allbackpts.extend(bplist)
        allbackpts=set(allbackpts)
        allbackpts=list(allbackpts)
        return allbackpts
            
    def process_refit(self, flags):
        ''' Parabolic fit over same range but possibly with points removed 
        list of lists with evrange and associated background points
        flags list passed by GUIoptions/ refitter 
        changes to EDXfile fitrange or backfitpts made prior to process_refit call
        changing fitrange w/o changing backfitpts won't change polynomial refit but 
        will change creation of piece-wise background 
        fitranges alterable using readback of display'''
        # Pass back any possibly modified params from GUIopts
        print('EDXfile process refit started')
        for i, mybool in enumerate(flags):
            # print('Fitflags are ', str(i), mybool)
            if mybool:
                # Decide on parabolic or linear
                self.refit_poly(i)
                print('Start refitting of region', str(i),)
        if any(flags): # If any sub-regions were refit the regenerate full background 
            self.recreate_backfit()
        # need to update plot and reset flags
        print('Finished with EDXfile process refit')

    def change_backrange(self, fitnum, newrange):
        # TODO move these alterations from GUIrefitter to here
        ''' Make change in backptrange (min, max) and readjust backpoints 
        contained in the range
        newrange is [min, max]'''
        # Pass back any possibly modified params from GUIopts
        if self.backptrange[fitnum]!=newrange:
            self.backptrange[fitnum]=newrange
            allbackpts=self.get_allbackpts()
            newpts=[i ]
        print('Change backpoint range and points')
        for i, mybool in enumerate(flags):
            self.refit_poly(i)
            print('Start refitting of region', str(i),)
        self.recreate_backfit()
        # need to update plot and reset flags
        print('Finished with EDXfile process refit')
        
    def refit_all(self, flags):
        ''' Parabolic fit over same range but possibly with points removed 
        list of lists with evrange and associated background points
        unlike process_refit this ignores flags and does complete refit
        changes to EDXfile fitrange or backfitpts made prior to process_refit call
        changing fitrange w/o changing backfitpts won't change polynomial refit but 
        will change creation of piece-wise background 
        fitranges alterable using readback of display'''
        # Pass back any possibly modified params from GUIopts
        print('EDXfile process refit started')
        for i, mybool in enumerate(flags):
            self.refit_poly(i)
            print('Start refitting of region', str(i),)
        self.recreate_backfit()
        # need to update plot and reset flags
        print('Finished with EDXfile process refit')
        
    def refit_poly(self, i):
        ''' Call polynomial refit of parabola for given region (i value indicates which one in list)
        needed from self:  counts, fitranges, backfitpoints, i value
        fitting over points range (not fitrange boundary which may be a wider region (used for cross-
        fading across boundaries  '''
        print('starting refit_poly for region', str(i))
        thisreg=self.EDXdf['Counts'][self.EDXdf['Counts'].index.isin(self.backfitpts[i])]
        thisrange=self.energy[self.energy.index.isin(self.backfitpts[i])]
        # print('length of self.backfitpts[i] is', len(self.backfitpts[i]))
        print('Polynomial fit of order', self.fitorders[i])
        newparams=np.polyfit(thisrange, thisreg, self.fitorders[i])
        newparams=np.ndarray.tolist(newparams)
        tempstr=[str(i) for i in newparams]
        print('New fit values are:',','.join(tempstr))
        while len(newparams)<4:
            newparams.append(np.nan) # pad with nan values for linear, parabolic
        self.backfitparams[i]=newparams
            
    def recreate_backfit(self):
        ''' Regenerate full range background fit from existing piece-wise fits
        energy, counts,backfit (all pd series); '''
        print('Proceeding with recreation of backfit')
        for i, [lower, upper] in enumerate(self.fitranges):
            # Get index numbers of this fit's boundaries (indexes not eV, right)
            if self.backfitparams[i][1]!='n/a': # Check for prior failed fit
                for j in range(lower, upper):
                    xval=self.energy[j]
                    if self.fitorders[i]==2: # parabola
                        [A,B,C]=[self.backfitparams[i][0], self.backfitparams[i][1],self.backfitparams[i][2]]
                        self.EDXdf['Backfit']=self.EDXdf['Backfit'].set_value(j, A * xval**2 + B * xval + C)
                    elif self.fitorders[i]==1: # linear
                        [A,B]=[self.backfitparams[i][0], self.backfitparams[i][1]]
                        self.EDXdf['Backfit']=self.EDXdf['Backfit'].set_value(j, A * xval+ B)
                    elif self.fitorders[i]==3: # cubic
                        [A,B,C, D]=[self.backfitparams[i][0], self.backfitparams[i][1],self.backfitparams[i][2],self.backfitparams[i][3]]
                        self.EDXdf['Backfit']=self.EDXdf['Backfit'].set_value(j, A * xval**3 + B * xval**2 + C* xval + D)
                        # Now find and cross-fade overlapping regions
        overlaps=[] # entry exists for each region boundary (even if not overlapping)
        for i in range(1, len(self.fitranges)):
            [start1, end1]=self.fitranges[i-1]        
            [start2, end2]=self.fitranges[i]
            # Lower of the overlapped regions and overlap range 
            overlaps.append([i-1, start2, end1])
        # Cross-fade background regions (assuming parabolic fits)
        for i, [regnum, start, end] in enumerate(overlaps):
            # overlap num not necessarily same as backfitregion 
            # Regnum is number of lower of two overlapping regions 
            # start and end of overlap region 
            print('Handling overlapping region: ', start, '-', end, 'between regs', regnum, 
                  'and ', regnum+1)
            # Check to ensure both are parabolic fits
            if self.fitorders[regnum]==1: # linear
                C0=self.backfitparams[regnum][0]
                D0=self.backfitparams[regnum][1]
                B0=0
                A0=0
            elif self.fitorders[regnum]==2: # parabola
                A0=0
                B0=self.backfitparams[regnum][0]
                C0=self.backfitparams[regnum][1]
                D0=self.backfitparams[regnum][2]
            elif self.fitorders[regnum]==3: # cubic
                A0=self.backfitparams[regnum][0]
                B0=self.backfitparams[regnum][1]
                C0=self.backfitparams[regnum][2]
                D0=self.backfitparams[regnum][3]
            else:
                print('Unknown fit type')
            if self.fitorders[regnum+1]==1:
                A1=0
                B1=0
                C1=self.backfitparams[regnum+1][0]
                D1=self.backfitparams[regnum+1][1]
            elif self.fitorders[regnum+1]==2:
                A1=0
                B1=self.backfitparams[regnum+1][0]
                C1=self.backfitparams[regnum+1][1]
                D1=self.backfitparams[regnum+1][2]
            elif self.fitorders[regnum+1]==3:
                A1=self.backfitparams[regnum+1][0]
                B1=self.backfitparams[regnum+1][1]
                C1=self.backfitparams[regnum+1][2]
                D1=self.backfitparams[regnum+1][3]
            else:
                print('Unknown fit type')
            thisrange=abs(end-start) # total for crossfading
            for j in range(start, end):
                xval=self.energy[j]
                yval=(1-(j-start)/thisrange)*(D0+C0*xval+B0*xval**2+A0*xval**3)+((j-start)/thisrange)*(D1+C1*xval+B1*xval**2+A1*xval**3)
                self.EDXdf['Backfit']=self.EDXdf['Backfit'].set_value(j, yval)
            # Immediate update of EDXdf and subdata (also done i save_csvfile)
            self.EDXdf['Subdata']=self.EDXdf['Counts']-self.EDXdf['Backfit']

    def open_csvfile(self):
        ''' read single edx csv file '''
        # needs to handle emsa or psmsa
        self.EDXdf=pd.read_csv(self.filename.replace('.emsa','.csv').replace('.psmsa','.csv'))
        self.energy=self.EDXdf['Energy']
        print('EDXfile ', self.filename,' loaded.')
    
    def save_csvfile(self):
        ''' Direct save of open csv file (after any modifications) 
        saves to current working directory (counts and energy unchanged)'''
        self.EDXdf['Subdata']=self.EDXdf['Counts']-self.EDXdf['Backfit'] # recompute subtracted data
        if '.emsa' in self.filename:
            self.EDXdf.to_csv(self.filename.replace('.emsa','.csv'), index=False)
        elif '.psmsa' in self.filename:
            self.EDXdf.to_csv(self.filename.replace('.psmsa','.csv'), index=False)
        print(self.filename.split('.')[0],'.csv saved.')
    
    def save_train(self):
        ''' Take lists of xvals (energies), get associated yvals and 
        save to backfitpts (along with filename) ... called by GUIfitter button
        pandas df with columns containing original x and y, removed x and y and 
        added x and y
        TODO FINISH ME
        '''
        print('Starting EDXfile.save_train()')
        # Ensure that some changes were made 
        print('length of removed and added points is', len(self.removedpts), 'and', len(self.addedpts))
        if not self.addedpts and not self.removedpts:
            print('no changes recorded')
            return
        # Single row dataframe with correct cols for this modified EDXfile
        newser=pd.Series()
        newser['Filename']=self.filename
        newser['Beamkv']=self.beamkv
        newser['Deadfraction']=self.deadfraction
        newser['Xvals']=self.origbackfitpts
        newser['Xrem']=self.removedpts # should be list of ints
        newser['Xadd']=self.addedpts
        # Get Yvals for these indices (not energy vals which are 0.01)
        newser['Yvals']=self.EDXdf[self.EDXdf.index.isin(self.origbackfitpts)]['Counts'].tolist()
        newser['Yrem']=self.EDXdf[self.EDXdf.index.isin(self.removedpts)]['Counts'].tolist()
        newser['Yadd']=self.EDXdf[self.EDXdf.index.isin(self.addedpts)]['Counts'].tolist()
        self.EDXdataset.update_training(self.filename, newser)
                
    def save_backfits(self):
        ''' Save changes to Backfitparamslog file (incl. fit ranges, backfitpts, etc.
        just put here and not in EDXdataset as changes are made on per EDXfile basis'''
        # TODO peak integrations will also be affected 
        match=self.EDXdataset.Backfitlog[self.EDXdataset.Backfitlog['Filename']==self.filename]
        # Save altered subset
        for i, evrange in enumerate(self.fitranges):
            thisind=match.index[i]
            # Set values on same index in EDXdataset backfitlog
            self.EDXdataset.Backfitlog=self.EDXdataset.Backfitlog.set_value(thisind, 
                'Fitrange', '-'.join([str(i) for i in self.fitranges[i]]))
            self.EDXdataset.Backfitlog=self.EDXdataset.Backfitlog.set_value(thisind, 'Backfitpts', self.backfitpts[i])
            if self.fitorders[i]==1:
                fittype='linear'
            elif self.fitorders[i]==2:
                fittype='parabola'
            elif self.fitorders[i]==3:
                fittype='cubic'
            self.EDXdataset.Backfitlog=self.EDXdataset.Backfitlog.set_value(thisind, 'Fittype', fittype)
            self.EDXdataset.Backfitlog=self.EDXdataset.Backfitlog.set_value(thisind, 'A', self.backfitparams[i][0])
            self.EDXdataset.Backfitlog=self.EDXdataset.Backfitlog.set_value(thisind, 'B', self.backfitparams[i][1])
            self.EDXdataset.Backfitlog=self.EDXdataset.Backfitlog.set_value(thisind, 'C', self.backfitparams[i][2])
            self.EDXdataset.Backfitlog=self.EDXdataset.Backfitlog.set_value(thisind, 'D', self.backfitparams[i][3])
        # Save entire modified csv file 
        self.EDXdataset.Backfitlog.to_csv('Backfitparamslog.csv', index=False)
        # On backfit save, redo gaussfit, integ counts calc and save peakfitlog and integlog
        self.fitpeaks()
        print('Modified Backfitlog and updated peak and integlogs saved.')
    
    def fitpeaks(self):
        ''' Rerun of fitpeaks after background fit adjustments 
        EDXdf contains loaded edxfile (Energy, Counts, Backfit, Subdata, Gauss
        run automatically if new backfits are saved  '''
        
        # gets Elemdata for this EDXfile from EDXdataset (Elements not used)
        Elements, self.Elemdata=self.EDXdataset.findelemregions(self.filename)

        Peakfits=pd.DataFrame(columns=self.EDXdataset.Peakfitlog.columns) # blank df for this spectrum's peak fits
        Integresults=pd.DataFrame(columns=self.EDXdataset.Integlog.columns) # blank df for this spectrum's integration results

        for i, (elem, idealindex, maxshift, halfwidth, kfact, errkfact, mass) in enumerate(self.Elemdata):
            Peakfitrow=pd.DataFrame(index=np.arange(0,1),columns=self.EDXdataset.Peakfitlog.columns) 
            Integresultrow=pd.DataFrame(index=np.arange(0,1),columns=self.EDXdataset.Integlog.columns)
            # linear fit below this elem's peak (shifts and adjustments already made)
            fitregion=self.EDXdf[idealindex-halfwidth-5:idealindex+halfwidth+6]
            if fitregion.empty==True: # skip if no data present (peak out of range problem)
                continue
            # Gaussian fit of subtracted data peaks > 50 cnts
            if fitregion['Subdata'].max()>50: # add flag to skip gaussian
                fitregion, fitparams, rsquared, ier = self.fitgauss(fitregion, halfwidth, 
                                            elem)
                # save Gaussian peaks as separate column by default
                if 'Gauss' not in self.EDXdf.columns: # add col if not already present                
                    self.EDXdf['Gauss']='' # add blank col for gaussian fit if not present
                # Copies new Gauss fit back to EDXdf
                self.EDXdf.loc[fitregion.index,fitregion.columns]=fitregion
                          # determination of peak shift 
            # If gaussian fit is successful set center integration channel to index nearest xc
            # ier flag of 1,2,3,4 if fit succeeds but rsquared threshold is better
                if rsquared!='n/a': # somewhat successful gaussian fit 
                    if rsquared>0.4:
                        xc=fitparams[0] # center of gaussian fit in keV
                        centerindex=int((xc+.01)*100)
                        shift= centerindex- idealindex # energy shift in channels
                        if abs(shift)>maxshift: # maxshift is element specific maximum move of integration window
                            # common problem with weak peaks... only use print for troubleshoot
                            # print('Warning: Gaussian shift of ', str(shift), ' channels indicated for ', elem, ' in ', EDXfileName)
                            if shift>0: # keep peak shift the same but only allow 3 channel shift in integration window
                                centerindex=idealindex+maxshift # set to max shift
                            else:
                                centerindex=idealindex-maxshift
                # TODO Maybe a better way of setting maximum allowable shift        
                    else: 
                        # common problem with mass fitting so skip print report
                        # print('Low quality gaussian fit for ', elem, ' in ', EDXfileName)
                        centerindex=idealindex # already stores index number of central peak (ideal - sm-diff shift value)
                        shift='n/a'
                # Write gaussian fit params to peakfit (eventually copied to peakfitlog)
                        
                else: # Fit attempted but failed result
                    print ('Fit attempted but result failed for ', elem, ' in ', self.filename)
                    fitparams=['n/a','n/a','n/a','n/a']            
                    rsquared='n/a'
                    
            else: # indication of failed Gaussian fit (use prior knowledge of peak position)
                # common problem with weak peaks... only use print for troubleshoot
                # print('Skip gaussian fit of tiny ', elem, ' peak in ', EDXfileName)
                # set center integration channel to value passed by integpeak 
                # this is ideal energy value but adjusted by shift found using smooth-diff quant method
                centerindex=idealindex # already stores index number of central peak (ideal - sm-diff shift value)
                shift='n/a'
                fitparams=['n/a','n/a','n/a','n/a']            
                rsquared='n/a'
            # Perform integration over peak center channel + integwidth on either side 
            EDXpeak=self.EDXdf[centerindex-halfwidth:centerindex+halfwidth+1]
            # TODO problem... subdata isn't updated w/ correct subtracted background
            integcounts=EDXpeak['Subdata'].sum() # get counts sum 
            backgroundcnts=EDXpeak['Backfit'].sum() # sum counts over identical width in background fit
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
            Integresultrow.iloc[0]['Rawcounts']=EDXpeak['Counts'].sum() 
            Integresultrow.iloc[0]['Backcounts']=backgroundcnts
            Integresultrow.iloc[0]['Subtractedcounts']=integcounts
            # Adjusted counts must be determined later for pathological overlaps
            # 2 sigma err due to counting statistics
            Integresultrow.iloc[0]['% err']=round(2/np.sqrt(integcounts),3)
            Integresultrow.iloc[0]['Significance']=round(integcounts/(np.sqrt(backgroundcnts)),3)
    		# TODO add 2/sqrt(n) calc of associated percent error (also can calculate later)
            Integresultrow.iloc[0]['Correctedcounts']=integcounts*kfact/mass
            # Calculated combined error for 2sig counting stats + loaded k-factor error
            comberr=np.sqrt(errkfact**2+(2/np.sqrt(integcounts))**2)
            # calculate error in Correctedcounts for given elemental peak
            Integresultrow.iloc[0]['Errcorrcnts']=(integcounts*kfact/mass)*comberr
            Integresultrow.iloc[0]['Kfact']=kfact
            Integresultrow.iloc[0]['Fullwidth']=2*halfwidth
            Integresults=pd.concat([Integresults,Integresultrow], ignore_index=True)
            
        # assign params that are common to this spectrum (all elemental peaks)
        for index,row in Peakfits.iterrows(): 
            # need to replace logmatch w/ correct row from EDXlog
            Peakfits.loc[index]['Filenumber']=self.row.Filenumber   
            Peakfits.loc[index]['Basename']=self.row.Basename
            Peakfits.loc[index]['Filename']=self.row.Filename
            Peakfits.loc[index]['Point']=self.row.Point
            Peakfits.loc[index]['Filepath']=self.row.FilePath
            Peakfits.loc[index]['Sample']=self.row.Sample
            Peakfits.loc[index]['Comments']=self.row.Comments
        for index,row in Integresults.iterrows(): # assign
            Integresults.loc[index]['Filenumber']=self.row.Filenumber   
            Integresults.loc[index]['Filename']=self.row.Filename
            Integresults.loc[index]['Basename']=self.row.Basename
            Integresults.loc[index]['Point']=self.row.Point
            Integresults.loc[index]['Filepath']=self.row.FilePath
            Integresults.loc[index]['Sample']=self.row.Sample
            Integresults.loc[index]['Comments']=self.row.Comments
        Peakfits=Peakfits[self.EDXdataset.Peakfitlog.columns] # put back in original order
        Integresults=Integresults[self.EDXdataset.Integlog.columns] # put back in original order
        # now write/replace values in EDXdataset.Peakfitlog and Integlog
        self.EDXdataset.updatepeaks(self.filename, Peakfits)
        self.EDXdataset.updateinteg(self.filename, Integresults)
        
    def fitgauss(self, df, halfwidth, elem):
        ''' Gaussian fit of direct peaks (pass EDXfile just around peaks region
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
            print('Gaussian fitting error for', elem, ' peak in file ', self.filename)
            fitparams=('n/a','n/a','n/a','n/a') # return all n/a
            rsquared='n/a'
            ier='n/a'
            return df, fitparams, rsquared, ier
        if 'Gauss' not in df:
            df['Gauss']='' # add col for gaussian fit
        for index,row in df.iterrows():
            xval=df.loc[index]['Energy']
            yval=fitparams[3]+fitparams[2]/(fitparams[1]*np.sqrt(2*np.pi))*np.exp(-((xval-fitparams[0])**2/(2*fitparams[1]**2)))
            df=df.set_value(index,'Gauss',yval)
        return df, fitparams, rsquared, ier
    
class EDXdataset():
    ''' loads all dataframes with EDX parameters from current project folder '''
    def __init__(self, directory, **kwargs):
        # open files 
        os.chdir(directory)
        self.EDXlog, self.Backfitlog, self.Integlog, self.Peakfitlog, self.EDXquantparams, self.Interferences, self.Backtraining = self.open_main_files()
		# Filenumbers are often not unique (using different basenames
        # self.filelist=np.ndarray.tolist(self.EDXlog.Filenumber.unique())
        self.numfiles=len(self.EDXlog)
        self.spectype= None # SEM or TEM 
        # Autoload first file
        print(str(self.numfiles),' loaded from EDXdataset.')

    def updatepeaks(self, EDXfilename, Peakfits):
        ''' After fitpeak rerun (in EDXfile), update and save Peakfits log for  
        EDXdataset 
        usually called at end of EDXfile fitpeaks method '''
        # Remove old values
        self.Peakfitlog=self.Peakfitlog[self.Peakfitlog['Filename']!=EDXfilename]
        self.Peakfitlog=self.Peakfitlog.append(Peakfits, ignore_index=True)
        # Save to file
        self.Peakfitlog.to_csv('Peakfitlog.csv', index=False)
        print('Updated peakfitlog for', EDXfilename)

    def updateinteg(self, EDXfilename, Integresults):
        ''' After fitpeak rerun (in EDXfile), update and save Integresults 
        for EDXdataset 
        usually called at end of EDXfile fitpeaks method '''
        # Remove old values
        self.Integlog=self.Integlog[self.Integlog['Filename']!=EDXfilename]
        self.Integlog=self.Integlog.append(Integresults, ignore_index=True)
        # Save to file
        self.Integlog.to_csv('Integquantlog.csv', index=False)
        print('Updated integquant for', EDXfilename)

    def update_training(self, EDXfilename, trainser):
        ''' After interactive modification of background fitting  points, update 
        and save background point training file (called by GUIrefitter backtraining
        button, through EDXfile '''
        # Remove old values
        self.Backtraining=self.Backtraining[self.Backtraining['Filename']!=EDXfilename]
        self.Backtraining=self.Backtraining.append(trainser, ignore_index=True)
        # Save to file
        self.Backtraining.to_csv('Backfit_training.csv', index=False)
        print('Updated background fitting training file with ', EDXfilename)
        
    def findelemregions(self, EDXfilename):
        ''' For active edx file, get prior elements list and then detailed 
        element data for each in elements list '''
        thislist=self.Integlog[self.Integlog['Filename']==EDXfilename]
        Elements=np.ndarray.tolist(thislist.Element.unique())
        Elemdata=[] # initialize as empty list
        for i, elem in enumerate(Elements):
            try:
                # find row in AESquantparams for this element
                thiselemdata=self.EDXquantparams[(self.EDXquantparams['element']==elem)]
                thiselemdata=thiselemdata.squeeze() # series with this elements params
                
                # integ peak position value is relative to negpeak in smooth-diff (i.e. -5 is 5 eV below ideal negpeak)
                idealindex=int((thiselemdata.energy+.01)*100) # ideal index value of EDX-EDX peak from energy in keV
                kfact=thiselemdata.kfactor # typical sensitivity k-factor associated with element for integration
                errkfact=thiselemdata.errkfact 
                mass=thiselemdata.mass
                maxshift=int(thiselemdata.maxshift) # on indices so must be int
                # full peak width in keV from EDXquantparams (usually 0.15keV or 15 channels at 0.1eV/chan)
                # integration width in channels for direct integration for this element
                width=int(((thiselemdata.fullwidth*100)-1)/2) 
                # total # of channels in AESquantparams but include n-1/2 channels on either side of peak center (usually width is 8 channels)
                #Elemdata is a list (of length number of elements) containing length5 tuples
                elemtuple=(elem, idealindex, maxshift, width, kfact, errkfact, mass) # add tuple with info for this element
                Elemdata.append(elemtuple) # now contains proper limits on fitting regions 
            except:
                print('Quant parameters not properly loaded for', elem)
        return Elements, Elemdata
            
    
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
        if os.path.exists('Backfit_training.csv'):
            # open file for interactive background fit training
            Backtraining=pd.read_csv('Backfit_training.csv', encoding='cp437')
        else:
            Backtraining=pd.DataFrame(columns=['Filename', 'Beamkv','Deadfraction',
                'Xvals','Yvals','Xrem','Yrem','Xadd','Yadd'])
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
            self.spectype='TEM'
        else:
            print(EDXlog['Beamkv'].max(),'keV SEM spectra and quant params loaded.')
            EDXquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEMquantparams.csv', encoding='utf-8')
            Interferences=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\EDX\\SEM_interferences.csv', encoding='utf-8')
            self.spectype='SEM'
        return EDXlog, Backfitlog, Integlog, Peakfitlog, EDXquantparams, Interferences, Backtraining
