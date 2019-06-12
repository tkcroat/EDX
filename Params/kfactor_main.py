# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 12:59:54 2016
k-factor determination from C2010W crater data
@author: tkc
"""
import os, sys
import pandas as pd
import numpy as np
#%%
os.chdir('C:\\Users\\tkc\\Documents\\Python_Scripts\\kfactors')
# SEM basis > 20, AESbasis>200, no skips (SIMS damage), 
bestAEScounts=pd.read_csv('C2010W_bestAES_counts.csv', encoding='cp437') 

bestSEMcounts=pd.read_csv('C2010W_bestSEM_counts.csv', encoding='cp437') 
#%%
# sort dfs by Crater_name
bestAEScounts=bestAEScounts.sort_values(['Crater_name'], ascending=True)
bestSEMcounts=bestSEMcounts.sort_values(['Crater_name'], ascending=True)

#find subset with quality SEM and Auger spectra
myarray=bestAEScounts.Crater_name.unique()
AESlist=np.ndarray.tolist(myarray)
AESset=set(AESlist)
myarray=bestSEMcounts.Crater_name.unique()
SEMlist=np.ndarray.tolist(myarray)
SEMset=set(SEMlist)
set(SEMset) ^ set([3,4,5])
commonset = SEMset ^ AESset
commonset = SEMset & AESset
commonlist=list(commonset)
commonlist.sort(key=str.lower) # sort in alphabetical order

# Return subsets with both AES and SEM data
bestAES=bestAEScounts[bestAEScounts['Crater_name'].isin(commonlist)]
bestSEM=bestSEMcounts[bestSEMcounts['Crater_name'].isin(commonlist)]

# find # of craters with duplicated spectra

# keep best spectrum 
bestAES=bestAES.sort_values(['Crater_name','Grade'], ascending=True)

df=df.sort_values(['Crater_name','Grade'], ascending=True) 

bestAESnodup=bestAES.drop_duplicates(subset=['Crater_name'], keep='first')
bestSEMnodup=bestSEM.drop_duplicates(subset=['Crater_name'], keep='first')

# 
mycols=bestAEScounts.columns.tolist()
for i, colname in enumerate(mycols):
    if 'width' in colname:
        bestAEScounts.drop(colname,axis=1,inplace=True)

    