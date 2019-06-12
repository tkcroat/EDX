# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 17:11:04 2016
k-factor functions
@author: tkc
"""

def computecomposition(df, type='AES'):
    ''' Pass AES or SEM counts and recompute/return dataframe  
    types are AES or SEM'''
    # load df with current k-factors
    

def calculateAESbasis(df,AESkfactors):
    
    # load AES quant params
    AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\AESquantparams.csv', encoding='utf-8')
    