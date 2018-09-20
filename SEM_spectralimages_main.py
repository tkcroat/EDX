# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:11:29 2016

@author: tkc
"""

import pandas as pd
import sys

if 'C:\\Users\\tkc\\Documents\\Python_Scripts' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts')
    
import SEMEDX_batch_import_functions as SEMimport
import SEMEDX_plot_functions as SEMplot
import SEM_SI_functions as SI
from skimage.feature import blob_dog, blob_log, blob_doh # for image segmentation

#%%
# Read or create spectral images log for this folder
NSScsvparams=SEMimport.batchNSSreader()  # option 1 which recreates param log after files are added
NSScsvparams=pd.read_csv('NSScsvparams.csv', encoding='cp437') # option 2 ... just load existing parameters file

# NSScsvparams are extracted from si files.. stagepos generage .. .use same name for tif and SI
# Find stage coords for NSS subset (using matching string and APF file name)
NSScsvparams=getstagefromapf(NSScsvparams, '3x3Simap', 'JEOLSEMarr3by3at3000overlap20.apf')
NSScsvparams.to_csv('NSScsvparams.csv', index=False)



#load/create set of x-ray image maps
Elements=['Si','Mg','O','Pt','Ar']
Elements=['Si','Mg','O','Pt','Grey']
Elements=['all'] # retrieves all available elements maps for given file

# Batch plotting of maps from all in NSS param list 
SEMplot.reportmaps(NSScsvparams, Elements, PDFreport='NSSmaps_report.pdf')

# Select single x-ray map file set (change index #s from list)
thismaprow=NSScsvparams.loc[[0]]

elementmaps=SEMimport.getelemmaps(thismaprow, Elements)
elementmaps=getelemmaps(thismaprow, Elements)

# Plot set of element maps from single spectral image
SEMplot.plotmaps(elementmaps, thismaprow, savename='')

# Batch plot/report of all available sets of x-ray maps (extracted from si files)
SEMplot.reportmaps(NSScsvparams, Elements, PDFreport='NSSmaps_report.pdf')

# Return numpy arrays with single element map
Omap=SEMimport.getsinglemap(elementmaps,'O')
Simap=SEMimport.getsinglemap(elementmaps,'Si')
SE=SEMimport.getsinglemap(elementmaps,'Grey')

# Plot with single element map
plt.ion()
plt.figimage(Omap, cmap='hot')
plt.imshow(SE, cmap='gray')
plt.imshow(Omap, cmap='hot')

# Construct larger x-ray map using stagepositions w/ included pixel shifts
Stagepositions=pd.read_csv('Stagepositions_mont3x3overlap20.csv') 
thismap=SI.assemblexraymap(Stagepositions, elem='Si K', pixsizes=[512,384])
plt.imshow(thismap)
# Assemble stack of x-ray image maps


