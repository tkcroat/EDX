# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 13:27:36 2017

@author: tkc
"""

import hyperspy.api as hs

spec=hs.load(filelist[0], signal_type='EDS_SEM')
spec.metadata.Acquisition_instrument.SEM
spec.set_microscope_parameters() # should start GUI

spec.metadata # dictionary tree browser associated with spectrum
spec.plot()

