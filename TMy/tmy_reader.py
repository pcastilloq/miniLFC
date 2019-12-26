# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 09:16:52 2019

@author: pablo
"""

import pvlib 
import matplotlib.pyplot as plt

#%%
tmy_scl = pvlib.iotools.read_tmy3('OK_SCL.csv')

for i in range(tmy_scl[0].shape[0]):
    
    DNI = tmy_scl[0]['DNI'][i]
    T_amb = tmy_scl[0]['DryBulb'][i]
    v_wind = tmy_scl[0]['Wspd'][i]
    
    
