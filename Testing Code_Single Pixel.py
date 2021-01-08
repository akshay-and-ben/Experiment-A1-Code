# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 16:48:31 2021

@author: rr3
"""

import numpy as np
import scipy as sp
from scipy.special import wofz
import pylab as pl
import matplotlib.pyplot as plt
import scipy.optimize as opt
from astropy.io import fits
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry

#%%
#CREATE TESTING AREA
size = (500,500)
data = np.zeros(size)
data[:,:]=3400
#Single pixel
data[200,100]=30000

#%%
""" Locating the positions of each galaxy and appending these to a list
for later analysis """
y = []
x = []
maxval = []
stray_x = []
stray_y = []
max_val = np.amax(data)

# A galaxy is defined to be an object with maximum pixel brightness >= 3500
while max_val >= 3600:
    # Locates the maximum value within the dataset, and notes
    # the x and y position
    a = np.where(data == max_val)
    ypos = a[0][0]  
    xpos = a[1][0]
   
    # ensures that a proper galaxy is detected, as opposed to a single
    # bright pixel (due to a stray electron hitting the CCD)
    if any(i >= 3600 for i in data[ypos-1:ypos+2, xpos-1:xpos+2][0]) or \
    all(i >= 3600 for i in data[ypos-1:ypos+2, xpos-1:xpos+2][1]) or \
    any(i >= 3600 for i in data[ypos-1:ypos+2, xpos-1:xpos+2][2]):
       
        # Appends the x and y position of the galaxy to a list, and its maximum
        # value
        x.append(xpos)
        y.append(ypos)
        maxval.append(max_val)
       
        # After a galaxy is identified, a square region of radius 14 is masked,
        # so as not to detect the galaxy again.
        data[ypos-13:ypos+14, xpos-13:xpos+14] = 0
        max_val = np.amax(data)
        print("Galaxy Detected")
       
    else:
        # Masks and flags up single pixel events.
        data[ypos, xpos] = 0
        stray_x.append(xpos)
        stray_y.append(ypos)
       
        max_val = np.amax(data)
        print("Single pixel caused by stray electron hitting CCD")

positions = []
for j in range(len(x)):
    positions.append((x[j], y[j]))


