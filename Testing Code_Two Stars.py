# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 17:32:40 2021

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

#LOAD DATA
hdulist = fits.open('A1_mosaic.fits')
dataset = hdulist[0].data

dp = np.flipud(dataset) #Orient data to match ds9 file

data = dp[900:1060, 340:400]
data_unmod = data.copy()
sample = data

#%%
""" Locating the positions of each galaxy and appending these to a list
for later analysis """
y = []
x = []
maxval = []
stray_x = []
stray_y = []
max_val = np.amax(sample)

# A galaxy is defined to be an object with maximum pixel brightness >= 3500
while max_val >= 3600:
    # Locates the maximum value within the dataset, and notes
    # the x and y position
    a = np.where(sample == max_val)
    ypos = a[0][0]  
    xpos = a[1][0]
   
    # ensures that a proper galaxy is detected, as opposed to a single
    # bright pixel (due to a stray electron hitting the CCD)
    if any(i >= 3600 for i in sample[ypos-1:ypos+2, xpos-1:xpos+2][0]) or \
    all(i >= 3600 for i in sample[ypos-1:ypos+2, xpos-1:xpos+2][1]) or \
    any(i >= 3600 for i in sample[ypos-1:ypos+2, xpos-1:xpos+2][2]):
       
        # Appends the x and y position of the galaxy to a list, and its maximum
        # value
        x.append(xpos)
        y.append(ypos)
        maxval.append(max_val)
       
        # After a galaxy is identified, a square region of radius 14 is masked,
        # so as not to detect the galaxy again.
        sample[ypos-13:ypos+14, xpos-13:xpos+14] = 0
        max_val = np.amax(sample)
        print("Galaxy Detected")
       
    else:
        # Masks and flags up single pixel events.
        sample[ypos, xpos] = 0
        stray_x.append(xpos)
        stray_y.append(ypos)
       
        max_val = np.amax(sample)
        print("Single pixel caused by stray electron hitting CCD")

positions = []
for j in range(len(x)):
    positions.append((x[j], y[j]))

#%%
#FIXED APERTURE

#Create circular aperture
aperture = CircularAperture(positions, r=5)

#Create circular annulus
annulus_aperture = CircularAnnulus(positions, r_in=10., r_out=12.)

apers = [aperture, annulus_aperture]

#Perform aperture photometry
phot_table = aperture_photometry(data_unmod, apers)

for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output
print(phot_table)

#Calculate mean background value per pixel in annulus
bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area

#Calulate total mean background in aperture
bkg_sum = bkg_mean * aperture.area

#Calculate residual aperture pixel count
final_sum = phot_table['aperture_sum_0'] - bkg_sum
phot_table['residual_aperture_sum'] = final_sum

phot_table['residual_aperture_sum'].info.format = '%.8g'  # for consistent table output
print(phot_table['residual_aperture_sum']) 
#%%   
#VARIABLE APERTURE
"""First stage of variable aperture - determine radius at which to perform 
the aperture photometry with"""

#List of aperture radii to iterate over
radii_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
rad_arr = np.array(radii_list)

#Establish apertures at radii given in radii_list
apertures_test = [CircularAperture(positions, r=r) for r in radii_list] 

#Initial background annulus
annuli = CircularAnnulus(positions, r_in=12, r_out=14) 

#Perform aperture photometry for all given radii
phot_table = aperture_photometry(data_unmod, apertures_test) 
back_table = aperture_photometry(data_unmod,annuli)

for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output

for col in back_table.colnames:
    back_table[col].info.format = '%.8g'  # for consistent table output
   
bkg_mean = back_table['aperture_sum'] / annuli.area

areas = np.zeros_like(apertures_test) 
for j in range(len(apertures_test)):
    #Find areas of respective apertures
    areas[j]=apertures_test[j].area 

bkg_sum_test = bkg_mean*areas

# get headers
headers = np.array(phot_table.columns) 

#extract data from Astropy table
cols=phot_table.columns 
mylist = [cols[k] for k in range(len(cols))]

#remove unnecessary information from Astropy table
cutoff = 3
del mylist[:cutoff] 

#%%
extract_aperture = np.zeros_like(mylist)
for i in range(len(mylist)):
    #convert Astropy table into NumPy array
    extract_aperture[i] = mylist[i] 

extract_array = extract_aperture[:,:]

"""Method to find cutoff radius by determining intensity of source as a 
function of distance from centre"""

#Determine average pixel intensity as distance from centre increases
count_decrease = (extract_aperture - bkg_sum_test)/areas
count_array = np.stack(count_decrease).astype(None) 
percs = np.zeros_like(count_array)

#Determine central pixel count
peak_val = np.array(maxval) 
num_rows = len(count_array[0])
for o in range(num_rows):
    #Determine intensity of light emitted as fraction of central value
    percs[:,o]=count_array[:,o]/peak_val[o] 

percs_cols = len(percs[0])
cut_idx = np.zeros(percs_cols)

for p in range(percs_cols):
    #Establish percentage cutoff for size of final aperture radius 
    cut_idx[p] = np.argmax(percs[:,p]< 0.3) 
cut_idx = np.array(cut_idx)
cut_idx = cut_idx.astype(int)

#Find maximum radius that aperture should extend to
max_rad = rad_arr[cut_idx]

#Establish minimum aperture radius
idx_min = np.where(max_rad < 4) 
max_rad[idx_min] = 4
#%%
"""Second part of variable aperture - perform the aperture photometry at
the previously-determined radii for each source respectively"""

apertures_list = np.zeros_like(positions)
annuli_list = np.zeros_like(positions)

for q in range(len(max_rad)):
    #Perform the aperture photometry at radius determined for each source
    aperture = CircularAperture(positions[q], r=max_rad[q]) 
    phot_table_apers = aperture_photometry(data_unmod, aperture)
    apertures_list[q] = phot_table_apers['aperture_sum']
    #Find total background in annulus for each source
    annul = CircularAnnulus(positions[q], r_in=10, r_out=12) 
    phot_table_annul = aperture_photometry(data_unmod,annul)
    annuli_list[q] = phot_table_annul['aperture_sum']
    
aperture_list_final = apertures_list[:,0]
annuli_list_final = annuli_list[:,0]

#Find mean background value around each source
bckg_mean_list = annuli_list_final/annul.area
 
bckg_vals = np.zeros_like(positions)
for r in range(len(max_rad)):
    #Find total background in a particular aperture
    bckg_vals[r] = (np.pi*((max_rad[r])**2))*bckg_mean_list[r] 
bckg_vals_final = bckg_vals[:,0]

#Determine corrected pixel count
residuals_list = aperture_list_final - bckg_vals_final
print(residuals_list)
print(max_rad)