# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 05:19:15 2021

@author: rr3
"""

"""
A1 Astronomical Imaging Lab:
Authors: Ben Amroota; Akshay Robert
ICL Shortcodes: bra18, ar6918

The code completes the following actions:
    1. Read the datafile
    2. Manually mask blooming stars, and remove edges of image due to
    increased noise
    3. Identifies galaxies present within the image and notes the positions
    4. Uses the positions in aperture photometry using the Astropy module to
    determine the total flux through a variable aperture, the local background
    for each galaxy, and return a list of corrected aperture fluxes.
    5. Plots the data in the form log(N(<m)) vs m, and fits a linear fit to it
    using the scipy.optimize.curvefit module
    6. Writes the data to a csv file for external analysis.
   
"""
#%%
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import scipy.optimize as opt
from astropy.io import fits
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
#%%
""" Read FITS datafile """
hdulist = fits.open('A1_mosaic.fits')
dataset = hdulist[0].data
dp = np.flipud(dataset)
ds = dp.copy()
#%%
""" Mask blooming effects """
dp[480:540, 540:580] = 0
dp[1150:1430, 720:835] = 0
dp[1760:1920, 950:1000] = 0
dp[800:910, 2110:2160] = 0
dp[2450:2480, 1160:1240] = 0
dp[3150:3240, 2060:2120] = 0
dp[1660:4611, 1420:1500] = 0
dp[1620:1660, 1423:1450] = 0
dp[0:1620, 1400:1470] = 0
dp[2266:2333, 2100:2160] = 0
dp[1150:1250, 2425:2500] = 0
dp[2240:2500, 860:940] = 0
dp[4600:4611, 969:1720] = 0
dp[4150:4600, 1240:1420] = 0
dp[4460:4520, 1500:1550] = 0
dp[4550:4600, 1630:1650] = 0
dp[4280:4305, 1000:1250] = 0
dp[4170:4190, 1499:1650] = 0
dp[4250:4300 ,1499:1705] = 0
dp[1290:1510, 1325:1550] = 0
dp[4146:4192, 1018:1250] = 0
dp[1160:1230, 2450:2482] = 0
dp[540:620, 1420:1500] = 0
dp[2625:2700, 2450:2550] = 0
#%%
""" Masking noise defects present in the edges of image;
removing a width of 113 pixels from each side """
dp[:,2457:2570] = 0
dp[:,0:113] = 0
dp[4498:4611,:] = 0
dp[0:113,:] = 0

#Choose a sample size of the overall image 'dp' for analysis - 
#currently set sample to the whole image for full analysis
sample = dp
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
while max_val >= 3500:
    # Locates the maximum value within the dataset, and notes
    # the x and y position
    a = np.where(sample == max_val)
    ypos = a[0][0]  
    xpos = a[1][0]
   
    # ensures that a proper galaxy is detected, as opposed to a single
    # bright pixel (due to a stray electron hitting the CCD)
    if any(i >= 3500 for i in sample[ypos-1:ypos+2, xpos-1:xpos+2][0]) or \
    all(i >= 3500 for i in sample[ypos-1:ypos+2, xpos-1:xpos+2][1]) or \
    any(i >= 3500 for i in sample[ypos-1:ypos+2, xpos-1:xpos+2][2]):
       
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
""" HEAT MAP PLOT """
'''
pl.figure()
pl.imshow(dp)
pl.show()'''
#%%
""" First stage of variable aperture - determine radius at which to perform
the aperture photometry with """
# List of aperture radii to iterate over
radii_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
rad_arr = np.array(radii_list)
# Establish apertures at radii given in radii_list
apertures_initial = [CircularAperture(positions, r=r) for r in radii_list]

# Initial background annulus
annuli = CircularAnnulus(positions, r_in=12, r_out=14)

# Perform aperture photometry for all given radii
phot_table = aperture_photometry(ds, apertures_initial)
back_table = aperture_photometry(ds,annuli)

for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output
for col in back_table.colnames:
    back_table[col].info.format = '%.8g'  # for consistent table output
   
bkg_mean = back_table['aperture_sum'] / annuli.area

# Find areas of respective apertures
areas = np.zeros_like(apertures_initial)
for j in range(len(apertures_initial)):
    areas[j]=apertures_initial[j].area

bkg_sum_initial = bkg_mean*areas

# Get headers
headers = np.array(phot_table.columns)

# Extract data from Astropy table and remove unnecessary information from
# Astropy table
cols=phot_table.columns
mylist = [cols[k] for k in range(len(cols))]
cutoff = 3
del mylist[:cutoff]
extract_aperture = np.zeros_like(mylist)

# Convert Astropy table into NumPy array
for i in range(len(mylist)):
    extract_aperture[i] = mylist[i]
#%%
""" Method to find cutoff radius by determining intensity of source as a
function of distance from centre """

count_decrease = (extract_aperture - bkg_sum_initial)/areas
# Determine average pixel intensity as distance from centre increases
count_array = np.stack(count_decrease).astype(None)
percs = np.zeros_like(count_array)
# Determine central pixel count
peak_val = np.array(maxval)
num_rows = len(count_array[0])

# Determine intensity of light emitted as fraction of central pixel value
for o in range(num_rows):
    percs[:,o]=count_array[:,o]/peak_val[o]
percs_cols = len(percs[0])
cut_idx = np.zeros(percs_cols)

# Establish cutoff of how large aperture radius can get before selection
for p in range(percs_cols):
    # if flux denisty of light drops below 0.1 % of the central pixel value
    # then choose this as the radius for aperture photometry
    cut_idx[p] = np.argmax(percs[:,p]< 0.001)
cut_idx = np.array(cut_idx)
cut_idx = cut_idx.astype(int)
# Find maximum radius that aperture should extend to
max_rad = rad_arr[cut_idx]

# Establish minimum aperture radius
idx_min = np.where(max_rad < 6)
max_rad[idx_min] = 6

#%%
""" Second part of variable aperture - perform the aperture photometry at
the previously-determined radii for each source respectively """

apertures_list = np.zeros_like(positions)
annuli_list = np.zeros_like(positions)
 
#Perform the aperture photometry at the determined radius
for q in range(len(max_rad)):
    aperture = CircularAperture(positions[q], r=max_rad[q])
    phot_table_apers = aperture_photometry(ds, aperture)
    apertures_list[q] = phot_table_apers['aperture_sum']
    #Find background in fixed annulus of 2 around the variable aperture
    annul = CircularAnnulus(positions[q], r_in=max_rad[q]+1, r_out=max_rad[q]+3)
    phot_table_annul = aperture_photometry(ds,annul)
    annuli_list[q] = phot_table_annul['aperture_sum']

aperture_list_final = apertures_list[:,0]
annuli_list_final = annuli_list[:,0]

# Find mean background value around each source
bckg_mean_list = annuli_list_final/annul.area
bckg_vals = np.zeros_like(positions)

#Find total background in a particular aperture
for r in range(len(max_rad)):
    bckg_vals[r] = (np.pi*((max_rad[r])**2))*bckg_mean_list[r]

bckg_vals_final = bckg_vals[:,0]
# Determine corrected pixel count
residuals_list = aperture_list_final - bckg_vals_final

# Error determined from Poisson statistics
err_bkg_sum = np.sqrt(annuli_list_final)
err_final_sum = np.sqrt((annuli_list_final)+(aperture_list_final))
err_residuals = err_final_sum
#%%
""" From FITS headers read magnitude calibration details, and
calculate the magnitude of each galaxy and error associated.
Values are moderated to remove nans and infinities. """
hdr = hdulist[0].header

# Find calibrated magnitude i.e. instrumental zero-point value
mag_cal = hdr['MAGZPT']

# Find error in calibrated magnitude
err_cal = hdr['MAGZRR']

# Convert number counts to instrumental magnitudes
mag = -2.5*np.log10(residuals_list)
err_mag = 0.434*(err_residuals/residuals_list)

# Number counts is the value after cataloguing the stars in photometry section of lab script
m = mag_cal + mag
err_m = np.sqrt(np.square(err_cal) + np.square(err_mag))

# Removing infinity values
m_a = np.delete(m, np.where(np.isinf(m)))
err_m_a = np.delete(err_m, np.where(np.isinf(m)))

# Removing nans from magnitude array
for i in range(len(m_a)):
    if np.isnan(m_a[i]) == True:
        m_mod = m_a[np.logical_not(np.isnan(m_a))]
        err_m_mod = err_m_a[np.logical_not(np.isnan(m_a))]
#%%
""" Calculating the number count N and the logarithm from the magnitudes """
m_array = np.linspace(min(m_mod),max(m_mod),len(m_mod))
N = np.array([0]*len(m_array))
err_N = N.copy()

# Counts the number of m_mod values less than the ith value of m_array; this
# is the number count N
for i in range(len(m_mod)):
    N[i] = np.sum(m_mod < m_array[i])

# The error on N is described with Poisson statistics, hence is np.sqrt(N)
err_N = np.sqrt(N)
logN = np.log10(N)
#%%
""" Editing errorbars to remove unphysical errorbars, and to remove infinities
in logN """

m_array_mod = m_array

# Removes all x errorbars deemed too large
m_array_mod_mod = np.delete(m_array_mod, np.where(err_m_mod >1))
logN_mod = np.delete(logN, np.where(err_m_mod >1))
err_N_mod = np.delete(err_N, np.where(err_m_mod >1))
N_mod = np.delete(N, np.where(err_m_mod >1))
err_m_mod_mod = np.delete(err_m_mod, np.where(err_m_mod >1))

# Removes remaining infinities in logN_mod
m_array_mod_mod = np.delete(m_array_mod_mod, np.where(np.isinf(logN_mod)))
err_N_mod = np.delete(err_N_mod, np.where(np.isinf(logN_mod)))
N_mod = np.delete(N_mod, np.where(np.isinf(logN_mod)))
err_m_mod_mod = np.delete(err_m_mod_mod, np.where(np.isinf(logN_mod)))
logN_mod = np.delete(logN_mod, np.where(np.isinf(logN_mod)))

# Selects the relevant points for plotting
m_array_mod_mod = np.delete(m_array_mod_mod, np.where(logN_mod < 0.4))
err_N_mod = np.delete(err_N_mod, np.where(logN_mod < 0.4))
N_mod = np.delete(N_mod, np.where(logN_mod < 0.4))
err_m_mod_mod = np.delete(err_m_mod_mod, np.where(logN_mod < 0.4))
logN_mod = np.delete(logN_mod, np.where(logN_mod < 0.4))

# The error on a log base 10 function
err_logN = 0.434*(err_N_mod/N_mod)

#%%
pl.figure()
def Fit_line(b, m, c):
    y = m*b + c
    return y

# x_array defines the range of values used for the linear fit
x_array = np.arange(11.9972, 16.3131, 0.001)
p0 = np.array([0.5, -2])                  
p, cov = opt.curve_fit(Fit_line, m_array_mod_mod[155:606], logN_mod[155:606], p0)
plt.plot(x_array, Fit_line(x_array, p[0], p[1]), zorder=10,color = 'red')
pl.errorbar(m_array_mod_mod, logN_mod, yerr = err_logN , color = "forestgreen", fmt='o', mew=1, ms=0.2, capsize=3)
pl.errorbar(m_array_mod_mod,logN_mod, xerr = err_m_mod_mod, color = "orange", fmt='o', mew=1, ms=0.2, capsize=3)
plt.plot(m_array_mod_mod,logN_mod, "x", label = "Variable Aperture")

#plt.xlabel("Apparent magnitude of galaxy", size = "13")
#pl.ylabel("Number of sources detected \n brighter than the apparent magnitude", size = "13")
plt.xlabel("Magnitude [$mag$]", size = "13")
pl.ylabel("log(N(<m))", size = "13")
plt.legend()
plt.grid()
pl.show()
#%%
""" Writing the data to a csv file """
import pandas as pd

indx = np.arange(1, len(x)+1, 1)
positions_array = np.column_stack((indx, x,y, max_rad, aperture_list_final,\
                                   bckg_mean_list, bckg_vals_final, residuals_list, err_residuals))
names = ["id", 'x position', 'y position', 'Aperture radius', 'Total aperture flux', \
         'Local mean background per pixel', 'Total background in aperture',\
         'Corrected Aperture flux', 'Error on corrected flux']
df = pd.DataFrame(positions_array, columns=names)
df.to_csv('Catalogue.csv', index=False)