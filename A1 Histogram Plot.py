# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:44:23 2020

@author: rr3
"""
#import modules
from astropy.io import fits
import scipy as sp
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#LOADING DATA
hdulist = fits.open('A1_mosaic.fits')

headers = hdulist[0].header
dataset = hdulist[0].data
d = np.flipud(dataset)
ds = np.flipud(dataset)

#Get rid of blooming
d[480:540, 540:580] = 0
d[1150:1400, 750:800] = 0
d[1760:1920, 950:1000] = 0
d[800:900, 2120:2160] = 0
d[2450:2480, 1160:1240] = 0
d[3160:3210, 2070:2100] = 0
d[1660:4611, 1420:1500] = 0
d[1620:1660, 1423:1450] = 0
d[0:1620, 1400:1470] = 0
d[2266:2333, 2120:2145] = 0
d[1162:1231, 2452:2475] = 0
d[2240:2390, 890:920] = 0
d[4600:4611, 969:1720] = 0
d[4150:4600, 1250:1420] = 0
d[4460:4520, 1500:1550] = 0
d[4550:4600, 1630:1650] = 0
d[4280:4305, 1000:1250] = 0
d[4170:4190, 1499:1650] = 0
d[4250:4300 ,1499:1705] = 0
d[1300:1500, 1350:1525] = 0
print(np.amax(d))
#%%
#isolate background
data = []
#remove edges from image
for y in range(550,4000):
    for x in range(300, 2300):
        if d[y][x]<4000 and d[y][x]>2000:
            data.append(d[y][x])
        else:
            continue
        
#%%
#extract data from the histogram
freq, count, q = plt.hist(data, bins = 1000)

#find the centre of bins
bins_centre = []
for i in range(1,len(count)):
    x = (count[i]+count[i-1])/2
    bins_centre.append(x)

#define fitting function
def gaus(x,A,mu,sig):
    return A*np.exp(-((x-mu)**2)/(2*(sig**2)))

#do the curve-fit
#guess = [3000,100,800000]
popt,pcov = curve_fit(gaus,bins_centre,freq,p0=[1000000,3500,100])

#plot the data
background_range = sp.arange(min(bins_centre),max(bins_centre),0.5)
plt.plot(background_range,gaus(background_range,*popt))

print(popt[0],popt[1],popt[2])
plt.show()