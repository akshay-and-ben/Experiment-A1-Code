# Experiment-A1-Code - Our code for Experiment A1 of the 3rd Lab Cycle (Department of Physics, Imperial College London)
# Authors: Akshay Robert (AR6918) and Ben Amroota (BRA18)

Below is a list detailing the contents of the attached files in the 'Experiment-A1-Code' folder and their function in relation to the overarching goal of producing the galactic catalogue:

    1. A1 Histogram Plot - plots a histogram of pixel values in the image and calculates the mean background value 
    2. Manually mask blooming stars, and remove edges of image due to
    increased noise
    3. Identifies galaxies present within the image and notes the positions
    4. Uses the positions in aperture photometry using the Astropy module to
    determine the total flux through a variable aperture, the local background
    for each galaxy, and return a list of corrected aperture fluxes.
    5. Plots the data in the form log(N(<m)) vs m, and fits a linear fit to it
    using the scipy.optimize.curvefit module
    6. Writes the data to a csv file for external analysis.

