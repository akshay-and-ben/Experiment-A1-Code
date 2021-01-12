# Experiment-A1-Code - Our code for Experiment A1 of the 3rd Lab Cycle (Department of Physics, Imperial College London)
# Authors: Akshay Robert (AR6918) and Ben Amroota (BRA18)

Below is a list detailing the contents of the attached files in the 'Experiment-A1-Code' folder and their function in relation to the overarching goal of producing the galactic catalogue:

    1. A1 Histogram Plot.py: Plots a histogram of pixel values in the image and calculates the mean background value. 
    2. A1 Astronomical Imaging Code:
        - Manually mask blooming stars, and remove edges of image due to increased noise.
        - Identifies galaxies present within the image and notes the positions.
        - Uses the positions to perform aperture photometry from tools in the Asttropy module to determine the total flux 
          through a variable aperture, the local background for each galaxy, and return a list of corrected aperture fluxes.
        - Plots a number count of log(N(<m)) vs m, and fits a linear function to it using the scipy.optimize.curvefit module.
        - Writes the pertaining data to a CSV file for further analysis.
    3. Testing Code - Single Pixel     : Assesses if the A1 Astronomical Imaging Code can identify single pixels (caused by 
                                         single photons incident on the CCD)
    4. Testing Code - Single Star      : Assesses if the A1 Astronomical Imaging Code can identify sources in the image
    5. Testing Code - Two Stars        : Assesses if the A1 Astronomical Imaging Code can identify two independent sources in the image
    6. Testing Code - Two Stars & Pixel: Assesses if the A1 Astronomical Imaging Code can identify two independent sources as in (5)
                                         and differentiate these from single pixel events as in (3)
    7. Catalogue.csv -  a catalogue of galactic sources including data such as coordinate positons, total flux emitted from source, local
                        background values and the size of the aperture used to perform photometry

The FITS data file analysed (A1_mosaic.fits) is not enclosed in this folder due to file size considerations

