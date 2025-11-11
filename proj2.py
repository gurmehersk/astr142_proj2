'''
proj2.py
---------

This module crates a 3-image mosaic of the Hubble Ultra Deep Field (HUDF) using optical ACS/WFC images.

I located 3 "drizzle" fits files of different wavelengths for filters  F435W, F606W, and F775W. These images were 
from 2003, but seem to work

https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html


'''

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import make_lupton_rgb
import logging

# configuring the logger, might put this in the "main" if i end up creating one
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")

def fits_image_loader(filename):
'''
    Parameters
    -----------
    
    filename : string
        Filepath of the fits file that needs to be loaded
        
    Returns
    ----------
    data: integer array
        It contains the data from the fits file
        
    wcs : object
        Contains the wcs information extracted from the fits file 
'''
    if not os.path.exists(filename):
        raise FileNotFoundError("The fits file input was not found, please retry")
    
    data, hdr = fits.getdata(filename, header = True)
    
    if data is None:
        raise ValueError("There is no data inside this fits file!")
    
    logging.info("Successfully loaded the data and information for FITS File")
    
    # get the wcs information
    wcs = WCS(hdr)
    
    return data, wcs


if __name__ == "__main__":
    logging.info("Creating the 3-color mosaic for the Hubble Ultra Deep Field")
    
    # 435 is for Blue light
    f_435 = "/Users/gurmeher/Downloads/proj2_astr142/435fw/HST/hst_10086_7d_acs_wfc_f435w_j8wc7d/hst_10086_7d_acs_wfc_f435w_j8wc7d_drc.fits"
    
    # 606 is for Green Light
    f_606 = "/Users/gurmeher/Downloads/proj2_astr142/606fw/HST/hst_10086_87_acs_wfc_f606w_j8wc87/hst_10086_87_acs_wfc_f606w_j8wc87_drc.fits"
    
    # 775 is for Red Light
    f_775 = "/Users/gurmeher/Downloads/proj2_astr142/775fw/HST/hst_10086_7e_acs_wfc_f775w_j8wc7e/hst_10086_7e_acs_wfc_f775w_j8wc7e_drc.fits"
    
    



