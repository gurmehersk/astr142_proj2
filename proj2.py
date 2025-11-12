'''
proj2.py
---------

This module crates a 3-image mosaic of the Hubble Ultra Deep Field (HUDF) using optical ACS/WFC images.

I located 3 "drizzle" fits files of different wavelengths for filters  F435W, F606W, and F775W. These images were 
from 2003, but seem to work

https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html

To create the RGB Color composites, I used the astropy website
 
https://docs.astropy.org/en/stable/visualization/rgb.html

"RGB images using the Lupton et al (2004) scheme"

For image scaling,

https://docs.astropy.org/en/stable/visualization/normalization.html

We query using 
the 2017A&A...608A...3B bibcode query from SIMBAD 
'''

# Color mapping documentation
#
# The three filters used in this composite correspond to optical bandpasses:
#   - F435W (~0.43 µm): Blue channel (B-band)
#   - F606W (~0.59 µm): Green channel (V-band)
#   - F775W (~0.77 µm): Red channel (I-band)
#
# These assignments follow the natural wavelength order (short → long)
# and approximate the colors our eyes would perceive if we could see
# the field directly.

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
# New import, taken from astropy documentation as mentioned in module docstring
from astropy.visualization import make_lupton_rgb
import logging
from astropy.visualization import MinMaxInterval, PercentileInterval, AsinhStretch, ImageNormalize
from astropy.visualization import ZScaleInterval
from astropy.visualization import quantity_support
import IPython


from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astroquery.simbad import Simbad


# configuring the logger, might put this in the "main" if i end up creating one
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

quantity_support()

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
    
    #### To keep the data alignment consistent, I'm replacing the NaN values to 0 values
    data = np.nan_to_num(data, nan=0.0)
    
    # get the wcs information
    wcs = WCS(hdr)
    
    return data, wcs


def scaling(data):
    norm = ImageNormalize(data, interval=ZScaleInterval(), stretch=AsinhStretch())
    return norm(data)
    
def rgb_img(green,blue,red,wcs):
    '''
    Used to create actual 3-color composite
    
    using Asinhstretch for scaling
    I'm not actually using this stretch anymore, im using the Lupton 2004 
    package 
    
    https://docs.astropy.org/en/stable/api/astropy.visualization.AsinhStretch.html
    
    Apparently good for rgb color images in astronomy 
    '''
    try:
        if not (red.shape == green.shape == blue.shape):
            logger.error("Input files have misaligned shapes, error..")
            print("Misaligned shapes")
            print(red.shape)
            print(green.shape)
            print(blue.shape)
            raise ValueError("All input images must have the same dimensions.")
    except:
        print("Exiting code...")
        sys.exit()
        
    red_scaled = scaling(red)
    green_scaled = scaling(green)
    blue_scaled = scaling(blue)
    
    
    rgb_image = np.dstack((red_scaled, green_scaled, blue_scaled))
    
    ## Trying this Lupton thing because the stacking is oversaturating image
    rgb_image = make_lupton_rgb(red_data, green_data, blue_data, stretch=0.01, Q=10)
    
    logger.info("Plotting RGB composite")
    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot(projection=wcs) if wcs else plt.subplot()
    ax.imshow(rgb_image, origin='lower')
    ax.set_title(f"RGB Composite")
    ax.set_xlabel("RA (J2000)")
    ax.set_ylabel("Dec (J2000)")
    plt.tight_layout()
    plt.show()

    logger.info("RGB mosaic created successfully.")
    return rgb_image
    
if __name__ == "__main__":

    logging.info("Creating the 3-color mosaic for the Hubble Ultra Deep Field")
    
    # 435 is for Blue light
    #f_435 = "/Users/gurmeher/Downloads/proj2_astr142/435fw/HST/hst_10086_7d_acs_wfc_f435w_j8wc7d/hst_10086_7d_acs_wfc_f435w_j8wc7d_drc.fits"
    f_435 = "/Users/gurmeher/Downloads/h_udf_wfc_b_drz_img.fits"
    #f_435 = "/Users/gurmeher/Downloads/proj2_astr142/hst_10086_7d_acs_wfc_f435w_j8wc7d_drc.fits"
    # 606 is for Green Light
    #f_606 = "/Users/gurmeher/Downloads/proj2_astr142/606fw/HST/hst_10086_87_acs_wfc_f606w_j8wc87/hst_10086_87_acs_wfc_f606w_j8wc87_drc.fits"
    f_606 = "/Users/gurmeher/Downloads/h_udf_wfc_v_drz_img.fits"
    
    # 775 is for Red Light
    #f_775 = "/Users/gurmeher/Downloads/proj2_astr142/775fw/HST/hst_10086_7e_acs_wfc_f775w_j8wc7e/hst_10086_7e_acs_wfc_f775w_j8wc7e_drc.fits"
    #f_775 = "/Users/gurmeher/Downloads/proj2_astr142/hst_10086_7e_acs_wfc_f775w_j8wc7e_drc.fits"
    f_775 = "/Users/gurmeher/Downloads/h_udf_wfc_i_drz_img.fits"
    
    #f_850 = "/Users/gurmeher/Downloads/proj2_astr142/hst_10086_7c_acs_wfc_f850lp_j8wc7c_drc.fits"
    
    blue_data, wcs = fits_image_loader(f_435)
    green_data, _  = fits_image_loader(f_606)  # WCS same
    red_data, _    = fits_image_loader(f_775)
    
    rgb_image = rgb_img(green_data, blue_data, red_data, wcs)
    
    from scipy.ndimage import rotate
    rotated_img = rotate(rgb_image, angle=45, reshape=False)
    
     # Optionally save the image as PNG
    plt.imsave("HST_RGB_composite.png", rgb_image, origin='lower')
    plt.savefig("HST_RGB_composite.pdf", bbox_inches='tight')
    plt.show()
    logging.info("Saved RGB composite as 'HST_RGB_composite.png' and pdf version")
    
    



