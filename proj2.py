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
https://ui.adsabs.harvard.edu/abs/2015AJ....150...31R/abstract
this is for 

photometric  redshifts: https://vizier.cds.unistra.fr/viz-bin/VizieR-4
spectroscopic redshifts: https://vizier.cds.unistra.fr/viz-bin/VizieR-4
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
from astroquery.vizier import Vizier


# configuring the logger, might put this in the "main" if i end up creating one
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

quantity_support()

def query_spectroscopy_catalog():
    '''
    Query Vizier for MUSE spectrscopic redshift catalog
    
    Returns
    -------
    c : Table
        Spectroscopic catalog with RA, DEC, z_spec
    '''
    logger.info("Querying Vizier for (Inami+ 2017) MUSE catalog")
    
    v = Vizier(columns=["*"], row_limit=-1)
    
    try:
        catalog_list = v.get_catalogs("J/A+A/608/A2")
        
        if len(catalog_list) == 0:
            logger.error("No catalog found!")
            return None
        
        # Get the combined table
        c = catalog_list[0]
        
        logger.info(f"Got {len(c)} objects with spectroscopic redshifts")
        logger.info(f"Columns: {c.colnames}")
        
        # I noticed that the rafelski and MUSE ones are different only in the decimal points since there is rounding rafelski, so im gonna round them to be the same.
        if 'RAJ2000' in c.colnames:
            c['RAJ2000'] = np.round(c['RAJ2000'], 6)
        if 'DEJ2000' in c.colnames:
            c['DEJ2000'] = np.round(c['DEJ2000'], 6)
            
        return c
        
    except Exception as e:
        logger.error(f"VizieR query failed: {e}")
        return None



def query_photometry_catalog():
    """
    Query Vizier for Rafelski et al. 2015 photometric redshift catalog.
    
    Returns
    -------
    c : Table
        Photometric catalog with RA, DEC, z_phot
    """
    
    ### We will be putting a magnitude limit, otherwise getting too crowded
    logger.info("Querying VizieR for catalog (Rafelski+ 2015)")
    
    v = Vizier(columns=["*"], row_limit=-1)
    
    try:
        # Query catalog by the name written
        catalog_lst = v.get_catalogs("J/AJ/150/31")
        
        if len(catalog_lst) == 0:
            logger.error("No catalog found!")
            return None
        
        # Get the main table thats labelled table5 on Vizier
        c = catalog_lst[0]
        
        logger.info(f"Got {len(c)} objects with photometric redshifts")
        logger.info(f"Columns: {c.colnames}")
        
        return c
        
    except Exception as e:
        logger.error(f"VizieR query failed: {e}")
        return None

def cross_match_catalogs(photo_cat, spec_cat, ra_col_photo, dec_col_photo, ra_col_spec, dec_col_spec, match_radius):
    '''
    We are going to check the RA and DEC of objects with extreme precision 
    
    Parameters
    -----------
    photo_cat : csv file
        Contains the rafelski catalog csv file, with ra, dec and everything we will need to do the ra dec verification
    
    spec_cat : csv file
        Contains the muse spectroscopic data csv file, with ra, dec and everything we will need to the verification, and with the ra and dec rounded to the 6th decimal place to ensure same precision in both catalogs
        
    ra_col_photo : string
        Name of the column in the photometric catalog that contains the RA
    
    dec_col_photo : string
        Name of the column in the photometric catalog that contains the DEC
        
    ra_col_spec : string
        Name of the column in the spectroscopic catalog that contains teh RA
        
    dec_col_spec : string
        Name of the column in the spectroscopic catalog that contains teh DEC
                
    match_radius : int (arcsecond units)
        The tolerance for matching 
        
    Returns
    -----------
    photo_indices : ndarray
        Indices of matched sources in photometric catalog
    spec_indices : ndarray
        Indices of matched sources in spectroscopic catalog
    '''
    
    # Let's create SkyCoord objects for both catalogs
    photo_coords = SkyCoord(ra=photo_cat[ra_col_photo].value*u.deg, dec=photo_cat[dec_col_photo].value*u.deg)
    
    spec_coords = SkyCoord(ra=spec_cat[ra_col_spec].value*u.deg, dec=spec_cat[dec_col_spec].value*u.deg)
    
    
    # Okay, this is a new thing i had to search up for the cross-verification cuz i had no idea how to do it
    # The match_to_catalog_sky returns 3 values, the index of spec_coords that is closest to the photo_coords, its separation, and 3D separation, but we dont need 3D separation
    idx_spec, sep, _ = photo_coords.match_to_catalog_sky(spec_coords)
    
    # Only keep matches within the specified radius
    matched_mask = sep < match_radius
    photo_indices = np.where(matched_mask)[0]
    spec_indices = idx_spec[matched_mask]
    
    ### The above method might be a little hard to grasp, but essentially we see if the sep < match_radius, and if that's true, we make the photo_indices wherever the mask is true, and spec index becomes the indexes of true. So this way, the order and alignment is preserved for the photometric and sepectroscopic galaxies
    match_fraction = 100.* len(photo_indices) / len(photo_cat)
    logger.info(f"{len(photo_indices)} matched sources ({match_fraction:.1f}% of photo catalog)")
    
    
    ## Without this line, life is becoming too complicated
    photo_cat['has_spec'] = False
    photo_cat['has_spec'][photo_indices] = True
    
    
    logger.info(f"Added 'has_spec' column to photometric catalog")
    return photo_indices, spec_indices
    
def overlayer(green, blue, red, wcs, photo_cat, matched_indices,
                          ra_col='RAJ2000', dec_col='DEJ2000'):
    '''
    Creates RGB composite with catalog overlays [the other functions to just plot the thing might be irrelevant now]
    
    Uses the 'has_spec' column added by cross_match_catalogs to color sources.
    
    Parameters
    ----------
    green, blue, red : ndarray
        Image data for each channel
    wcs : WCS
        World coordinate system
    photo_cat : Table
        Photometric catalog (with 'has_spec' column added)
    matched_indices : tuple
        (photo_idx, spec_idx) of matched sources
    ra_col, dec_col : str
        Column names for coordinates
        
    Returns
    -------
    fig : Figure
        Matplotlib figure
    rgb_image : ndarray
        RGB composite image
    '''
    if not (red.shape == green.shape == blue.shape):
        raise ValueError("All input images must have the same dimensions.")
    
    rgb_image = make_lupton_rgb(red, green, blue, stretch=0.01, Q=10)
    
    fig = plt.figure(figsize=(16, 14))
    ax = plt.subplot(projection=wcs)
    ax.imshow(rgb_image, origin='lower')
    
    ra_all = photo_cat[ra_col]
    dec_all = photo_cat[dec_col]
    # this is because we need to make sure the circle is plotted at the right place, and not in the wrong area over the wrong galaxies
    pixel_coords = wcs.world_to_pixel_values(ra_all, dec_all)
    
    # Plot photo-z only, i.e. places where the (has_spec == False) in cyan
    mask_photo_only = ~photo_cat['has_spec'] # using logical not for quicker and concise code D:
    
    n_photo_only = np.sum(mask_photo_only)
    logger.info(f"Plotting {n_photo_only} photo-z only sources")
    
    ax.scatter(pixel_coords[0][mask_photo_only], pixel_coords[1][mask_photo_only],
              s=25, facecolors='none', edgecolors='cyan',
              linewidths=0.6, alpha=0.4,
              label=f'Photo-z only (N={n_photo_only})')
    
    # Plot spec-z available (has_spec == True) in red
    mask_has_spec = photo_cat['has_spec']
    n_spec = np.sum(mask_has_spec)
    logger.info(f"Plotting {n_spec} sources with spec-z")
    
    ax.scatter(pixel_coords[0][mask_has_spec], pixel_coords[1][mask_has_spec],
              s=60, facecolors='none', edgecolors='red',
              linewidths=1.2, alpha=0.8,
              label=f'Spec-z available (N={n_spec})')
    
    ax.set_title("Hubble Ultra Deep Field: Photometric vs Spectroscopic Redshifts\n" +
                "RGB: F775W (Red) / F606W (Green) / F435W (Blue)",
                fontsize=15, fontweight='bold', pad=20)
    ax.set_xlabel("Right Ascension (J2000)", fontsize=13)
    ax.set_ylabel("Declination (J2000)", fontsize=13)
    ax.legend(loc='upper right', fontsize=11, framealpha=0.95)
    ax.grid(color='white', ls='dotted', alpha=0.25)
    
    plt.tight_layout()
    logger.info("Catalog overlay complete")
    
    return fig, rgb_image

                
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

    
def rgb_img(green,blue,red,wcs):
    '''
    Used to create actual 3-color composite
    
    using Asinhstretch for scaling
    I'm not actually using this stretch anymore, im using the Lupton 2004 
    package, but i was using it before, so i've kept it in here
    
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
    
    ## Trying this Lupton thing because the stacking is oversaturating image
    rgb_image = make_lupton_rgb(red, green, blue, stretch=0.01, Q=10)
    
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
    return fig, rgb_image
    
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
    
    fig, rgb_image = rgb_img(green_data, blue_data, red_data, wcs)
    
    '''from scipy.ndimage import rotate
    rotated_img = rotate(rgb_image, angle=45, reshape=False)'''
    
     # Optionally save the image as PNG
    fig.savefig("HST_RGB_composite.png", dpi=300, bbox_inches='tight')
    fig.savefig("HST_RGB_composite.pdf", bbox_inches='tight')
    plt.show()
    logging.info("Saved RGB composite as 'HST_RGB_composite.png' and pdf version")
    
    logger.info("\n[STEP 2] Querying catalogs from VizieR...")
    
    photo = query_photometry_catalog()
    spec = query_spectroscopy_catalog()
    
    if photo is None or spec is None:
        logger.error("failed to get catalogs, try again")
        sys.exit(1)
    
    
    photo.write('rafelski_photometric_catalog.csv', format='csv', overwrite=True)
    spec.write('muse_spectroscopic_catalog.csv', format='csv', overwrite=True)
    
    
    matched_indices = cross_match_catalogs(photo, spec, ra_col_photo='RAJ2000', dec_col_photo='DEJ2000',ra_col_spec='RAJ2000', dec_col_spec='DEJ2000', match_radius=0.04*u.arcsec) # I have faith in the match_radius being this low by eye seeing how 2-3 sources were pretty close in RA and DEC on this thing.

    
    fig_cats, _ = overlayer(green_data, blue_data, red_data, wcs, photo, matched_indices)
    
    fig_cats.savefig("hudf_with_catalogs.png", dpi=300, bbox_inches='tight')
    fig_cats.savefig("hudf_with_catalogs.pdf", bbox_inches='tight')
    logger.info("Saved overlayed mosaic")
    plt.show()
