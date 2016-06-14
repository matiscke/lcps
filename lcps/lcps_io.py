# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:22:51 2016

@author: mschleck
"""
import numpy as np
from astropy.table import Table
from astropy.io import fits
import warnings
from astropy.utils.exceptions import AstropyUserWarning


def open_fits(filename):
    """ Open a light curve file in the usual Kepler FITS format and extract
    the PDCSAP light curve.
    
    Parameters
    ----------
    filename : str
        file name of the FITS file containing the light curve data
        
    Returns
    -------
    EPICno : int
        EPIC number ('KEPLERID' in the hdu header)
    photometry : Astropy table
        Columns are time in BJD - 2454833, PDCSAP flux, PDCSAP flux error
        
    Example
    -------
    >>> filename = 'tests/ktwo205919993-c03_llc.fits'
    >>> EPICno, photometry = open_fits(filename)
    """
    try:
        hdulist = fits.open(filename)
    except IOError:
        warnings.warn("Could not open FITS file.", AstropyUserWarning)       
        return None  
    EPICno = hdulist[1].header['KEPLERID']
    tbdata = hdulist[1].data
    hdulist.close()

    # extract light curve data from hdu
    time = tbdata['TIME']
    flux = tbdata['PDCSAP_FLUX']
    flux_err = tbdata['PDCSAP_FLUX_ERR']
    photometry = Table([time, flux, flux_err], names = ('TIME', 'FLUX','FLUX_ERR'))
    
    # remove nans
    photometry = photometry[~np.isnan(photometry['FLUX'])]
    return EPICno, photometry
  
  
if __name__ == "__main__":
    import doctest
    doctest.testmod()