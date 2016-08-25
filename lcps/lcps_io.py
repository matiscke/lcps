# -*- coding: utf-8 -*-
""" Auxiliary functions for light curve file handling. 

For now, contains only a function "open_fits" to extract the Kepler PDCSAP 
light curve from a FITS file.
"""

import numpy as np
from astropy.table import Table
from astropy.io import fits, ascii
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

def open_csv(filename):
    """ Open a light curve file in csv format and extract from it flux time
    series.
    
    Parameters
    ----------
    filename : str
        file name of the ascii file containing the photometry
    
    Returns
    -------
    filename : str
        The filename serves as a unique identifier for the object
    photometry : Astropy table
        Columns are named after the file header and contain time, flux
    """
    photometry = ascii.read(filename, format='csv')
    return filename, photometry
 
def open_k2sff(filename):
    """ Extract a light curve from a 'K2SFF' ascii file. The default aperture
    light curve data of this product are not strictly 'comma-separated' and 
    lead to crashes when opened by standard Astropy ascii I/O functions. 
    
    Parameters
    ----------
    filename : str
        file name of the ascii file containing the photometry
    
    Returns
    -------
    filename : str
        The filename serves as a unique identifier for the object
    photometry : Astropy table
        Columns contain time, flux 
        
    Example
    -------
    >>> filename = 'tests/220132548'
    >>> filename, photometry = open_k2sff(filename)
    """
    with open(filename, 'r') as infile:
        lines = infile.readlines() 
        phot = np.zeros([len(lines) - 1, 2])
        for i, line in enumerate(lines[1:]):
            # strip trailing comma and '\n' and save to a table
            line = line.rstrip(',\n')
            line = line.split(',')
            phot[i][:] = line
        photometry = Table(phot, names = ('TIME', 'FLUX'))    
    
    # remove nans
    photometry = photometry[~np.isnan(photometry['FLUX'])]   
    return filename.split('/')[-1], photometry
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    