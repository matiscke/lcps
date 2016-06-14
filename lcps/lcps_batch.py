# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:14:43 2016

@author: mschleck
"""

import numpy as np
import os
from astropy.table import Table, vstack
from lcps_io import open_fits
import slidingWindow


def batchjob(path, winSize=10, stepSize=1, Nneighb=2, minDur=2, maxDur=5,
        detectionThresh=0.995):
    """ Check all light curve files in a folder for transit signatures.
    
    batchjob forwards all FITS files in the `path` to the dip search of the 
    `slidingWindow` module. Any detected dips are stored in an Astropy table
    `candidates` together with the EPIC number of the target.
    
    Parameters
    ----------
    path : str
        folder which is scanned for fits files
    winSize : int
        Size of a sliding window
    stepSize : int
        steps per slide (Default = 1, i.e. slide one data point per iteration).
    Nneighb : int
        Number of neighboring windows to be considered for the local median (At
        the boundaries of the time series, the considered data extends to the
        beginning or end of the array, respectively)    
    minDur : int
        minimum dip duration in # of data points
    maxDur : int
        maximum dip duration in # of data points
    detectionThresh : float
        fraction of flux, below which a deviation is registered
    
    Returns
    -------    
    candidates : Astropy table
        table with detected dips
        
    Example
    -------
    >>> path = './tests/'
    >>> candidates = batchjob(path)
    """
    filelist = [file for file in os.listdir(path) if file.endswith('fits')]
    candidates = Table(names=('EPIC','t_egress','minFlux'), dtype=['S8',float,float])
    for file in filelist:
        # extract photometry from fits file
        EPICno, photometry = open_fits(path + file)
        
        # Search for transit signatures via sliding window algorithm
        dips = slidingWindow.dipsearch(EPICno, photometry, winSize, stepSize,\
            Nneighb, minDur, maxDur, detectionThresh)
#        candidates.add_row([EPICno, dips['t_egress'], dips['minFlux']])
        candidates = vstack([candidates, dips], join_type='outer')
    return candidates
    
    
#if __name__ == "__main__":
#    import doctest
#    doctest.testmod()
    
    
### DEBUGGING
path = './tests/'
candidates = batchjob(path)


