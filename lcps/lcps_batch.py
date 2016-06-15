# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:14:43 2016

@author: mschleck
"""

import numpy as np
import os
from astropy.table import Table, vstack
from astropy.io import ascii
from lcps_io import open_fits
from astropy import log
import slidingWindow

def lcps_output(logtable, logfilename):
    """ Write table with transit candidates to file."""
    if os.path.isfile(logfilename):
        # append to existing logfile 
        oldlog = ascii.read(logfilename, format='csv')
        logtable = vstack([oldlog, logtable])
    logtable.write(logfilename, format='csv')
        

def batchjob(path, logfilename='./dips.log', winSize=10, stepSize=1,\
        Nneighb=2, minDur=2, maxDur=5, detectionThresh=0.995):
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
    INFO: Scanning target 1/1: EPIC 205919993 [__main__]
    """
    filelist = sorted([file for file in os.listdir(path) if file.endswith('fits')])
    candidates = Table(names=('EPIC','t_egress','minFlux'),\
        dtype=['i8',float,float])
    for i, file in enumerate(filelist):
        # extract photometry from fits file
        EPICno, photometry = open_fits(path + file)
        log.info('Scanning target {}/{}: EPIC {}'.format(i + 1,\
            len(filelist),EPICno))
        
        # Search for transit signatures via sliding window algorithm
        dips = slidingWindow.dipsearch(EPICno, photometry, winSize, stepSize,\
            Nneighb, minDur, maxDur, detectionThresh)
#        candidates.add_row([EPICno, dips['t_egress'], dips['minFlux']])
        candidates = vstack([candidates, dips], join_type='outer')
        
    # write transit candidates to file
    lcps_output(candidates, logfilename)
    return candidates
    
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
    
#### DEBUGGING
#path = '/run/media/mschleck/scratch2/KeplerData/DADS_20160517/'
#logfile = '/run/media/mschleck/scratch2/KeplerData/DADS_20160517/lcps.log'
#candidates = batchjob(path, logfile)
#print candidates

