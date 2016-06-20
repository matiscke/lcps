# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:14:43 2016

@author: mschleck
"""

import os
from astropy.table import Table, vstack
from lcps_io import open_fits
from astropy import log
import slidingWindow

def lcps_output(logtable, logfile):
    """ Write table with transit candidates to file."""
#    if os.path.isfile(logfile):
#        # append to existing logfile 
#        oldlog = ascii.read(logfile, format='csv')
#        logtable = vstack([oldlog, logtable])
    logtable.write(logfile, format='csv')
        

def batchjob(path, logfile='./dips.log', winSize=10, stepSize=1,\
        Nneighb=1, minDur=2, maxDur=5, detectionThresh=0.995):
    """ Check all light curve files in a folder for transit signatures.
    
    batchjob forwards all FITS files in the `path` to the dip search of the 
    `slidingWindow` module. Any detected dips are stored in an Astropy table
    `candidates` together with the EPIC number of the target.
    
    Parameters
    ----------
    path : str
        folder which is scanned for fits files
    logfile : str
        output file for transit candidates
    winSize : int
        Size of a sliding window
    stepSize : int
        steps per slide (Default = 1, i.e. slide one data point per iteration)
    Nneighb : int
        Number of neighboring windows per side to be considered for the local 
        median (At the boundaries of the time series, the considered data
        extends to the beginning or end of the array, respectively)    
    minDur : int
        minimum dip duration in # of data points
    maxDur : int
        maximum dip duration in # of data points
    detectionThresh : float
        fraction of flux below which a dip is registered
    
    Returns
    -------    
    candidates : Astropy table
        table with detected dips
        
    Example
    -------
    >>> path = './tests/'
    >>> candidates = batchjob(path)
    INFO: Scanning target 1/1: EPIC 205919993 [__main__]
    INFO: 12 transit candidates found in 1 light curves. [__main__]
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
        candidates = vstack([candidates, dips], join_type='outer')
        
        ####### DEBUGGING
        if dips:
            print 'transit-like event detected'
        
        #####################################
        
        
        
        
        
        
        # Every 50th iteration, write intermediate results to file
        if i % 50 == 0:
            lcps_output(candidates, logfile + '.part')
        
    # write transit candidates to file
    lcps_output(candidates, logfile)
    try:
        os.remove(logfile + '.part')
    except OSError:
        pass
    
    log.info('{} transit candidates found in {} light curves.'.format(\
        len(candidates), len(set(candidates['EPIC']))))
    return candidates
    
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
    # parse parameters
    import argparse
    parser = argparse.ArgumentParser(\
        description='pre-select light curves with possible transit signatures')
    parser.add_argument('path',\
        help='path containing light curve (FITS) files', type=str)
    parser.add_argument('--logfile', default='./dips.log',\
        help='name of log file that will contain transit candidates', type=str)  
    parser.add_argument('--winSize', default=50,\
        help='Size of a sliding window', type=int)
    parser.add_argument('--stepSize', default=10,\
        help='steps per slide (Default = 1, i.e. slide one data point per iteration)', type=int)
    parser.add_argument('--Nneighb', default=1,\
        help='Number of neighboring windows to be considered for the local median', type=int)
    parser.add_argument('--minDur', default=2,\
        help='minimum dip duration in # of data points', type=int)
    parser.add_argument('--maxDur', default=49,\
        help='maximum dip duration in # of data points', type=int)
    parser.add_argument('--detectionThresh', default=0.98,\
        help='fraction of flux below which a dip is registered', type=float)
    args = parser.parse_args()
    
    batchjob(args.path, args.logfile, args.winSize, args.stepSize,\
        args.Nneighb, args.minDur, args.maxDur, args.detectionThresh)

    
##### DEBUGGING
#path = '/run/media/mschleck/scratch2/KeplerData/C7/'
#logfile = '/run/media/mschleck/scratch2/KeplerData/C7/lcps_debugging01.log'
#winSize = 500
#stepSize = 10
#Nneighb=1
#minDur = 10
#maxDur = 499
#detectionThresh = 0.96
#candidates = batchjob(path, logfile,winSize,stepSize,Nneighb,minDur,maxDur,detectionThresh)
#print candidates


