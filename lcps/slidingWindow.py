# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 18:10:47 2016

@author: mschleck
"""

import numpy as np
from astropy.table import Table, vstack

def findDip(fluxWindow, minDur=1, maxDur=5, localMedian=1.00, detectionThresh=0.995):
    """ Search for negative deviations, i.e. dips, in an array.
    
    findDip counts the number of entries in fluxWindow that fall short of a 
    threshold in relative flux specified in detectionThresh. If this number
    lies between minDur and maxDur, a dip is detected and findDip returns the
    time of egress and the minimum flux relative to localMedian.
    
    Parameters
    ----------
    fluxWindow : Astropy Table
        A table containing columns "t", "flux"
    minDur : int
        minimum dip duration in # of data points
    maxDur : int
        maximum dip duration in # of data points
    localMedian : float
        local median flux
    detectionThresh : float
        fraction of flux, below which a deviation is registered
    
    Returns
    -------
#    dip : bool
#        True, if dip of length between [minDur, maxDur] detected
    
    t_egress : float
        time at end of detected dip 
    minFlux : float
        Minimum flux relative to localMedian
        
    Example
    -------
    >>> fluxWindow = Table([[0.,1.,2.,3.,4.,5.], [1.00,1.01,0.99,0.80,0.75,0.95]],\
        names=["t","flux"], dtype=[float, float])
    >>> findDip(fluxWindow, detectionThresh=0.90)
    (5.0, 0.75)
    """
    
    ###
    # put parameters check out of this function to speed things up    
    #
    # Check if parameters are consistent
    if minDur > maxDur:
        raise ValueError('min dip duration larger than max dip duration')
    if len(fluxWindow) <= maxDur:
        raise ValueError('max dip duration greater than or equal Window size')
     
    fluxThresh =  detectionThresh*localMedian
    if len(fluxWindow[fluxWindow["flux"] < fluxThresh]) >= minDur:
        # There are low fluxes, check for coherence
            NloFlux = 0
            for i, datum in enumerate(fluxWindow["flux"]):
                if datum < fluxThresh:
                    NloFlux += 1
                else:                    
                    # End of dip, check if length falls between limits
                    if minDur <= NloFlux <= maxDur:
                        # return time of egress and min. flux rel. to median
                        return fluxWindow["t"][i], \
                            np.min(fluxWindow["flux"][:i])/localMedian
                    else:
                        # look for additional dips in the Window
                        NloFlux = 0
                        continue    
    # No dips found
    return None, None
                    

def localMedian(photometry, iWinStart, winSize, Nneighb=2):
    """ Find the local median of fluxes, ignoring the current window.
    
    localMedian computes the median flux in the neighboring windows. The fluxes
    in the current window is ignored.
    
    Parameters
    ----------
    photometry : Astropy Table
        A table with the whole photometric data containing columns "t", "flux"
    iWinStart : int
        Index of the first datum in the current window
    winSize : int
        Size of a window
    Nneighb : int
        Number of neighboring windows to be considered for the averaging (At
        the boundaries of the time series, the considered data extends to the
        beginning or end of the array, respectively)
    
    Returns
    -------
    localMedian : float
        Median of the flux in the windows neighboring the current one
        
    Example
    -------
    >>> photometry = Table([[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.],\
        [1.00,1.01,0.99,0.80,0.75,0.95,0.99,0.99,1.00,0.80,1.01]],\
        names=["t","flux"], dtype=[float, float])
    >>> localMedian(photometry, 4, 4, Nneighb=1)
    1.0
    """
    
    # At the boundaries, expand neighborhood towards center
    if (iWinStart < winSize) or (iWinStart > len(photometry) - winSize):
        Nneighb *= 2

    # construct neighborhood without current window        
    iMin = max(0, iWinStart - Nneighb*winSize)
    iMax = min(len(photometry), iWinStart + (1 + Nneighb)*winSize)
    neighborhood = vstack([photometry[iMin:iWinStart],\
        photometry[iWinStart + winSize:iMax + 1]])
        
    return np.median(neighborhood['flux'])
    
   
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
##############
#    #TEST
#fluxWindow = Table([[0.,1.,2.,3.,4.,5.], [1.00,1.01,0.99,0.80,0.75,0.95]],\
#        names=["t","flux"], dtype=[float, float])
#t_egress, minFlux = findDip(fluxWindow, detectionThresh=0.90) 
#print t_egress, minFlux
    
#photometry = Table([[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.],\
#        [1.00,1.01,0.99,0.80,0.75,0.95,0.99,0.99,1.00,0.80,1.01]],\
#        names=["t","flux"], dtype=[float, float])
#median = localMedian(photometry, 4, 4, Nneighb=1)