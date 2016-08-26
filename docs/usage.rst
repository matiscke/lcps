Command Line Usage
==================

You can use lcps as a command line tool, too. Just run lcps_batch.py from the shell and add your desired parameters. You can show a help screen with ::

   
   $ python lcps_batch.py --help
  
The minimum information lcps needs is the path containing the light curves to be processed. The command ::
   
   $ python lcps_batch.py /lightcurves/
  
will tell lcps to use default parameters (see help screen), search for dips in all FITS or ascii files in ``/lightcurves/`` and save the results in the default log file ``./dips.log``.

Arguments
---------
You can change the behavior of lcps's dipsearch algorithm by changing one or several of the following parameters:

=====================   =======================================================
positional arguments:
=====================   =======================================================
  path                  path containing light curve (FITS or ascii) files
=====================   =======================================================
  
  
===================   =======================================================
optional arguments:
===================   =======================================================
  -h, --help          show help message and exit
  logfile             name of log file that will contain dips
  winSize             Size of a sliding window
  stepSize            steps per slide (Default = 1, i.e. slide one data
                      point per iteration)
  Nneighb             Number of neighboring windows to be considered for the
                      local median
  minDur              minimum dip duration in # of data points
  maxDur              maximum dip duration in # of data points
  detectionThresh     fraction of flux below which a dip is registered                       
===================   =======================================================


Notation
--------
Here's how *lcps* is commanded from the shell::

   $ python lcps_batch.py [-h] [--logfile LOGFILE] [--winSize WINSIZE]
                               [--stepSize STEPSIZE] [--Nneighb NNEIGHB]
                               [--minDur MINDUR] [--maxDur MAXDUR]
                               [--detectionThresh DETECTIONTHRESH]
                               path
