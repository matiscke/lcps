lcps: *light curve pre-selection*
=================================


Introduction
------------
lcps is a tool to search for transit-like features (i.e. dips) in photometric data. 

Its main purpose is to restrict large sets of light curves to a number of files that show interesting behavior (i.e. drops in flux). While lcps is adaptable to any format of time series, its current I/O module was designed specifically for photometry of the *Kepler* spacecraft. It extracts the pre-conditioned PDCSAP data from light curves files created by the standard `Kepler pipeline <https://archive.stsci.edu/k2/download_options.html>`_.

lcps uses a **sliding window** technique to compare a section of flux time series with its surroundings. A dip is detected if the flux within the window is lower than a threshold fraction of the surrounding fluxes. 

Installation
------------
The latest version of lcps can be downloaded from github or installed from PyPI via

  $ pip install lcps



Quick Start Guide
-----------------
Use the lcps_batch module to run lcps on a set of Kepler long cadence photometry, say, a complete K2 campaign: 
   >>> import lcps   
   >>> path = 'K2C8/'		  # path to light curve files
   >>> logfile = 'K2C8/lcps.log'  # dip detections are written to this file  
   >>> lcps.lcps_batch.batchjob(path, logfile)
 
Copyright and License
---------------------
Copyright 2016 Martin Schlecker

lcps is free software made available under the MIT License.



Module Reference
----------------

.. toctree::
  :maxdepth: 1

  api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`