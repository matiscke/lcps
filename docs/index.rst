*********************************
lcps: *light curve pre-selection*
*********************************


Introduction
------------
lcps is a tool to search for transit-like features (i.e. dips) in photometric data.

Its main purpose is to restrict large sets of light curves to a number of files that show interesting behavior (i.e. drops in flux). While lcps is adaptable to any format of time series, its current I/O module was designed specifically for photometry of the *Kepler* spacecraft. It extracts the pre-conditioned PDCSAP data from light curves files created by the standard `Kepler pipeline <https://archive.stsci.edu/k2/download_options.html>`_. It can also handle csv-formatted ascii files.

lcps uses a **sliding window** technique to compare a section of flux time series with its surroundings. A dip is detected if the flux within the window is lower than a threshold fraction of the surrounding fluxes.


Installation
------------

pip
^^^

The easiest and recommended way of installing lcps is via pip. To install the latest released version of lcps from PyPI::

   $ pip install lcps

From Source
^^^^^^^^^^^

If you prefer to use the most current development version, you can download lcps from `GitHub <https://github.com/matiscke/lcps>`_.

After unpacking the package, go to its root directory and run the setup script:
::

   $ sudo python setup.py install


Quick Start Guide
-----------------
Use the lcps_batch module to run lcps on a set of Kepler long cadence photometry, say, a complete K2 campaign::

   >>> import lcps
   >>> path = 'K2C8/'		  # path to light curve files
   >>> logfile = 'K2C8/lcps.log'  # dip detections are written to this file
   >>> lcps.lcps_batch.batchjob(path, logfile)

Usage
-----
.. toctree::
  :maxdepth: 2

  usage.rst


License and Attribution
-----------------------
Copyright 2018 Martin Schlecker

lcps is free software made available under the MIT License.

If you make use of this package in your scientific work, please acknowledge my work by citing `Schlecker, 2016 <https://zenodo.org/record/221659#.Wutd_oAiEQ8>`_.


Module Reference
----------------

.. toctree::
  :maxdepth: 1

  api


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
