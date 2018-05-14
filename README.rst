lcps (light curve pre-selection)
================================

**A tool to search for dips in photometric time series**

.. image:: https://img.shields.io/badge/ascl-1805.003-blue.svg?colorB=262255
    :target: http://ascl.net/1805.003


lcps searches for transit-like features in light curve files.



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
   

Usage Example
-------------
::

  $ python lcps_batch.py /KeplerData/C7/ --logfile /KeplerData/C7/lcps_long+deep01.log --winSize 700 --stepSize 20 --minDur 20 --maxDur 698 --detectionThresh 0.90


Documentation
-------------
A full documentation of lcps is available under `<http://lcps.readthedocs.io>`_.


License
-------
Copyright 2016 Martin Schlecker

lcps is free software made available under the MIT License. For details see
the LICENSE file.
