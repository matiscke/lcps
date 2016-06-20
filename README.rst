lcps (light curve pre-selection)
================================

**A tool to search for dips in photometric time series**

lcps searches for transit-like features in light curve files.


Usage Example
-------------
python lcps_batch.py /KeplerData/C7/ --logfile /KeplerData/C7/lcps_long+deep01.log --winSize 700 --stepSize 20 --minDur 20 --maxDur 698 --detectionThresh 0.90


License
-------

Copyright 2016 Martin Schlecker

lcps is free software made available under the MIT License. For details see
the LICENSE file.