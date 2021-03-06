# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A tool to search for dips in photometric time series
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
#    from .example_mod import *
    import slidingWindow
    import lcps_batch
    import lcps_io
  
