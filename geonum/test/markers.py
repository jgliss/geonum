# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
from __future__ import (absolute_import, division)
import pytest
from geonum import SRTM_AVAILABLE

# custom skipif marker that is used below for test functions that 
# require SRTM.py to be installed
srtm_avail = pytest.mark.skipif(SRTM_AVAILABLE==False,
                   reason='Skipping SRTM database tests. srtm.py library is '
                   'not installed')
                    #allow_module_level=True)