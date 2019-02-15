# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
from __future__ import (absolute_import, division)
import numpy.testing as npt

import geonum.topodata as tp


def test_etopo1_init():
    acc = tp.Etopo1Access(check_access=False, search_database=False)
    from geonum import LOCAL_TOPO_PATH
    assert LOCAL_TOPO_PATH == acc.local_path
    assert acc.file_name == 'ETOPO1_Ice_g_gmt4.grd'
    
def test_srtm():
    acc = tp.SRTMAccess()
    d = acc.get_data(acc._TESTLAT, acc._TESTLON)
    npt.assert_equal(d.data.shape, (2, 2))
    vals = [d.data.std(), d.data.mean(), d.data.min(), d.data.max()]
    
    npt.assert_array_almost_equal(vals, [4.0, 857.0, 853.0, 861.0])

    
if __name__=='__main__':
    test_etopo1_init()
    test_srtm()
