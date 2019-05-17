# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
from __future__ import (absolute_import, division)
import numpy.testing as npt
from geonum.test.markers import srtm_avail, netcdf4_avail
import geonum.topodataaccess as tp

@netcdf4_avail
def test_etopo1_init():
    acc = tp.Etopo1Access(check_access=False, search_database=False)
    from geonum import LOCAL_TOPO_PATH
    assert LOCAL_TOPO_PATH == acc.local_path
    assert acc.file_name == 'ETOPO1_Ice_g_gmt4.grd'
    
@srtm_avail
def test_srtm():
    acc = tp.SRTMAccess()
    d = acc.get_data(acc._TESTLAT, acc._TESTLON)
    npt.assert_equal(d.data.shape, (2, 2))
    vals = [d.data.std(), d.data.mean(), d.data.min(), d.data.max()]
    
    npt.assert_array_almost_equal(vals, [4.0, 857.0, 853.0, 861.0])

def test_topoaccess_invalid_mode():
    all_ok = True
    try:
        tp.TopoDataAccess('bla')
    except tp.InvalidTopoMode:
        pass
    except:
        all_ok = False
    
    acc = tp.TopoDataAccess()
    try:
        acc.get_data(0,0, mode='bla')
    except tp.InvalidTopoMode:
        pass
    except:
        all_ok = False    
    assert all_ok
    
if __name__=='__main__':
    test_etopo1_init()
    test_srtm()
    test_topoaccess_invalid_mode()