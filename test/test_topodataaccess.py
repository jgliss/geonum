# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
import pytest
import numpy.testing as npt
from geonum.conftest import skip_srtm, skip_netcdf4
import geonum.topodataaccess as mod

@pytest.fixture
def acc():
    return mod.TopoDataAccess()

@skip_netcdf4
def test_etopo1_init():
    acc = mod.Etopo1Access(check_access=False, search_database=False)
    from geonum import LOCAL_TOPO_DIR
    assert LOCAL_TOPO_DIR == acc.local_path
    assert acc.file_name == 'ETOPO1_Ice_g_gmt4.grd'

@skip_srtm
def test_srtm():
    acc = mod.SRTMAccess()
    d = acc.get_data(acc._TESTLAT, acc._TESTLON)
    npt.assert_equal(d.data.shape, (1, 1))
    assert d.data[0][0] == 861

@skip_srtm
@pytest.mark.parametrize('lat0, lon0, lat1, lon1, dsh, dmin, dmax, dmean', [
    (60.052, 7.414, None, None, (2, 2), 1153., 1183., 1167.5),
    (-18.55, -69.2, -18.35, -69.0, (241, 121), 4084., 6057., 4690.5),

    ])
def test_srtm_access(lat0, lon0, lat1, lon1, dsh, dmin, dmax, dmean):
    acc = mod.TopoDataAccess(mode='srtm')
    data = acc.get_data(lat0, lon0, lat1, lon1)
    npt.assert_array_equal(data.shape, dsh)
    npt.assert_allclose([data.min, data.max, data.mean()],
                                  [dmin, dmax, dmean], atol=0.1)

def test_topoaccess_invalid_mode():
    all_ok = True
    try:
        mod.TopoDataAccess('bla')
    except mod.InvalidTopoMode:
        pass
    except Exception:
        all_ok = False

    acc = mod.TopoDataAccess()
    try:
        acc.get_data(0,0, mode='bla')
    except mod.InvalidTopoMode:
        pass
    except Exception:
        all_ok = False
    assert all_ok

def test_modes(acc):
    assert acc.modes == ['srtm', 'etopo1']

if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)