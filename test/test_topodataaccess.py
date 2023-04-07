# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
import pytest
import numpy.testing as npt
from geonum.conftest import skip_srtm, skip_netcdf4
import geonum.topodataaccess as tp

@skip_netcdf4
def test_etopo1_init():
    acc = tp.Etopo1Access(check_access=False, search_database=False)
    from geonum import LOCAL_TOPO_DIR
    assert LOCAL_TOPO_DIR == acc.local_path
    assert acc.file_name == 'ETOPO1_Ice_g_gmt4.grd'

@skip_srtm
def test_srtm():
    acc = tp.SRTMAccess()
    d = acc.get_data(acc._TESTLAT, acc._TESTLON)
    npt.assert_equal(d.data.shape, (2, 2))
    vals = [d.data.std(), d.data.mean(), d.data.min(), d.data.max()]

    npt.assert_array_almost_equal(vals, [4.0, 857.0, 853.0, 861.0])

@skip_srtm
@pytest.mark.parametrize('lat0, lon0, lat1, lon1, dsh, dmin, dmax, dmean', [
    (60.052, 7.414, None, None, (2, 2), 1153., 1183., 1167.5),
    (-18.55, -69.2, -18.35, -69.0, (242, 122), 4084., 6057., 4690.692827),

    ])
def test_srtm_access(lat0, lon0, lat1, lon1, dsh, dmin, dmax, dmean):
    acc = tp.TopoDataAccess(mode='srtm')
    data = acc.get_data(lat0, lon0, lat1, lon1)
    npt.assert_array_equal(data.shape, dsh)
    npt.assert_array_almost_equal([data.min, data.max, data.mean()],
                                  [dmin, dmax, dmean])

def test_topoaccess_invalid_mode():
    all_ok = True
    try:
        tp.TopoDataAccess('bla')
    except tp.InvalidTopoMode:
        pass
    except Exception:
        all_ok = False

    acc = tp.TopoDataAccess()
    try:
        acc.get_data(0,0, mode='bla')
    except tp.InvalidTopoMode:
        pass
    except Exception:
        all_ok = False
    assert all_ok

if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)