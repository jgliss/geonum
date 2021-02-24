# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
import numpy.testing as npt
import pytest
from geonum.conftest import skip_srtm
from geonum import TopoDataAccess

@skip_srtm
@pytest.mark.parametrize('lat0, lon0, lat1, lon1, dsh, dmin, dmax, dmean', [
    (-18.55, -69.2, -18.35, -69.0, (242, 122), 4084., 6057., 4690.692827)
    ])
def test_srtm_access(lat0, lon0, lat1, lon1, dsh, dmin, dmax, dmean):
    acc = TopoDataAccess(mode='srtm')
    data = acc.get_data(lat0, lon0, lat1, lon1)
    npt.assert_array_equal(data.shape, dsh)
    npt.assert_array_almost_equal([data.min, data.max, data.mean()],
                                  [dmin, dmax, dmean])


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)



