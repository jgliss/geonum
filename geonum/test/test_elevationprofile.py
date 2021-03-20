# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
import numpy.testing as npt
import pytest
from geonum.conftest import skip_srtm
import numpy as np
from geonum import GeoPoint, ElevationProfile

@pytest.mark.filterwarnings("ignore:Failed to compute elevation profile.")
def test_ElevationProfile_wrong_input():
    ElevationProfile('bla', 1, 'blub')

@skip_srtm
def test_ElevationProfile():

    lat_taranaki = -39.296571
    lon_observer = 173.9224

    obs = GeoPoint(lat_taranaki, lon_observer,
                   auto_topo_access=True)

    npt.assert_almost_equal(obs.altitude, 324)

    prof = obs.get_elevation_profile(azimuth=90, dist_hor=50)

    prof_nans = obs.get_elevation_profile(azimuth=45, dist_hor=50)
    prof2 = obs.get_elevation_profile(azimuth=45, dist_hor=50,
                                      order=1)

    npt.assert_array_equal([len(prof.dists), len(prof.profile),
                            np.sum(np.isnan(prof_nans.profile))],
                            [10092,10092,10263])
    npt.assert_array_almost_equal([prof.profile.max(), prof.profile.min(),
                                   prof.profile.mean(), prof.profile.std(),
                                   prof.dists.max(),
                                   prof2.profile.max(), prof2.profile.min()],
                                  [2482.598405, 176.008484, 583.012202,
                                   486.57815, 50.033677, 1005.208932, 0.])

if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)
