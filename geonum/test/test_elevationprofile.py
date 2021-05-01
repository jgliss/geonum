# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
import numpy.testing as npt
import pytest
from geonum.conftest import skip_srtm, does_not_raise_exception
import numpy as np
from geonum import GeoPoint, ElevationProfile, TopoData

LAT0, LON0 = -39.296571, 173.9224
LAT1, LON1 = -39.3538, 174.4383

@pytest.fixture(scope='module')
def profile():
    return ElevationProfile(GeoPoint(LAT0, LON0), GeoPoint(LAT1, LON1),
                            calc_on_init=False)

@pytest.mark.filterwarnings("ignore:Failed to compute elevation profile.")
def test_ElevationProfile_wrong_input():
    with pytest.raises(ValueError):
        ElevationProfile('bla', 1, 'blub')

@pytest.mark.parametrize('args,raises', [
    ({'observer': 1}, pytest.raises(TypeError)),
    ({'observer': 1, 'endpoint': 2}, pytest.raises(ValueError)),
    ({'observer': GeoPoint(46,7),
      'endpoint': GeoPoint(46.1,7.1),
      'calc_on_init': False}, does_not_raise_exception())
    ])
def test_ElevationProfile___init__(args,raises):
    with raises:
        prof = ElevationProfile(**args)
        assert isinstance(prof, ElevationProfile)

def test_ElevationProfile_observer(profile):
    # test getter
    val = profile.observer
    assert isinstance(val, GeoPoint)
    # test setter
    profile.observer = val
    assert profile.observer is val

def test_ElevationProfile_endpoint(profile):
    # test getter
    val = profile.endpoint
    assert isinstance(val, GeoPoint)
    # test setter
    profile.endpoint = val
    assert profile.endpoint is val

def test_ElevationProfile_topo_data(profile):
    # test getter
    val = profile.topo_data
    assert isinstance(val, TopoData)
    # test setter
    profile.topo_data = val
    assert profile.topo_data is val

def test_ElevationProfile_observer_topogrid(profile):
    profile._coords_topo = None
    val = profile.observer_topogrid
    assert isinstance(val, GeoPoint)

def test_ElevationProfile_endpoint_topogrid(profile):
    profile._coords_topo = None
    val = profile.endpoint_topogrid
    assert isinstance(val, GeoPoint)

def test_ElevationProfile_azimuth(profile):
    val = profile.azimuth
    npt.assert_allclose(val, 98.3, atol=0.1)

def test_ElevationProfile_profile_unavail(profile):
    with pytest.raises(AttributeError):
        profile.profile

def test_ElevationProfile_dists_unavail(profile):
    with pytest.raises(AttributeError):
        profile.dists

@pytest.mark.parametrize('args,num,avg', [
    ({}, 9222, 576),
    ({'interpolate' : False},319, 575),
    ({'interpolate' : True, 'resolution' : 1000},319, 575),
    ])
def test_ElevationProfile_det_profile(profile, args, num, avg):
    val = profile.det_profile(**args)
    assert isinstance(val, np.ndarray)
    assert len(val) == num
    mean = np.nanmean(val)
    npt.assert_allclose(mean, avg, atol=1)

def test_ElevationProfile_resolution(profile):
    npt.assert_allclose(profile.resolution, 0.14, atol=0.1)

def test_ElevationProfile_gradient(profile):
    grad = profile.gradient
    assert len(grad) == 319
    npt.assert_allclose(np.mean(grad), -0.43, atol=0.1)

@pytest.mark.parametrize('elev_angle,view_above_topo_m,num,avg', [
    (0, 0, 319, 318),
    (0, 10, 319, 328),
    (2, 10, 319, 1114),
    ])
def test_ElevationProfile_get_altitudes_view_dir(profile, elev_angle,
                                                 view_above_topo_m,num,avg):
    alts = profile.get_altitudes_view_dir(elev_angle, view_above_topo_m)
    assert isinstance(alts, np.ndarray)
    assert len(alts) == num
    mean = np.mean(alts)
    npt.assert_allclose(mean, avg, atol=1)

@skip_srtm
def test_ElevationProfile_via_GeoPoint():


    obs = GeoPoint(LAT0, LON0,
                   auto_topo_access=True)

    npt.assert_almost_equal(obs.altitude, 324)

    prof = obs.get_elevation_profile(azimuth=90, dist_hor=50)

    prof_nans = obs.get_elevation_profile(azimuth=45, dist_hor=50)
    prof2 = obs.get_elevation_profile(azimuth=45, dist_hor=50, order=1)

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
