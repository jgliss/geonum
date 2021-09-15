# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
import numpy.testing as npt
import pytest
from geonum.conftest import skip_srtm, does_not_raise_exception
import numpy as np
from geonum import GeoPoint, ElevationProfile, TopoData
from matplotlib.axes import Axes

LAT0, LON0 = -39.296571, 173.9224
LAT1, LON1 = -39.3538, 174.4383

P0 = GeoPoint(LAT0, LON0)
P1 = GeoPoint(LAT1, LON1)

@pytest.fixture(scope='module')
def profile():
    return ElevationProfile(P0, P1,
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

def test_ElevationProfile_dists_fail():
    with pytest.raises(AttributeError):
        p = ElevationProfile(P0,P1,calc_on_init=False)
        p.dists
    p.det_profile()
    assert isinstance(p.dists, np.ndarray)

def test_ElevationProfile_resolution(profile):
    profile.det_profile()
    res = profile.resolution
    assert isinstance(res, float)
    npt.assert_allclose(res, 0.005, atol=0.01)

def test_ElevationProfile_gradient(profile):
    profile.det_profile()
    grad = profile.gradient
    assert isinstance(grad, np.ndarray)
    assert len(grad) == 9222
    npt.assert_allclose(np.mean(grad), -0.015, atol=0.01)

def test_ElevationProfile_slope(profile):
    profile.det_profile()
    sl = profile.slope
    assert isinstance(sl, np.ndarray)
    npt.assert_allclose(np.mean(sl), -0.003, atol=0.001)


@pytest.mark.parametrize('elev_angle,view_above_topo_m,num,avg', [
    (0, 0, 319, 318),
    (0, 10, 319, 328),
    (2, 10, 319, 1114),
    ])
def test_ElevationProfile_get_altitudes_view_dir(profile, elev_angle,
                                                 view_above_topo_m,num,avg):
    profile.det_profile(**{'interpolate': True, 'resolution': 1000})
    alts = profile.get_altitudes_view_dir(elev_angle, view_above_topo_m)
    assert isinstance(alts, np.ndarray)
    assert len(alts) == num
    mean = np.mean(alts)
    npt.assert_allclose(mean, avg, atol=1)

@pytest.mark.parametrize(
    'elev_angle,view_above_topo_m,min_dist,local_tolerance,plot,d,derr,alt', [
    (0,0,0,3,True,31.1,0.3,319.2),
    (0,0,0,10,True,31.1,0.9,319.2),
    (6,0,0,10,True,10.8,0.7,1399.7),
    (10,0,0,10,True,np.nan,np.nan,None)
    ]
    )
def test_ElevationProfile_get_first_intersection(profile,
    elev_angle,view_above_topo_m,min_dist,local_tolerance,plot,
    d,derr,alt):
    profile.det_profile(**{'interpolate': True, 'resolution': 1000})
    val = profile.get_first_intersection(elev_angle,view_above_topo_m,min_dist,
                                         local_tolerance,plot)
    dist, dist_err, intersect, view_elevations, ax = val

    assert isinstance(view_elevations, np.ndarray)
    if np.isnan(d):
        assert intersect is None
        assert np.isnan(dist)
        assert np.isnan(dist_err)
    else:
        assert isinstance(intersect, GeoPoint)
        npt.assert_allclose([dist,dist_err,intersect.altitude], [d, derr, alt],
                        atol=0.1)
    if plot:
        assert isinstance(ax, Axes)

@pytest.mark.parametrize(
    'elev_start,elev_stop,step_deg,raises,num,elev', [
    (0,10,0.1,does_not_raise_exception(),69,6.9),
    (6.8,7,0.01,does_not_raise_exception(),10,6.9)
    ])
def test_ElevationProfile_find_horizon_elev(profile,
    elev_start,elev_stop,step_deg,raises,num,elev):
    with raises:
        (_elev,
         elev_sects,
         dists_sects) = profile.find_horizon_elev(elev_start, elev_stop,
                                                  step_deg)

        npt.assert_allclose(_elev, elev, atol=0.1)
        assert len(elev_sects) == len(dists_sects) == num

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
