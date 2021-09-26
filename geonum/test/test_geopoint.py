#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 13:55:05 2021

@author: jonasg
"""

import pytest
import numpy as np
from geonum.conftest import does_not_raise_exception
from LatLon23 import LatLon
from geonum.geopoint import GeoPoint
from geonum.geovector3d import GeoVector3D
from geonum.exceptions import OutOfDomain
from geonum import TopoData

FAKE_TOPO = TopoData(lats=[-0.5,0.5],
                     lons=[-0.5,0.5],
                     data=[[10,10],[10,10]])

def test_GeoPoint__ALTERR_DEFAULT():
    assert GeoPoint._ALTERR_DEFAULT == 1e9

def test_GeoPoint__init__empty():
    gp = GeoPoint()
    assert isinstance(gp, LatLon)
    assert gp.name == 'undefined'
    assert gp.topo_access_mode == None
    assert gp.topo_data == None
    assert gp.local_topo_path == None
    assert gp.longitude == 0
    assert gp.latitude == 0
    assert gp.altitude == 0 # no SRTM over seas
    assert gp.altitude_err == GeoPoint._ALTERR_DEFAULT


@pytest.mark.parametrize('args,eval_args,raises', [
    ({},{}, does_not_raise_exception()),
    ({'altitude' : 10, 'auto_topo_access' : True}, {},
     pytest.raises(ValueError)),
    ({'altitude' : 10}, {'altitude' : 10},
     does_not_raise_exception()),
    ({'topo_data' : 43}, {}, pytest.raises(ValueError)),
    ({'topo_data' : FAKE_TOPO, 'auto_topo_access' : True},
     {'altitude' : 10}, does_not_raise_exception())
    ])
def test_GeoPoint__init__(args, eval_args, raises):
    with raises:
        gp = GeoPoint(**args)
        assert isinstance(gp, GeoPoint)
        for arg, val in eval_args.items():
            assert getattr(gp, arg) == val

@pytest.mark.parametrize('topo_data,gp,raises', [
    ('blaa', GeoPoint(), pytest.raises(ValueError)),
    (FAKE_TOPO, GeoPoint(), does_not_raise_exception()),
    (FAKE_TOPO, GeoPoint(10,10), pytest.raises(OutOfDomain))
    ])
def test_GeoPoint_set_topo_data(topo_data, gp, raises):
    with raises:
        gp.set_topo_data(topo_data)
        assert gp.topo_data == topo_data

@pytest.mark.parametrize('points,ll,tr', [
    ((), [44.994,14.991],[45.006,15.009]),
    ([GeoPoint(46, 16)], [44.913,14.878],[46.086,16.124]),
    ([GeoPoint(46, 16), GeoPoint(-20, -80)], [-27.464, -88.592],
     [53.016, 28.711])
])
def test_GeoPoint_range_borders(points,ll,tr):
    p0 = GeoPoint(45, 15)
    pll, ptr = p0.range_borders(*points)
    assert isinstance(pll, GeoPoint)
    assert isinstance(ptr, GeoPoint)
    _ll = [pll.latitude, pll.longitude]
    _tr = [ptr.latitude, ptr.longitude]
    np.testing.assert_allclose(_ll, ll, atol=0.01)
    np.testing.assert_allclose(_tr, tr, atol=0.01)

@pytest.mark.parametrize('kwargs,pf', [
    ({}, (45,15)),
    (dict(lat1=46, lon1=16), (46,16)),
    (dict(geo_point=GeoPoint(46,16)), (46,16)),
    (dict(azimuth=45,dist_hor=111*np.sqrt(2)),(45.99,16.433))
])
def test_GeoPoint_get_topo_data(kwargs,pf):
    p0 = GeoPoint(45, 15)
    topo, p1 = p0.get_topo_data(**kwargs)
    assert isinstance(topo, TopoData)
    assert isinstance(p1, GeoPoint)
    np.testing.assert_allclose([p1.latitude, p1.longitude], pf, atol=1e-2)


@pytest.mark.parametrize('topo,lat1,lon1,result', [
    (None,None,None,False),
    (FAKE_TOPO,None,None,True),
    (FAKE_TOPO,0.1,0.1,True),
    (FAKE_TOPO,10,20,False),
])
def test_GeoPoint_check_topo(topo,lat1,lon1,result):
    pt = GeoPoint(-0.1,-0.1,topo_data=topo)
    val = pt.check_topo(lat1,lon1)
    assert result == val

@pytest.mark.parametrize('kwargs,raises', [
    ({}, pytest.raises(ValueError)),
    (dict(geo_point=GeoPoint(46,15)), does_not_raise_exception()),
])
def test_GeoPoint_get_elevation_profile(kwargs,raises):
    from geonum.elevationprofile import ElevationProfile
    p = GeoPoint(45,15)
    with raises:
        ep = p.get_elevation_profile(**kwargs)
        assert isinstance(ep, ElevationProfile)

def test_GeoPoint__sub_geo_vector_2d():
    from LatLon23 import GeoVector
    pt = GeoPoint(0,1)
    v = GeoVector(dx=1,dy=0)
    result = pt._sub_geo_vector_2d(v)
    assert isinstance(result, GeoPoint)

def test_GeoPoint__add_geo_vector_2d():
    from LatLon23 import GeoVector
    pt = GeoPoint(0,1)
    v = GeoVector(dx=1,dy=0)
    result = pt._add_geo_vector_2d(v)
    assert isinstance(result, GeoPoint)

def test_GeoPoint__sub_geo_vector_3d():
    pt = GeoPoint(0,1, altitude=1000)
    v = GeoVector3D(dx=1,dy=1, dz=100)
    result = pt._sub_geo_vector_3d(v)
    assert isinstance(result, GeoPoint)

def test_GeoPoint__add_geo_vector_3d():
    pt = GeoPoint(0,1, altitude=1000)
    v = GeoVector3D(dx=1,dy=1, dz=100)
    result = pt._add_geo_vector_3d(v)
    assert isinstance(result, GeoPoint)

def test_GeoPoint__sub_latlon():
    from LatLon23 import LatLon
    pt = GeoPoint(0,1, altitude=1000)
    p2 = LatLon(-1, 0)
    result = pt._sub_latlon(p2)
    assert isinstance(result, GeoVector3D)
