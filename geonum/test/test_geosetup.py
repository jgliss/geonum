#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 21:15:34 2021

@author: jonasg
"""
import pytest
from geonum.conftest import does_not_raise_exception
from geonum import GeoPoint, GeoVector3D, TopoDataAccess
from geonum import geosetup as gs
from geonum.conftest import skip_srtm

@pytest.mark.parametrize(
    'points,vectors,lat_ll,lon_ll,lat_tr,lon_tr,id,topo_access_mode, local_topo_path,cmap_vecs,raises',[
        (None,None,None,None,None,None,None,None,None,None, does_not_raise_exception()),
        (GeoPoint(),None,None,None,None,None,None,None,None,None, does_not_raise_exception()),
        ('Blaaaa',None,None,None,None,None,None,None,None,None, pytest.raises(ValueError)),
        (None,GeoVector3D(0,0),None,None,None,None,None,None,None,None, does_not_raise_exception()),
        (None, 'Blaaaa',None,None,None,None,None,None,None,None, pytest.raises(ValueError)),

        ])
def test_GeoSetup__init__(points,vectors,lat_ll,lon_ll,lat_tr,lon_tr,
                          id,topo_access_mode, local_topo_path,
                          cmap_vecs,raises):
    with raises:
        stp=gs.GeoSetup(points,vectors,lat_ll,lon_ll,lat_tr,lon_tr,
                              id,topo_access_mode, local_topo_path,
                              cmap_vecs)
        if id is None:
            id ='MyGeoSetup'
        if topo_access_mode is None:
            topo_access_mode = "srtm"
        if cmap_vecs is None:
            cmap_vecs = 'Greens'
        if points is None:
            points = []
        if vectors is None:
            vectors = []
        assert stp.id == id
        assert stp.topo_access_mode == topo_access_mode
        assert stp.local_topo_path == local_topo_path

@skip_srtm
def test_GeoSetup_create_test_setup():
    stp = gs.GeoSetup.create_test_setup()
    assert isinstance(stp, gs.GeoSetup)
    assert len(stp.points) == 4
    assert len(stp.vectors) == 2

@pytest.fixture(scope='module')
def empty_setup():
    return gs.GeoSetup()

def test_GeoSetup_cmap_vecs(empty_setup):
    assert empty_setup._cmap_vecs == 'Greens'
    from matplotlib.colors import LinearSegmentedColormap
    assert isinstance(empty_setup.cmap_vecs, LinearSegmentedColormap)

def test_GeoSetup_topo_access(empty_setup):
    assert isinstance(empty_setup.topo_access, TopoDataAccess)

def test_has_points(empty_setup):
    assert not empty_setup.has_points()
    empty_setup.add_geo_point(GeoPoint(name='null'), assert_in_domain=False)
    assert empty_setup.has_points()
    empty_setup.delete_geo_point('null')
    assert not empty_setup.has_points()


if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)