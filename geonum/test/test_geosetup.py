#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 21:15:34 2021

@author: jonasg
"""
import pytest
from geonum import geosetup as gs
from geonum.conftest import skip_srtm

@pytest.mark.parametrize(
    'points,vectors,lat_ll,lon_ll,lat_tr,lon_tr,id,topo_access_mode, local_topo_path,cmap_vecs',[
        (None,None,None,None,None,None,None,None,None,None)

        ])
def test_GeoSetup__init__(points,vectors,lat_ll,lon_ll,lat_tr,lon_tr,
                          id,topo_access_mode, local_topo_path,
                          cmap_vecs):
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
def test_GeoSetup_create_test_data():
    stp = gs.GeoSetup()
    stp.create_test_data()
    assert len(stp.points) == 4
    assert len(stp.vectors) == 2

if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)