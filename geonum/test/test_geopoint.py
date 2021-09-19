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

if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)