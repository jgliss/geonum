#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 21:15:34 2021

@author: jonasg
"""
import pytest
import numpy as np
import numpy.testing as npt

from geonum import topodata as td
from geonum.conftest import does_not_raise_exception
arr = np.array([[647., 647.], [651., 651.]])
lats = np.array([48., 48.00083333])
lons = np.array([11., 11.00083333])

@pytest.fixture(scope='module')
def topodata():
    return td.TopoData(
        np.array([48., 48.00083333]),
        np.array([11., 11.00083333]),
        np.array([[647., 647.], [651., 651.]]),None,False)

@pytest.mark.parametrize('lats,lons,data,data_id,repl_nan_minval,raises', [
        (np.array([48., 48.00083333]), np.array([11., 11.00083333]),
         np.array([[647., 647.], [651., 651.]]),None,False,
         does_not_raise_exception()),
        (np.array([48., 48.00083333]), np.array([11., 11.00083333]),
         np.array([[647., 647.], [651., np.nan]]),None,False,
         does_not_raise_exception()),
         (np.array([48., 48.00083333]), np.array([11., 11.00083333]),
         np.array([[647., 647.], [651., np.nan]]),None,True,
         does_not_raise_exception()),
        (np.array([48., 48.00083333]), np.array([11., 11.00083333]),
         np.array([[647., 647.], [651., np.nan]]),'great_data',True,
         does_not_raise_exception()),
        ([1,2,3],[1,2],[[1,2],[1,2]],None,False,pytest.raises(ValueError)),
        ([1,2],[1,2,3],[[1,2,3],[1,2,3]],None,False,does_not_raise_exception()),
        ([1,2,3],[1,2],[[1,2,3],[1,2,3]],None,False,pytest.raises(ValueError))
        ])
def test__init__(lats, lons, data, data_id, repl_nan_minval,
                          raises):
    with raises:
        topo = td.TopoData(lats, lons, data, data_id, repl_nan_minval)
        npt.assert_array_equal(topo.lats, lats)
        npt.assert_array_equal(topo.lons, lons)

        if data_id is None:
            data_id='undefined'
        assert topo.data_id == data_id
        if repl_nan_minval:
            data[np.isnan(data)] = np.nanmin(arr)
        npt.assert_array_equal(topo.data, data)

def test_latitude(topodata):
    assert topodata.latitude is topodata.lats

def test_longitude(topodata):
    assert topodata.longitude is topodata.lons

def test_dims(topodata):
    assert topodata.dims == ['latitude', 'longitude']

def test___str__(topodata):
    assert str(topodata).startswith('geonum.TopoData')

def test___repr__(topodata):
    assert repr(topodata).startswith('geonum.TopoData')

def test_replace_nans():
    topo = td.TopoData(
        np.array([48., 48.00083333]),
        np.array([11., 11.00083333]),
        np.array([[647., 647.], [np.nan, 651.]]),None,False)
    topo.replace_nans(1000.)
    assert not np.any(np.isnan(topo.data))

def test_alt_range(topodata):
    assert topodata.alt_range == 4

def test_resolution_ERR():
    d = td.TopoData([1],[1],[[1]])
    with pytest.raises(ValueError):
        d.resolution
def test_resolution(topodata):
    np.testing.assert_allclose(topodata.resolution, (0.09,0.06), atol=0.01)

def test_mean(topodata):
    assert topodata.mean() == 649

def test_std(topodata):
    assert topodata.std() == 2
def test_delta_lon(topodata):
    np.testing.assert_allclose(topodata.delta_lon, 0.000833, atol=1e-5)

def test_delta_lat(topodata):
    np.testing.assert_allclose(topodata.delta_lat, 0.000833, atol=1e-5)

def test_center_coordinates(topodata):
    np.testing.assert_allclose(
        topodata.center_coordinates,
        (48+0.000833/2, 11+0.000833/2),
        atol=1e-5)

def test_init_mesh(topodata):
    X,Y = topodata.init_mesh()
    assert X.shape == Y.shape == (2,2)

def test_closest_index_OUTOFBOUNDS(topodata):
    with pytest.raises(ValueError):
        topodata.closest_index(0,0)

def test_closest_index(topodata):
    assert topodata.closest_index(48,11) == (0,0)

def test___call__ERR(topodata):
    with pytest.raises(ValueError):
        topodata(0,0)

def test___call__(topodata):
    assert topodata(48,11) == 647

if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)