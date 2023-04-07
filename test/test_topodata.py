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
        ('blaa',[1,2],[[1,2],[1,2]],None,False,pytest.raises(ValueError)),
        ([1,2,3],[1,2],[[1,2],[1,2]],None,False,pytest.raises(ValueError)),
        ([1,2],[1,2,3],[[1,2,3],[1,2,3]],None,False,does_not_raise_exception()),
        ([1,2,3],[1,2],[[1,2,3],[1,2,3]],None,False,pytest.raises(ValueError))
        ])
def test_TopoData__init__(lats, lons, data, data_id, repl_nan_minval,
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

def test_TopoData_latitude(topodata):
    assert topodata.latitude is topodata.lats

def test_TopoData_longitude(topodata):
    assert topodata.longitude is topodata.lons

def test_TopoData_dims(topodata):
    assert topodata.dims == ['latitude', 'longitude']

def test_TopoData___str__(topodata):
    assert str(topodata).startswith('geonum.TopoData')

def test_TopoData___repr__(topodata):
    assert repr(topodata).startswith('geonum.TopoData')

def test_TopoData_replace_nans():
    topo = td.TopoData(
        np.array([48., 48.00083333]),
        np.array([11., 11.00083333]),
        np.array([[647., 647.], [np.nan, 651.]]),None,False)
    topo.replace_nans(1000.)
    assert not np.any(np.isnan(topo.data))

if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)