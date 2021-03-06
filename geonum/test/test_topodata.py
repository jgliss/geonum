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
arr = np.array([[647., 647.], [651., 651.]])
lats = np.array([48., 48.00083333])
lons = np.array([11., 11.00083333])

@pytest.fixture(scope='module')
def topodata():
    return td.TopoData(
        np.array([48., 48.00083333]),
        np.array([11., 11.00083333]),
        np.array([[647., 647.], [651., 651.]]),None,False)

@pytest.mark.parametrize(
    'lats, lons, data, data_id, repl_nan_minval', [
        (np.array([48., 48.00083333]), np.array([11., 11.00083333]),
         np.array([[647., 647.], [651., 651.]]),None,False),
        (np.array([48., 48.00083333]), np.array([11., 11.00083333]),
         np.array([[647., 647.], [651., np.nan]]),None,False),
         (np.array([48., 48.00083333]), np.array([11., 11.00083333]),
         np.array([[647., 647.], [651., np.nan]]),None,True),
        (np.array([48., 48.00083333]), np.array([11., 11.00083333]),
         np.array([[647., 647.], [651., np.nan]]),'great_data',True),
        ])
def test_TopoData__init__(lats, lons, data, data_id, repl_nan_minval):
    topo = td.TopoData(lats, lons, data, data_id, repl_nan_minval)
    npt.assert_array_equal(topo.lats, lats)
    npt.assert_array_equal(topo.lons, lons)

    if data_id is None:
        data_id=''
    assert topo.data_id == data_id
    if repl_nan_minval:
        data[np.isnan(data)] = np.nanmin(arr)
    npt.assert_array_equal(topo.data, data)

def test_TopoData_latitude(topodata):
    assert topodata.latitude is topodata.lats

def test_TopoData_longitude(topodata):
    assert topodata.longitude is topodata.lons

def test_TopoData_replace_nans():
    topo = td.TopoData(
        np.array([48., 48.00083333]),
        np.array([11., 11.00083333]),
        np.array([[647., 647.], [np.nan, 651.]]),None,False)
    topo.replace_nans(1000.)
    assert not np.any(np.isnan(topo.data))

if __name__=='__main__':

    data = td.TopoData(lats, lons, arr)


    import sys
    pytest.main(sys.argv)