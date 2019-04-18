# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
from __future__ import (absolute_import, division)
import numpy.testing as npt
from geonum import GeoPoint, GeoVector3D

def test_GeoPoint():
    """Test basic arithmetic operations on GeoPoints."""
    p = GeoPoint(lat=0, lon=0, auto_topo_access=False)
    assert p.latitude == 0
    assert p.longitude == 0

def test_GeoVector3D():
    """Test basic arithmetic operations on GeoVectors."""
    v = GeoVector3D(dx=1, dy=1, dz=100)
    
    npt.assert_array_equal([v.dx, v.dy, v.dz], [1,1,100])
    npt.assert_array_almost_equal([v.elevation, 
                                   v.magnitude, 
                                   v.azimuth, 
                                   v.dist_hor], 
                                  [4.044691, 1.417745, 45., 1.414214])
def test_diffvector():
    p1 = GeoPoint(lat=37.751005, lon=14.993435, altitude=3264.0,
                  auto_topo_access=False)
    p2 = GeoPoint(37.765755,  15.016696, altitude=2820.0,
                  auto_topo_access=False)
    connection_vector = p2 - p1
    
    assert connection_vector.anchor is p1
    
    actual = [connection_vector.azimuth,
              connection_vector.elevation, 
              connection_vector.magnitude]
    
    
    npt.assert_allclose(actual=actual,
                        desired=[51.378677249653983,
                                 -9.6064174085658465,
                                 2.6606074796318557], rtol=1e-7)
    
if __name__ == '__main__':
    test_GeoPoint()
    test_GeoVector3D()
    test_diffvector()