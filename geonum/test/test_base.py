# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
from __future__ import (absolute_import, division)
from geonum.base import GeoPoint, GeoVector3D


def test_GeoPoint():
    """Test basic arithmetic operations on GeoPoints."""
    assert(GeoPoint(1, 1))


def test_GeoVector3D():
    """Test basic arithmetic operations on GeoVectors."""
    assert(GeoVector3D(1, 1))
