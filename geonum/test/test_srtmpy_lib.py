#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:53:56 2021

@author: jonasg
"""
import os
import pytest
import srtm

def test_filehandler_access():
    fh = srtm.utils.FileHandler()

def test_filehandler_local_cache_dir():
    fh = srtm.utils.FileHandler()
    assert 'local_cache_dir' in fh.__dict__
    assert not fh.exists('nonexistingfile.hgt')
    assert os.path.exists(fh.local_cache_dir)

def test_get_data():
    acc = srtm.get_data()
    assert isinstance(acc, srtm.data.GeoElevationData)

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
