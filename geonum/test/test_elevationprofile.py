# -*- coding: utf-8 -*-
"""Test environment for base.py module."""
from __future__ import (absolute_import, division)
import numpy.testing as npt
import pytest
import numpy as np
from geonum import SRTM_AVAILABLE
from geonum import GeoPoint

@pytest.mark.skipif(SRTM_AVAILABLE==False,
                   reason='Skipping Elevation profile tests. srtm.py library is '
                   'not installed')
    
def test_ElevationProfile():
    
    lat_taranaki = -39.296571
    lon_observer = 173.9224
    
    obs = GeoPoint(lat_taranaki, lon_observer)
    
    npt.assert_almost_equal(obs.altitude, 324)
    
    prof = obs.get_elevation_profile(azimuth=90, dist_hor=50)
    
    prof_nans = obs.get_elevation_profile(azimuth=45, dist_hor=50)
    prof2 = obs.get_elevation_profile(azimuth=45, dist_hor=50,
                                      order=1)
    
    npt.assert_array_equal([len(prof.dists), len(prof.profile),
                            np.sum(np.isnan(prof_nans.profile))],
                            [10092,10092,10263])
    npt.assert_array_almost_equal([prof.profile.max(), prof.profile.min(), 
                                   prof.profile.mean(), prof.profile.std(),
                                   prof.dists.max(),
                                   prof2.profile.max(), prof2.profile.min()],
                                  [2482.598405, 176.008484, 583.012202,  
                                   486.57815, 50.033677, 1005.208932, 0.])
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    plt.close('all')
    test_ElevationProfile()
    
# =============================================================================
#     prof.line.plot_line_on_grid(prof.topo_data.data)
#     
#     
#     d = prof.topo_data.data
#     
#     # all NaN
#     
#     z1 = prof.line.get_line_profile(d, prefilter=False)
#     z2 = prof.line.get_line_profile(d, order=1, prefilter=True)
#     z3 = prof.line.get_line_profile(d, order=1, prefilter=False)
#     
#     z0 = prof.line.get_line_profile(d)
#     
#     
#     d[np.isnan(d)] = 0
#     
#     zi = prof.line.get_line_profile(d)
#     
#     fig = plt.figure(figsize=(16, 8))
#     
#     
#     plt.plot(z1, label='prefilter=False')
#     plt.plot(z2, label='order=1, prefilter=True')
#     plt.plot(z3, ' x', label='order=1')
#     plt.plot(zi, label='default, after NaNs -> 0')
#     ax = fig.axes[0]
#     ax.set_yscale("log", nonposy='clip')
#     ax.set_ylim([0.1, 3000])
#     
#     #
#     
#     plt.legend()
# =============================================================================
