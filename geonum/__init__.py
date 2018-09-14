# Geonum is a Python library for geographical calculations in 3D
# Copyright (C) 2017 Jonas Gliss (jonasgliss@gmail.com)
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License a
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

def check_requirements():
    try:
        from LatLon23 import LatLon
    except:
        try:
            from LatLon import LatLon
        except:
            raise ImportError('Cannot import geonum. Require either LatLon23 '
                              'or LatLon library. Use\n\npip install LatLon23\n\n'
                              'to install.')
    BASEMAP_AVAILABLE = True
    CV2_AVAILABLE = True
    NETCDF_AVAILABLE = True
    try:
        from mpl_toolkits.basemap import Basemap
    except:
        print('Plotting of maps etc. is deactivated, please install Basemap')
        BASEMAP_AVAILABLE = False
        
    
    try:
        from cv2 import pyrUp
    except:
        CV2_AVAILABLE = False
    try:
        from netCDF4 import Dataset
    except:
        NETCDF_AVAILABLE = False
    return (BASEMAP_AVAILABLE, CV2_AVAILABLE, NETCDF_AVAILABLE)
        
from os.path import abspath, dirname, join
from pkg_resources import get_distribution

BASEMAP_AVAILABLE, CV2_AVAILABLE, NETCDF_AVAILABLE = check_requirements()

__dir__ = abspath(dirname(__file__))
__version__ = get_distribution('geonum').version

_LIBDIR = __dir__ #from older version

LOCAL_TOPO_PATH = join(_LIBDIR, "local_topo_data")

from .base import GeoPoint, GeoVector3D
from .geosetup import GeoSetup
from .topodata import TopoData, TopoDataAccess
from .processing import LineOnGrid, ElevationProfile  

if BASEMAP_AVAILABLE:
    from .mapping import Map
    
from . import helpers
from . import atmosphere