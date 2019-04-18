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
    from warnings import warn
    
    BASEMAP_AVAILABLE = True
    CV2_AVAILABLE = True
    NETCDF_AVAILABLE = True
    LATLON_AVAILABLE = True
    SRTM_AVAILABLE = True
    
    try:
        from LatLon23 import LatLon
    except:
        warn('Neither LatLon23 nor LatLon are available. Many basic features '
             'will not be available (e.g. objects GeoPoint or GeoVector ')
        LATLON_AVAILABLE = False
    try:
        import srtm
    except:
        SRTM_AVAILABLE = False
    try:
        from mpl_toolkits.basemap import Basemap
    except:
        warn('Plotting of maps etc. is deactivated, please install Basemap')
        BASEMAP_AVAILABLE = False
    
    try:
        from cv2 import pyrUp
    except:
        CV2_AVAILABLE = False
    try:
        from netCDF4 import Dataset
    except:
        NETCDF_AVAILABLE = False
        
    return (LATLON_AVAILABLE, 
            SRTM_AVAILABLE, 
            BASEMAP_AVAILABLE, 
            CV2_AVAILABLE, 
            NETCDF_AVAILABLE)
        
from os.path import abspath, dirname, join
from pkg_resources import get_distribution

(LATLON_AVAILABLE, 
 SRTM_AVAILABLE, 
 BASEMAP_AVAILABLE, 
 CV2_AVAILABLE, 
 NETCDF_AVAILABLE) = check_requirements()

__dir__ = abspath(dirname(__file__))
__version__ = get_distribution('geonum').version

_LIBDIR = __dir__ #from older version

LOCAL_TOPO_PATH = join(_LIBDIR, "local_topo_data")

TOPO_INFO_FILE = join(LOCAL_TOPO_PATH,  "LOCAL_TOPO_PATHS.txt")

from . import exceptions
from . import helpers
from . import atmosphere
from .topodata import TopoData
from .topodataaccess import TopoDataAccess
from .topoaccessbase import delete_all_local_srtm_files

if LATLON_AVAILABLE:
    from .geopoint import GeoPoint
    from .geovector3d import GeoVector3D
    from .geosetup import GeoSetup
    from .elevationprofile import ElevationProfile
    from .processing import LineOnGrid  

if BASEMAP_AVAILABLE:
    from .mapping import Map
    
