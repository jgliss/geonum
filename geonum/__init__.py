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
    BASEMAP_AVAILABLE = True
    CV2_AVAILABLE = True
    NETCDF_AVAILABLE = True
    try:
        from mpl_toolkits.basemap import Basemap
    except:
        BASEMAP_AVAILABLE = False

    try:
        from cv2 import pyrUp
    except: # pragma: no cover
        CV2_AVAILABLE = False
    try:
        from netCDF4 import Dataset
    except: # pragma: no cover
        NETCDF_AVAILABLE = False

    return (BASEMAP_AVAILABLE,
            CV2_AVAILABLE,
            NETCDF_AVAILABLE)

def _init_local_topodir():
    import os
    home = os.path.expanduser('~')
    LOCAL_TOPO_DIR = os.path.join(home, '.geonum')
    if not os.path.exists(LOCAL_TOPO_DIR):
        os.mkdir(LOCAL_TOPO_DIR)
    TOPO_INFO_FILE = os.path.join(LOCAL_TOPO_DIR,  "LOCAL_TOPO_PATHS")
    if not os.path.exists(TOPO_INFO_FILE):
        with open(TOPO_INFO_FILE,'w') as f:
            f.write(f'{LOCAL_TOPO_DIR}\n')
    return (LOCAL_TOPO_DIR, TOPO_INFO_FILE)

try:
    LOCAL_TOPO_DIR, TOPO_INFO_FILE = _init_local_topodir()
except Exception as e: # pragma: no cover
    print('Failed to create local topo directory for geonum '
          f'{LOCAL_TOPO_DIR}')
    LOCAL_TOPO_DIR, TOPO_INFO_FILE =  None, None

(BASEMAP_AVAILABLE,
 CV2_AVAILABLE,
 NETCDF_AVAILABLE) = check_requirements()

def init_dir_and_version():
    import os
    from pkg_resources import get_distribution
    return (os.path.abspath(os.path.dirname(__file__)),
            get_distribution('geonum').version)

__dir__, __version__ = init_dir_and_version()

from . import exceptions
from . import helpers
from . import atmosphere
from .topodata import TopoData
from .topodataaccess import TopoDataAccess
from .topoaccessbase import delete_all_local_srtm_files


from .geopoint import GeoPoint
from .geovector3d import GeoVector3D
from .geosetup import GeoSetup
from .elevationprofile import ElevationProfile
from .processing import LineOnGrid

if BASEMAP_AVAILABLE: # pragma: no cover
    from .mapping import Map

