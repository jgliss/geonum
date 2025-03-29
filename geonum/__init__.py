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
from pathlib import Path as _Path
from importlib import metadata
from ._init_helpers import _check_requirements, _init_local_topodir

LOCAL_TOPO_DIR, TOPO_INFO_FILE = _init_local_topodir()

(BASEMAP_AVAILABLE,
 CV2_AVAILABLE,
 NETCDF_AVAILABLE) = _check_requirements()

__dir__ = str(_Path(__file__).parent)
__version__ = metadata.version("geonum")

from . import exceptions
from . import helpers
from . import atmosphere
from . import plot_helpers

from .topodata import TopoData
from .topodataaccess import TopoDataAccess

from .geopoint import GeoPoint
from .geovector3d import GeoVector3D
from .geosetup import GeoSetup
from .elevationprofile import ElevationProfile
from .lineongrid import LineOnGrid

from .utils import delete_local_srtm_files

if BASEMAP_AVAILABLE: # pragma: no cover
    from .mapping import Map

