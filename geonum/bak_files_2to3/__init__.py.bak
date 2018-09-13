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

from os.path import abspath, dirname, join
from pkg_resources import get_distribution

__dir__ = abspath(dirname(__file__))
__version__ = get_distribution('geonum').version

_LIBDIR = __dir__ #from older version

LOCAL_TOPO_PATH = join(_LIBDIR, "local_topo_data")

from .base import GeoPoint, GeoVector3D
from .geosetup import GeoSetup
from .topodata import TopoData, TopoDataAccess
from .processing import LineOnGrid, ElevationProfile  
from .mapping import Map
import helpers
import atmosphere