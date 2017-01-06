from os.path import abspath, dirname
from pkg_resources import get_distribution

__version__ = get_distribution('geonum').version

_LIBDIR = abspath(dirname(__file__))

from .base import GeoPoint, GeoVector3D
from .geosetup import GeoSetup
from .topodata import TopoData, TopoDataAccess
from .processing import LineOnGrid, ElevationProfile  
from .mapping import Map
import helpers