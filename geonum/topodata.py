# -*- coding: utf-8 -*-
#
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

"""
Access and handling of topography data

.. note::

    The SRTM dataset has no global coverage, it is therefore recommended
    to download one of the Etopo1 data set files as a backup if SRTM data
    cannot be accessed.
    
"""

from numpy import argmin, mod, ceil, log2, poly1d, polyfit,\
    arange, linspace, empty, nan, asarray, nanmax, nanmin, isnan

from os.path import basename, exists, join, dirname, normpath
from os import listdir
from warnings import warn
# python 2 and 3 support see e.g.
# http://python-future.org/compatible_idioms.html#metaclasses
from six import with_metaclass

from abc import ABCMeta, abstractmethod
 
from geonum import (NETCDF_AVAILABLE, CV2_AVAILABLE, LATLON_AVAILABLE, 
                    SRTM_AVAILABLE)
if LATLON_AVAILABLE:
    try:
        from LatLon23 import LatLon
    except:
        from LatLon import LatLon
        
class TopoFile(with_metaclass(ABCMeta, object)):
    """Abstract base class for topgraphy file implementations"""
    local_path = None
    topo_id = None
    
    @abstractmethod
    def get_data(self, lat0, lon0, lat1=None, lon1=None):
        """Declaration of data access method 
        It is obligatory to implement this method into inhereted classes
        """
        pass
    
    def _prep_borders(self, lat0, lon0, lat1, lon1):
        """Sort by longitudes and determines LL and TR coordinates
        
        :param float lat0: latitude of first point
        :param float lon0: longitude of first point
        :param float lat1: latitude of first point
        :param float lon1: longitude of first point
        :return: 
            - float, smallest latitude 
            - float, smallest longitude
            - float, largest latitude
            - float, largest longitude
        """
        lats, lons = asarray([lat0, lat1]), asarray([lon0, lon1])
        return nanmin(lats), nanmin(lons), nanmax(lats), nanmax(lons)
        
    def _init_lons_lats(self, lats_all, lons_all, lat0, lon0, lat1=None, 
                        lon1=None):
        """Get all latitudes and longitudes on a topodata grid 
        
        :param array lats_all: numpy array with all latitudes of the topo data 
            file        
        :param array lons_all: numpy array with all longitudes of the topo data 
            file
        :param float lat0: latitude of point of interest
        :param float lon0: longitude of point of interest
        :param float lat1: latitude of an (optional) second point
        :param float lon1: longitude of an (optional) second point
        :return: 
            - ndarray, latitudes
            - ndarray, longitudes
        """
        if any([x is None for x in [lat1, lon1]]):
            lat1, lon1 = lat0, lon0
        if lon0 > lon1:
            lon0, lon1 = lon1, lon0
            lat0, lat1 = lat1, lat0
        #print lat0, lon0, lat1, lon1     
        #closest indices
        idx_lons = [argmin(abs(lons_all - lon0)), argmin(abs(lons_all - lon1))]
        idx_lats = [argmin(abs(lats_all - lat0)), argmin(abs(lats_all - lat1))]
        #Make sure that the retrieved indices actually INCLUDE the input ranges
        if idx_lons[0] == 0 and lons_all[0] > lon0:
            warn("Error: Lon0 smaller than range covered by file, using first"
                " available index in topodata..")
            lon0 = lons_all[0]
            idx_lons[0] = 0
        elif lons_all[idx_lons[0]] > lon0:
            idx_lons[0] -= 1
        if idx_lons[1] == len(lons_all) - 1 and lons_all[-1] < lon1:
            warn("Error: Lon1 larger than range covered by file, using last"
                " available index in topodata..")
            lon1 = lons_all[-1]
            idx_lons[1] = len(lons_all) - 1
        elif lons_all[idx_lons[1]] < lon1:
            idx_lons[1] += 1
        if idx_lats[0] == 0 and lats_all[0] > lat0:
            warn("Error: Lat0 smaller than range covered by file, using first"
                " available index in topodata..")
            lat0 = lats_all[0]
            idx_lats[0] = 0
        elif lats_all[idx_lats[0]] > lat0:
            idx_lats[0] -= 1
        if idx_lats[1] == len(lats_all) - 1 and lats_all[-1] < lat1:
            warn("Error: Lat1 larger than range covered by file, using last"
                " available index in topodata..")
            lat1 = lats_all[-1]
            idx_lats[1] = len(lats_all) - 1
        elif lats_all[idx_lats[1]] < lat1:
            idx_lats[1] += 1
        #make sure that no odd array lengths occur
        if not (idx_lats[1] - idx_lats[0] + 1) %2 == 0:
            #try append index at the end
            if not idx_lats[1] == len(lats_all) - 1:
                idx_lats[1] += 1
            elif not idx_lats[0] == 0:
                idx_lats[0] -= 1
            else:
                raise ValueError("Fatal error, odd length of latitude array")
        if not (idx_lons[1] - idx_lons[0] + 1) %2 == 0:
            #try append index at the end
            if not idx_lons[1] == len(lons_all) - 1:
                idx_lons[1] += 1
            elif not idx_lons[0] == 0:
                idx_lons[0] -= 1
            else:
                raise ValueError("Fatal error, odd length of longitude array")
        if idx_lats[0] > idx_lats[1]:
            return (lats_all[idx_lats[1] : idx_lats[0] + 1], 
                    lons_all[idx_lons[0] : idx_lons[1] + 1],
                    idx_lats, idx_lons)
        else:
            return (lats_all[idx_lats[0] : idx_lats[1] + 1],
                    lons_all[idx_lons[0] : idx_lons[1] + 1],
                    idx_lats, idx_lons)
        
class Etopo1Access(TopoFile):
    """A class representing netCDF4 data access of Etopo1 data
    
    See `here <https://github.com/jgliss/geonum#supported-etopo1-files>`_ for
    instructions on the data access.
    
    """
    def __init__(self, local_path=None, file_name=None):
        """Class initialisation
        
        :param str local_path: directory where Etopo data files are stored
        :param str file_name: file name of etopo data file
        
        """
        if not NETCDF_AVAILABLE:
            raise ModuleNotFoundError("Etopo1Access class cannot be initiated. "
                                      "Please install netCDF4 library first")
        from netCDF4 import Dataset
        self.loader = Dataset
        self.topo_id = "etopo1"
        
        self.local_path = local_path
        self.file_name = file_name

        self.supported_topo_files = ["ETOPO1_Ice_g_gmt4.grd",
                                     "ETOPO1_Bed_g_gmt4.grd"]
        #checks input and if applicable adds to database
        self._check_topo_path(local_path) 
        if local_path is None or not exists(local_path):
            from geonum import LOCAL_TOPO_PATH
            self.local_path = LOCAL_TOPO_PATH
        
        # make sure the file is supported
        if not file_name in self.supported_topo_files:
            # this function only returns valid file names, else None
            self._search_topo_file() #searches current local_path for valid file
        if not exists(self.file_path):
            self.search_topo_file_database()
        
        # check if file exists
        if not exists(join(self.local_path, self.file_name)):
            raise TopoAccessError("File %s could not be found in local "
                "topo directory: %s" %(self.file_name, self.local_path))
    
    @property
    def file_path(self):
        """Return full file path of current topography file"""
        return join(self.local_path, self.file_name)
    
    def _check_topo_path(self, path):
        """Check if path exists and if it is already included in database"""
        if path is None or not exists(path):
            return 
        if not path in self._get_all_local_topo_paths():
            from geonum import LOCAL_TOPO_PATH
            with open(join(LOCAL_TOPO_PATH, "LOCAL_TOPO_PATHS.txt"), "a") as f:
                f.write("\n" + path  + "\n")
                print(("Adding new default local topo data path to "
                    "file LOCAL_TOPO_DATA.txt: %s" %path))
            f.close()
        
    def _get_all_local_topo_paths(self):
        """Get all search paths for topography files"""
        from geonum import LOCAL_TOPO_PATH
        paths = [LOCAL_TOPO_PATH]
        with open(join(LOCAL_TOPO_PATH, "LOCAL_TOPO_PATHS.txt"), "r") as f:
            lines = f.readlines()
            for line in lines:
                p = line.split("\n")[0]
                if exists(p):
                    paths.append(normpath(p))
        f.close()
        return paths
        
    def _check_local_topo_path(self, path):
        """Checks if local topo path exists and is part of database"""
        if path is None or not exists(path):
            return False
        
    def _search_topo_file(self, path=None):
        """Checks if a valid etopo data file can be found in local folder
        
        Searches in ``self.local_path`` for any of the file names specified 
        in ``supported_topo_files``
        
        """
        if path is None:
            path = self.local_path
        print(("Searching valid topo file in folder: %s" %path))
        fnames = listdir(path)
        for name in fnames:
            if name in self.supported_topo_files:
                self.file_name = name
                self.local_path = path
                print(("Found match, setting current filepath: %s" 
                %self.file_path))
                return True
        return False
        
    def _find_supported_files(self):
        """Look for all supported files in ``self.local_path```and return list"""
        files = listdir(self.local_path)
        lst = []
        for name in files:
            if name in self.supported_topo_files:
                lst.append(name)
        return lst
    
    def search_topo_file_database(self):
        """Checks if a valid topo file can be found in database"""
        all_paths = self._get_all_local_topo_paths()
        for path in all_paths:
            if self._search_topo_file(path):
                return True
        return False
        
    def set_file_location(self, full_path):
        """Set the full file path of a topography data file
        
        :param str full_path: file path with topo data
        """
        if not exists(full_path):
            raise TopoAccessError("Input file location %s does not exist"
                                                                %full_path)
        if not basename(full_path) in self.supported_topo_files:
            raise TopoAccessError("Invalid topography data file, please use "
                "one of the supported files from the Etopo1 data set\n%s" 
                %self.supported_topo_files)
        self.local_path = dirname(full_path)
        self.file_name = basename(full_path)
    
    def get_data(self, lat0, lon0, lat1=None, lon1=None):
        """Retrieve data from topography file 
        
        :param float lon0: start longitude for data extraction
        :param float lat0: start latitude for data extraction
        :param float lon1: stop longitude for data extraction (default: None). 
            If None only data around lon0, lat0 will be extracted
        :param float lat1: stop latitude for data extraction (default: None). 
            If None only data around lon0, lat0 will be extracted
        
        :returns:
            - ndarray (2D), map containing retrieved elevations
            - ndarray, array containing corresponding latitudes
            - ndarray, array containing corresponding longitudes
        
            
        """
#==============================================================================
#         print ("Trying topo data access Etopo1 within borders (Lon0 | Lat0), : "
#             "(Lon1 | Lat1): (%s | %s), (%s | %s)" %(lat0, lon0, lat1, lon1))
#==============================================================================
    
        etopo1 = self.loader(self.file_path)

        lons = etopo1.variables["x"][:]
        lats = etopo1.variables["y"][:]
        
        lats, lons, idx_lats, idx_lons = self._init_lons_lats(lats, lons, lat0,
                                                              lon0, lat1, lon1)

        vals = asarray(etopo1.variables["z"][idx_lats[0] : idx_lats[1] + 1,
                                             idx_lons[0] : idx_lons[1] + 1], 
                                             dtype = float)
        etopo1.close()

        return TopoData(lats, lons, vals, data_id=self.topo_id)
         
class SRTMAccess(TopoFile):
    """Class for SRTM topographic data access (uses 
    `SRTM <https://pypi.python.org/pypi/SRTM.py/0.3.1>`_ library for 
    access).
    
    .. note::
    
        :mod:`SRTM` downloads the topo data from 
        `here <https://www.ngdc.noaa.gov/mgg/global/global.html>`_
        `the official web page <http://dds.cr.usgs.gov/srtm/version2_1/>`_ and
        stores a copy of the unzipped data files in the current home directory, 
        i.e., either::
            import os
            os.path.join(os.environ["HOMEPATH"], ".cache\\srtm")
            
        or::
        
            os.path.join(os.environ["HOME"], ".cache\\srtm")
        
        Whenever data access is requested, the SRTM library checks if the file
        already exists on the local machine and if not downloads it online. The
        online access is rather slow, so do not be surprised, if things take a
        while when accessing a specific location for the first time.   
          
    """
    def __init__(self):
        """Class initialisation"""
        if not SRTM_AVAILABLE:
            raise ModuleNotFoundError("SRTM access class cannot be initiated. "
                                      "Please install srtm.py library first")
        import srtm
        self.loader = srtm
        self.topo_id = "srtm"
        
    def coordinate_covered(self, access_obj, lat, lon):
        """Checks if SRTM data is available for input coordinate
        
        :param GeoElevationData access_obj: data access object from 
            :mod:`srtm` module (can be created calling ``srtm.get_data()``)
        :param float lat: latitude of point
        :param float lon: longitude of point
        :returns: - bool
        """
        if access_obj.get_file_name(lat, lon) is None:
            return False
        return True

        
    def get_data(self, lat0, lon0, lat1=None, lon1=None):
        """Data access wrapper for SRTM module
        
        Retrieves SRTM data in the speciefied input range.
        
        :param float lon0: start longitude for data extraction
        :param float lat0: start latitude for data extraction
        :param float lon1: stop longitude for data extraction (default: None). 
            If None only data around lon0, lat0 will be extracted
        :param float lat1: stop latitude for data extraction (default: None). 
            If None only data around lon0, lat0 will be extracted
        
        :returns:
            - ndarray (2D), map containing retrieved elevations
            - ndarray, array containing corresponding latitudes
            - ndarray, array containing corresponding longitudes
        
        """
        print("Retrieving SRTM data (this might take a while) ... ")
        # create GeoElevationData object for data access
        dat = self.loader.get_data()
        # check if second input point is specified and set equal first point if
        # not
        if any([x is None for x in [lat1, lon1]]):
            lat1, lon1 = lat0, lon0
        # Check if first point is covered by dataset
        if not self.coordinate_covered(dat, lat0, lon0):
            raise SRTMNotCoveredError("Point not covered by SRTM (Lat|Lon): "
                                                    "(%s | %s)" %(lat0, lon0))
        # check if second point is covered by dataset
        if not self.coordinate_covered(dat, lat1, lon1):
            raise SRTMNotCoveredError("Destination point not covered by SRTM ("
                "Lat | Lon): (%s | %s)" %(lat0, lon0))
        # prepare borders of covered lon / lat regime
        lat_ll, lon_ll, lat_tr,lon_tr = self._prep_borders(lat0, lon0,
                                                           lat1, lon1)
        # get SRTM file for lower left corner of regime
        f_ll = dat.get_file(lat_ll, lon_ll)
        # get SRTM file for top right corner of regime
        f_tr = dat.get_file(lat_tr, lon_tr)
        # create array of longitude values for regime
        lons_all = linspace(f_ll.longitude, f_tr.longitude + 1,
                            f_ll.square_side)
        # create array of latitude values for regime
        lats_all = linspace(f_ll.latitude, f_tr.latitude + 1,
                            f_ll.square_side)
        #prepare coordinates
        lats, lons, _, _= self._init_lons_lats(lats_all, lons_all,
                                               lat0, lon0, lat1, lon1)
        # Init data array 
        vals = empty((len(lats), len(lons))) * nan
        #loop over all coordinates and try access the elevation data
        for i in range(len(lats)):
            for j in range(len(lons)):
                #print "Lat: %s, Lon: %s" % (lats[i], lons[j])
                vals[i, j] = dat.get_elevation(lats[i], lons[j])
    
        return TopoData(lats, lons, vals, data_id=self.topo_id)
        
class TopoDataAccess(object):
    """Factory class for accessing topographic data
    
    This is a high level access class which handles the access of topography
    data. It is, for instance, implemented within :class:`geonum.base.GeoPoint` 
    or :class:`geonum.geosetup.GeoSetup` classes. 
    
    Default access mode is SRTM. 
    """
    def __init__(self, mode="srtm", local_path=None):
        """Class initialisation
        
        :param str mode: choose between one of the supported data
            access types (currently 2: etopo1, srtm)
        :param str local_path: local path to etopo data (only relevant for 
            etopo1 mode)
        """
        #mode and data access variables
        self.mode = None
        self.modes = ["srtm", "etopo1"]
        self.local_path = None
        self.topo_file = SRTMAccess()
        
        self.set_mode(mode, local_path)
            
    def init_default_mode(self):
        """Set default mode (e.g. called when object is initiated with etopo1
        but without a file specification)
        """
        self.mode = "srtm"
        self.topo_file = SRTMAccess()
    
    def set_mode(self, mode="srtm", local_path=None):
        """Change the current topography mode 
        
        :param str mode: the new mode
        :param str local_path (None): option to update the current local_path 
            variable (i.e. the location of etopo1 grid files)
            
        :raises InvalidTopoMode: if input mode is not suppored
        :raises TopoAccessError: if mode "etopo1" and corresponding data 
            file cannot be found at ``self.local_path``
        """
        if local_path is not None: 
            if not exists(local_path):
                raise ValueError("Failed to set local topography path, path "
                                 "does not exist: %s" %local_path)
            self.local_path = local_path
        
        
        if not mode in self.modes:
            raise InvalidTopoMode("Mode %s not supported...\nSupported modes:"
                "%s\nCurrent mode: %s " %(mode, self.modes, self.mode))
        if mode == "etopo1":
            if self.local_path is None or not exists(self.local_path):
                raise IOError("Could not change to mode Etopo1, local path to "
                    "topography dataset not set, please download")
            tf = Etopo1Access(self.local_path)
            if tf.file_path is not None:
                self.mode = mode
                self.topo_file = tf
                return
            self.init_default_mode()
            raise TopoAccessError("Could not access Etopo1 data on local "
                                            "machine, check local_path...")
        
        self.init_default_mode() #which means SRTM
    
    def get_data(self, lat0, lon0, lat1=None, lon1=None):
        """Retrieve data from topography file and, if applicable, increase the 
        resolution by upscaling (Gaussian pyramide approach) and interpolation
        
        Data access wrapper for SRTM module
        
        Retrieves SRTM data in the speciefied input range.
        
        :param float lon0: start longitude for data extraction
        :param float lat0: start latitude for data extraction
        :param float lon1: stop longitude for data extraction (default: None). 
            If None only data around lon0, lat0 will be extracted
        :param float lat1: stop latitude for data extraction (default: None). 
            If None only data around lon0, lat0 will be extracted
                
        :returns: :class:`TopoData` object containing the data
        :raises: :class:`TopoAccessError`  if access fails
        
        """
        topo_data = None
        try:
            lat0, lon0 = float(lat0), float(lon0)
        except Exception as e:
            raise TopoAccessError(repr(e))
        try:
            topo_data = self.topo_file.get_data(lat0, lon0, lat1, lon1)
        except Exception as e:
            print(("Topo retrieval failed using mode " + str(self.mode) +
                                ", try extracting data from Etopo1 dataset"))
            msg = repr(e)
            try:
                self.set_mode("etopo1")
            except (InvalidTopoMode, TopoAccessError) as err:
                msg = "Error loading topodata: " + msg + "\n" + repr(err)
                raise TopoAccessError(msg)

            try:
                topo_data = self.topo_file.get_data(lat0, lon0, lat1, lon1)
                
            except Exception as err:
                msg = "Error loading topodata: %s\n%s" %(msg, repr(err))
                raise TopoAccessError(msg)
        try:
            print("Important remarks while accessing topodata: %s" %msg)
        except:
            pass
        return topo_data
        
class TopoData(object):
    """Class representing topography data as 2D numpy grid specified by corner
    coordinates (lon0, lat0, lon1, lat1), normally created by 
    :class:`TopoDataAccess`.
    
    .. note::
    
        This object does not provide topography data access. For data access, 
        please use :class:`TopoDataAccess` which returns :class:`TopoData` 
        objects. This object is intended for storage and post processing of
        topographic data 
        
    """
    def __init__(self, lats, lons, data, data_id="", repl_nan_minval=True):
        """Class initialisation
        
        :param ndarray lats: numpy array with latitude coordinates of the 
            topographic dataset
        :param ndarray lons: numpy array with longitude coordinates of the 
            topographic dataset
        :param ndarray data: 2D numpy array containing elevation values
        :param str data_id: string identification of this data set
        :param bool repl_nan_minval: pixels containing NaN values are 
            replaced with the minimum altitude in the range
        """
        self.data_id = data_id #: ID of topodata file

        self.lats = lats #asarray(lats)
        self.lons = lons #asarray(lons)
        if repl_nan_minval:
            data[isnan(data)] = nanmin(data)
        
        self.data = data
        
    @property
    def lon0(self):
        """Returns first longitude (i.e. ``self.lons[0]``)"""
        return self.lons[0]
        
    @property
    def lat0(self):
        """Returns first latitude (i.e. ``self.lats[0]``)"""
        return self.lats[0]
        
    @property
    def lon1(self):
        """Returns last longitude (i.e. ``self.lons[-1]``)"""
        return self.lons[-1]
        
    @property
    def lat1(self):
        """Returns last latitude (i.e. ``self.lats[-1]``)"""
        return self.lats[-1]
        
    @property
    def max(self):
        """Returns maximum altitude of topo data"""
        return nanmax(self.data)
        
    @property
    def min(self):
        """Returns minimum altitude of topo data"""
        return nanmin(self.data)
        
    @property
    def alt_range(self):
        """Returns covered altitude range"""
        return self.max - self.min
        
    @property
    def resolution(self):
        """Returns tuple (lat, lon) of the current grid resolution in km
        
        .. note::
        
            The resolution is determined at the center of this grid
        """
        if not LATLON_AVAILABLE:
            raise ModuleNotFoundError('Feature disabled: Neither LatLon nor '
                                      'LatLon23 are installed')
        x_lon, x_lat = int(len(self.lons) / 2), int(len(self.lats) / 2)
        p0 = LatLon(self.lats[x_lat], self.lons[x_lon])
        r_lon = (p0 - LatLon(self.lats[x_lat], self.lons[x_lon + 1])).magnitude
        r_lat = (p0 - LatLon(self.lats[x_lat + 1], self.lons[x_lon])).magnitude
        return (r_lat, r_lon)
        
    def increase_grid_resolution(self, res=0.2, polyorder=2):
        """Gaussian pyramide based upscaling of topographic grid
        
        This function checks the current topographic resolution in the center
        of the grid. Based on this, an upsacling factor is determined and the
        :func:`cv2.pyrUp` is used to increase the resolution using
        interpolation. Note, that this does not increase the actual resolution
        of the topographic data grid.
        
        :param float res: desired grid resolution in km (default: 0.2)
        :param int polyorder: order of polynomial used for interpolation
            (default: 2)
        :returns: 
            - :class:`TopoData`, new object with desired grid resolution
            
        .. note::
        
            This functionality is only available, if :mod:`cv2` is installed
            
            
        """
        if not CV2_AVAILABLE or not LATLON_AVAILABLE:
            raise ModuleNotFoundError('Feature disabled: Require opencv and '
                                      'LatLon (or LatLon23) library to change '
                                      'grid resolution ')
        from cv2 import pyrUp
        lons = self.lons
        lats = self.lats
        vals = self.data
        if not all(mod(x, 2) == 0 for x in [len(lons), len(lats)]):
            print ("Fatal error, odd array size detected, no upscaling possible."
                " Return current data dict")
            return False
            
        c_lon = len(lons) / 2 - 1 #center longitude index
        c_lat = len(lats) / 2 - 1 #center latitude index
        p1 = LatLon(lats[c_lat], lons[c_lon]) 
        p2 = LatLon(lats[c_lat + 1], lons[c_lon + 1])
        dist = (p2 - p1).magnitude #distance between 2 points on grid
        res_fac = dist / res #factor for increasing the spatial resolution
        if res_fac <= 1:
            print(("No interpolation of topodata necessary: topo raw data "
                "already has desired resolution...\nCurrent resolution: %s" 
                + str(dist) + "km\nDesired resolution: %s km" %(dist, res)))
            return self
        fac = int(ceil(log2(res_fac))) #the corresponding up-factor for the scale space
        print(("Increasing spatial topography resolution by factor %s" %(2**fac)))
        for k in range(fac):
            vals = pyrUp(vals)
        p_lons = poly1d(polyfit(arange(len(lons)), lons, polyorder))
        p_lats = poly1d(polyfit(arange(len(lats)), lats, polyorder))
        lons_new = p_lons(linspace(0, len(lons) - 1, vals.shape[1]))
        lats_new = p_lats(linspace(0, len(lats) - 1, vals.shape[0]))
        return TopoData(lats_new, lons_new, vals, self.data_id + "_interp")

    @property
    def delta_lon(self):
        """Longitude range of data (in decimal degrees)"""
        return self.lons[-1] - self.lons[0]
        
    @property 
    def delta_lat(self):
        """Latitude range of data (in decimal degrees)"""
        return self.lats[-1] - self.lats[0]
        
    @property
    def center_coordinates(self):
        """Tuple (lat, lon) with center coordinates of data"""
        return (self.lats[0] + self.delta_lat / 2.,
                self.lons[0] + self.delta_lon / 2.)
    
    def plot_2d(self, ax=None):
        """Plot 2D basemap of topodata"""
        from geonum.mapping import Map
        latc, lonc = self.center_coordinates
        m = Map(self.lon0, self.lat0, self.lon1, self.lat1,
                projection="merc", lat_0=latc, lon_0=lonc, ax=ax)
        m.set_topo_data(self)
        m.draw_topo(insert_colorbar=True)
        return m
        
    def plot_3d(self, ax=None):
        """Creates 3D surface plot of data
        
        :param Axes3D ax (None): 3D axes object
        :return:
            - :class:`geonum.mapping.Map` object
            
        """
        from geonum.mapping import Map
        latc, lonc = self.center_coordinates
        m = Map(self.lon0, self.lat0, self.lon1, self.lat1, 
                projection="merc", lat_0=latc, lon_0=lonc, ax=ax)
        m.set_topo_data(self)
        m.draw_topo_3d()
        return m
        
    def includes_coordinate(self, lat=None, lon=None):
        """Checks if input coordinates are covered by this dataset"""
        if not (self.lat0 <= lat <= self.lat1 and 
                self.lon0 <= lon <= self.lon1):
            print("Input coordinates out of topo data range...")
            return False
        return True
    
    def get_altitude(self, lat, lon):
        """Get altitude value at input position
        
        :param float lat: latitude of point
        :param float lon: longitude of point
        """
        if not self.includes_coordinate(lat, lon):
            raise ValueError("Input values out of range...")
        idx_lon = argmin(abs(self.lons - lon))
        idx_lat = argmin(abs(self.lats - lat))
        dat = self.data[idx_lat, idx_lon]
        if isnan(dat):
            raise ValueError("Invalid value encountered in topodata...")
        return dat
        
    def __call__(self, lat, lon):
        """Get altitude value at input position"""
        try:
            return self.get_altitude(lat, lon)
        except ValueError as e:
            raise e
        
class SRTMNotCoveredError(Exception):
    pass

class TopoAccessError(Exception):
    pass

class InvalidTopoMode(Exception):
    pass

def delete_all_local_srtm_files():
    """Deletes all locally stored SRTM files"""
    raise NotImplementedError
    #from srtm.main import FileHandler
    #fh = FileHandler()
#==============================================================================
#     dir = fh.get_srtm_dir()
#     filelist = [ f for f in os.listdir(dir) if f.endswith(".bak") ]
#==============================================================================
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    topopath = r'D:/Dropbox/Python27/jgliss/modules_priv/pygeomapping/topodata/'
    
    etopo_cell = "ETOPO1_Ice_c_gmt4.grd"
    etopo_grid = "ETOPO1_Ice_g_gmt4.grd"
    
    lon0, lat0 = 14.85, 37.7
    lon1, lat1 = 15.2, 37.85
    cell = Etopo1Access(topopath, etopo_cell)
    grid = Etopo1Access(topopath, etopo_grid)
    
    srtm_access = SRTMAccess()
    
    d1 = cell.get_data(lat0, lon0, lat1, lon1)
    d2 = grid.get_data(lat0, lon0, lat1, lon1)
    d3 = srtm_access.get_data(lat0, lon0, lat1, lon1)    
    
    fig, ax = plt.subplots(1,3,figsize=(18,6))
    ax[0].imshow(d1[0])
    ax[1].imshow(d2[0])
    ax[2].imshow(d3[0])
    
    access = TopoDataAccess()
    td = access.get_data(lat0, lon0, lat1, lon1)