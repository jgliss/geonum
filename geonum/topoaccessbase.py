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
Access and handling of topographic data
"""
import numpy as np
import os
from warnings import warn
# python 2 and 3 support see e.g.
# http://python-future.org/compatible_idioms.html#metaclasses
from six import with_metaclass

import abc
 
from geonum import NETCDF_AVAILABLE, SRTM_AVAILABLE
from geonum.exceptions import (TopoAccessError, SRTMNotCoveredError, 
                               InvalidTopoMode)
   
class TopoAccessBase(with_metaclass(abc.ABCMeta, object)):
    """Abstract base class for topgraphy file implementations
    
    Defines minimum interface for derived access classes of different 
    topographic datasets. 
    """
    #local_path = None
    #topo_id = None
    
    #: A coordinate for which data should be available
    _TESTLAT = 45
    _TESTLON = 15
    
    def __init__(self, local_path=None, check_access=True):
        self.local_path = local_path
        self.topo_id = None
        
        if check_access:
            self.check_access()
    
    @abc.abstractmethod
    def get_data(self, lat0, lon0, lat1=None, lon1=None):
        """Declaration of data access method 
        
        It is obligatory to implement this method into derived classes.
        
        Parameters
        ----------
        lat0 : float
            first latitude coordinate of topographic range (lower left coord)
        lon0 : float
            first longitude coordinate of topographic range (lower left coord)
        lat1 : int or float, optional
            second latitude coordinate of topographic range (upper right 
            coord). If None only data around lon0, lat0 will be extracted.
        lon1 : int or float, optional
            second longitude coordinate of topographic range (upper right 
            coord). If None only data around lon0, lat0 will be extracted.
        
        Returns
        -------
        TopoData
            instance of TopoData class
        """
        pass
    
    def check_access(self):
        """Check if topography data can be accessed"""
        from geonum.topodata import TopoData
        try: 
            d = self.get_data(self._TESTLAT, self._TESTLON)
            if not isinstance(d, TopoData):
                raise ValueError('Invalid return type, expected instance '
                                 'of TopoData class, got {}'.format(type(d)))
            return True
        except Exception as e:
            print('Could not access topodata: {}'.format(repr(e)))
        return False
            
    def _prep_borders(self, lat0, lon0, lat1, lon1):
        """Sort by longitudes and determines LL and TR coordinates
        
        Parameters
        ----------
        lat0 : float
            first latitude coordinate of topographic range (lower left coord)
        lon0 : float
            first longitude coordinate of topographic range (lower left coord)
        lat1 : float
            second latitude coordinate of topographic range (upper right coord)
        lon1 : float
            second longitude coordinate of topographic range (upper right coord)
        
        Returns
        -------
        tuple
            4-element tuple, containing:
                
            - float, smallest latitude (LL corner) 
            - float, smallest longitude (LL corner)
            - float, largest latitude (TR corner)
            - float, largest longitude (TR corner)
        """
        lats, lons = np.asarray([lat0, lat1]), np.asarray([lon0, lon1])
        return (np.nanmin(lats), np.nanmin(lons), 
                np.nanmax(lats), np.nanmax(lons))
        
    def _init_lons_lats(self, lats_all, lons_all, lat0, lon0, lat1=None, 
                        lon1=None):
        """Get all latitudes and longitudes on a topodata grid 
        
        Parameters
        ----------
        lats_all : ndarray 
            numpy array containing available latitudes of the accessed topo 
            dataset
        lons_all : ndarray
            numpy array containing available longitudes of the accessed topo 
            dataset
        lat0 : float
            first latitude coordinate of topographic range (lower left coord)
        lon0 : float
            first longitude coordinate of topographic range (lower left coord)
        lat1 : float, optional
            second latitude coordinate of topographic range (upper right coord)
        lon1 : float, optional
            second longitude coordinate of topographic range (upper right coord)

        Returns
        -------
        tuple
            2-element tuple, containing
            
            - ndarray, topodata latitudes overlapping with input range
            - ndarray, topodata longitudes overlapping with input range
        """
        if any([x is None for x in [lat1, lon1]]):
            lat1, lon1 = lat0, lon0
        if lon0 > lon1:
            lon0, lon1 = lon1, lon0
            lat0, lat1 = lat1, lat0
        #print lat0, lon0, lat1, lon1     
        #closest indices
        idx_lons = [np.argmin(abs(lons_all - lon0)), 
                    np.argmin(abs(lons_all - lon1))]
        idx_lats = [np.argmin(abs(lats_all - lat0)), 
                    np.argmin(abs(lats_all - lat1))]
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
       
class Etopo1Access(TopoAccessBase):
    """A class representing netCDF4 data access of Etopo1 data
    
    See `here <https://github.com/jgliss/geonum#supported-etopo1-files>`_ for
    instructions on the data access.
    
    Attributes
    ----------
    loader 
        data loader (:class:`netCDF4.Dataset`)
    local_path : str
        directory where Etopo1 data files are stored
    file_name :  str
        file name of etopo data file
        
    Parameters
    ----------    
    local_path : str 
        directory where Etopo1 data files are stored
    file_name :  str
        file name of etopo data file
    check_access : bool
        if True, then access to topography data is checked on init and an
        error is raised if no dataset can be accessed
    search_database : bool
        if True and topodata file :attr:`file_path` does not exist, then
        a valid topography file is searched in all paths that are specified
        in file `LOCAL_TOPO_PATHS.txt` that is shipped with this library and 
        can be found in installation subdirectory local_topo_data.
        
    Raises
    ------
    TopoAccessError
        if input arg `check_access` is True and if no supported data file 
        can be found
    
    """
    #: ID of dataset
    topo_id = "etopo1"
    
    #: filenames of supported topographic datasets in preferred order
    supported_topo_files = ["ETOPO1_Ice_g_gmt4.grd",
                            "ETOPO1_Bed_g_gmt4.grd"]
    def __init__(self, local_path=None, file_name=None, check_access=False,
                 search_database=True):
        
        if not NETCDF_AVAILABLE:
            raise ImportError("Etopo1Access class cannot be initiated. "
                                      "Please install netCDF4 library first")
        self._local_path = None
        self._file_name = None
        
        from netCDF4 import Dataset
        
        self.loader = Dataset
        
        self.local_path = local_path
        self.file_name = file_name
        
        if not os.path.exists(self.file_path) and search_database:
            self.search_topo_file_database()
        
        # check if file exists
        if check_access:
            if not os.path.exists(self.file_path):
                raise TopoAccessError('File {} could not be found in local '
                                      'topo directory: {}'.format(self.file_name, 
                                                       self.local_path))
            elif not self.check_access():
                raise TopoAccessError('Failed to extract topography data for '
                                      'Etopo dataset')
    
    @property
    def local_path(self):
        """Directory containing ETOPO1 gridded data files"""
        return self._local_path
    
    @local_path.setter
    def local_path(self, val):
        if val is None or not os.path.exists(val):
            from geonum import LOCAL_TOPO_PATH
            print('Invalid input for local_path, setting default: '
                  '{}'.format(LOCAL_TOPO_PATH))
            val = LOCAL_TOPO_PATH
        self._check_topo_path(val)
        self._local_path = val
    
    @property
    def file_name(self):
        """File name of topographic dataset used"""
        return self._file_name
    
    @file_name.setter
    def file_name(self, val):
        if not val in self.supported_topo_files:
            
            _val = self.supported_topo_files[0]
            print('Invalid input for file_name ({})), setting default: '
                  '{}'.format(val, _val))
            val = _val
        self._file_name = val
        
    @property
    def file_path(self):
        """Return full file path of current topography file"""
        return os.path.join(self.local_path, self.file_name)
    
    @file_path.setter
    def file_path(self, val):
        self.set_file_location(val)
    
    def _check_topo_path(self, path):
        """Check if path exists and if it is already included in database
        
        Parameters
        ----------
        path : str
            path to be checked
        """
        from geonum.helpers import check_and_add_topodir
        check_and_add_topodir(path)
            
    def _get_all_local_topo_paths(self):
        """Get all search paths for topography files"""
        from geonum.helpers import all_topodata_search_dirs
        return all_topodata_search_dirs()
        
    def _search_topo_file(self, path=None):
        """Checks if a valid etopo data file can be found in local folder
        
        Searches in ``self.local_path`` for any of the file names specified 
        in ``supported_topo_files``
        
        """
        if path is None:
            path = self.local_path
        print(("Searching valid topo file in folder: %s" %path))
        fnames = os.listdir(path)
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
        files = os.listdir(self.local_path)
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
        
        Parameters
        ----------
        full_path : str
            full file path of topography file
            
        Raises
        ------
        TopoAccessError
            if filepath does not exist or if the provided file is not 
            supported by this interface.
        """
        if not os.path.exists(full_path):
            raise TopoAccessError('Input file location %s does not exist'
                                  .format(full_path))
        _dir = os.path.dirname(full_path)
        _f = os.path.basename(full_path)
        if not _f in self.supported_topo_files:
            raise TopoAccessError('Invalid topography data file name, please '
                                  'use either of the supported files from the '
                                  'Etopo1 data set: {}'
                                  .format(self.supported_topo_files))
        self.local_path = _dir
        self.file_name = _f
        
        if not os.path.basename(full_path) in self.supported_topo_files:
            raise TopoAccessError("Invalid topography data file, please use "
                "one of the supported files from the Etopo1 data set\n%s" 
                %self.supported_topo_files)
        self.local_path = os.path.dirname(full_path)
        self.file_name = os.path.basename(full_path)
    
    def get_data(self, lat0, lon0, lat1=None, lon1=None):
        """Retrieve data from topography file 
        
        Parameters
        ----------
        lat0 : float
            first latitude coordinate of topographic range (lower left coord)
        lon0 : float
            first longitude coordinate of topographic range (lower left coord)
        lat1 : int or float, optional
            second latitude coordinate of topographic range (upper right 
            coord). If None only data around lon0, lat0 will be extracted.
        lon1 : int or float, optional
            second longitude coordinate of topographic range (upper right 
            coord). If None only data around lon0, lat0 will be extracted.
        
        Returns
        -------
        TopoData
            instance of TopoData class
        """
        from geonum import TopoData
        etopo1 = self.loader(self.file_path)

        lons = etopo1.variables["x"][:]
        lats = etopo1.variables["y"][:]
        
        lats, lons, idx_lats, idx_lons = self._init_lons_lats(lats, lons, lat0,
                                                              lon0, lat1, lon1)

        vals = np.asarray(etopo1.variables["z"][idx_lats[0] : idx_lats[1] + 1,
                                                idx_lons[0] : idx_lons[1] + 1], 
                                                dtype = float)
        etopo1.close()

        return TopoData(lats, lons, vals, data_id=self.topo_id)
         
class SRTMAccess(TopoAccessBase):
    """Class for SRTM topographic data access
    
    Uses library `srtm.py <https://pypi.python.org/pypi/SRTM.py/0.3.1>`_ 
    for online access of data. 
    
    Note
    ----
    :mod:`srtm.py` downloads the topo data from `this source <http://
    dds.cr.usgs.gov/srtm/version2_1/>`_ and stores a copy of the unzipped data 
    files in the current cache directory found in home.
    
    Whenever data access is requested, the :mod:`srtm.py` checks if the file
    already exists on the local machine and if not downloads it online. The
    online access is rather slow, so do not be surprised, if things take a
    while when accessing a specific location for the first time. 
    
    **Deleting cached SRTM files**: 
        use :func:`geonum.topoaccess.delete_all_local_srtm_files`
      
    Parameters
    ----------
    check_access : bool
        check if data can be accessed on class initialisation
    **kwargs
        additional keyword arguments that are passed through (irrelevant for 
        this class but relevant for factory loader class 
        :class:`TopoDataAccess`, particularly :func:`set_mode` therein.
    """
    
    def __init__(self, check_access=False, **kwargs):
        """Class initialisation"""
        if not SRTM_AVAILABLE:
            raise ImportError("SRTM access class cannot be initiated. "
                                      "Please install srtm.py library first")
        import srtm
        self.loader = srtm
        self.topo_id = "srtm"
        
        if check_access:
            self.check_access()
        
    def _coordinate_covered(self, access_obj, lat, lon):
        """Checks if SRTM data is available for input coordinate
        
        Parameters
        ----------
        access_obj : GeoElevationData
            data access object from :mod:`srtm.py` module 
            (can be created calling ``srtm.get_data()``)
        lat : float
            latitude of point
        lon : float  
            longitude of point
        
        Returns
        -------
        bool
            True, if SRTM data is available for coordinate, else False.
        """
        if access_obj.get_file_name(lat, lon) is None:
            return False
        return True

        
    def get_data(self, lat0, lon0, lat1=None, lon1=None):
        """Load SRTM topographic subset for input range
        
        Parameters
        ----------
        lat0 : float
            first latitude coordinate of topographic range (lower left coord)
        lon0 : float
            first longitude coordinate of topographic range (lower left coord)
        lat1 : int or float, optional
            second latitude coordinate of topographic range (upper right 
            coord). If None only data around lon0, lat0 will be extracted.
        lon1 : int or float, optional
            second longitude coordinate of topographic range (upper right 
            coord). If None only data around lon0, lat0 will be extracted.
            
        Returns
        -------
        TopoData
            instance of TopoData class
        """
        from geonum import TopoData
        print("Retrieving SRTM data (this might take a while) ... ")
        # create GeoElevationData object for data access
        dat = self.loader.get_data()
        # check if second input point is specified and set equal first point if
        # not
        if any([x is None for x in [lat1, lon1]]):
            lat1, lon1 = lat0, lon0
        # Check if first point is covered by dataset
        if not self._coordinate_covered(dat, lat0, lon0):
            raise SRTMNotCoveredError('Point (lat={:.2f}, lon={:.2f}) not '
                                      'covered by SRTM'.format(lat0, lon0))
        # check if second point is covered by dataset
        if not self._coordinate_covered(dat, lat1, lon1):
            raise SRTMNotCoveredError('Endpoint coordinate (lat={:.2f}, '
                                      'lon={:.2f}) not covered by SRTM'
                                      .format(lat1, lon1))
        # prepare borders of covered lon / lat regime
        lat_ll, lon_ll, lat_tr,lon_tr = self._prep_borders(lat0, lon0,
                                                           lat1, lon1)
        # get SRTM file for lower left corner of regime
        f_ll = dat.get_file(lat_ll, lon_ll)
        # get SRTM file for top right corner of regime
        f_tr = dat.get_file(lat_tr, lon_tr)
        # create array of longitude values for regime
        lons_all = np.linspace(f_ll.longitude, f_tr.longitude + 1,
                            f_ll.square_side)
        # create array of latitude values for regime
        lats_all = np.linspace(f_ll.latitude, f_tr.latitude + 1,
                            f_ll.square_side)
        #prepare coordinates
        lats, lons, _, _= self._init_lons_lats(lats_all, lons_all,
                                               lat0, lon0, lat1, lon1)
        # Init data array 
        vals = np.ones((len(lats), len(lons))) * np.nan
        #loop over all coordinates and try access the elevation data
        for i in range(len(lats)):
            for j in range(len(lons)):
                #print "Lat: %s, Lon: %s" % (lats[i], lons[j])
                vals[i, j] = dat.get_elevation(lats[i], lons[j])

        return TopoData(lats, lons, vals, data_id=self.topo_id)

def delete_all_local_srtm_files():
    """Deletes all locally stored SRTM files"""
    if not SRTM_AVAILABLE:
        return
    import glob
    from srtm.main import FileHandler
    fh = FileHandler()
    for file in glob.glob('{}/*.hgt'.format(fh.get_srtm_dir())):
        print('Deleting SRTM data file at {}'.format(file))
        os.remove(file)
        
### Deprecated stuff
class TopoFile(TopoAccessBase):
    """Old name of :class:`TopoAccessBase`."""

    def __init__(self, *args, **kwargs):
        super(TopoFile, self).__init__(*args, **kwargs)
        warn(DeprecationWarning('Class TopoFile was renamed to TopoAccessBase'
                                'Please use the latter name'))