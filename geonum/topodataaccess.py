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
import os

from geonum.exceptions import InvalidTopoMode
from geonum.topoaccessbase import SRTMAccess, Etopo1Access

class TopoDataAccess(object):
    """Factory class for accessing topographic data
    
    This is a high-level factory class which handles the access of topography
    data. Registered topographic datasets can be found in :attr:`REGISTERED` 
    
    Default access mode is SRTM (topographic dataset from NASA 
    `Shuttle Radar Topography Mission <https://www2.jpl.nasa.gov/srtm/>`__). 
    
    Example
    -------
    
    >>> acc = TopoDataAccess(mode='srtm')
    >>> topo_data = acc.get_data(lat0=5, lon0=30, lat1=15, lon1=40)
        
    Note
    ----
    
    For developers: the registered access  classes (such as 
    :class:`SRTMAccess`) should be based on (or follow the API of) template 
    base class :class:`TopoAccessBase`.
    
    
    Attributes
    ----------
    mode : str
        current access mode (string specifying which topographic dataset is
        supposed to be used for access). 
    local_path : str
        local path to etopo data (only relevant for etopo1 mode)
        
    Parameters
    ----------
    mode : str
        one of the supported data access types (cf. keys of :attr:`REGISTERED`)
    local_path : str
        local path to etopo data (only relevant for etopo1 mode)
    """
    #: supported access modes (topographic datasets)
    REGISTERED = dict(srtm      = SRTMAccess,
                      etopo1    = Etopo1Access)
    
    def __init__(self, mode="srtm", local_path=None):
        #mode and data access variables
        if not mode in self.REGISTERED:
            raise InvalidTopoMode(mode)
        self.mode = mode
    
        self.local_path = local_path
        
    @property
    def modes(self):
        """List of supported topographic datasets"""
        return list(self.REGISTERED.keys())
    
    @property
    def supported(self):
        """List of supported datasets (wrapper for :attr:`modes`)"""
        return self.modes
    
    def __deepcopy__(self, memo):
        return TopoDataAccess(self.mode, self.local_path)
    
    def get_data(self, lat0, lon0, lat1=None, lon1=None, 
                 mode=None, local_path=None, **access_opts):
        """Retrieve data from topography file
        
        Parameters
        ----------
        lat0 : float 
            start longitude for data extraction
        lon0 : float 
            start latitude for data extraction
        lat1 : float 
            stop longitude for data extraction (default: None). If None only 
            data around lon0, lat0 will be extracted.
        lon1 : float
            stop latitude for data extraction (default: None). 
            If None only data around lon0, lat0 will be extracted
        mode : str, optional
            mode specifying the topographic dataset that is supposed to be used
        local_path : str, optional
            local path where topography data is stored (can be dictionary or
            filepath, is passed to corresponding access class and handled 
            as implemented there)
        **access_opts
            additional access options that may be specific for the mode 
            specified (e.g. search_database in case of etopo1)
        
        
        Returns
        -------
        TopoData
            object containing the data
        
        Raises
        ------
        TopoAccessError
            if access fails
        
        """
        if mode is not None:
            if not mode in self.REGISTERED:
                raise InvalidTopoMode(mode)
            self.mode = mode
        if local_path is not None and os.path.exists(local_path):
            self.local_path=local_path
        access = self.REGISTERED[self.mode](local_path=self.local_path,
                                             **access_opts)
        return access.get_data(lat0, lon0, lat1, lon1)