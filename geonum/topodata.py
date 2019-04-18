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

import numpy as np
from geonum import (CV2_AVAILABLE, LATLON_AVAILABLE)
if LATLON_AVAILABLE:
    from LatLon23 import LatLon
   
class TopoData(object):
    """Data class for topography data
    
    This object represents topographic data in a certain latitude and longitude 
    range, specified by the corner coordinates of the range 
    (:attr:`lon0, lat0, lon1, lat1`). It may be used for 2D and 3D plotting of 
    the topographic data (cf. :func:`plot_2d`, :func:`plot_3d`) or for further
    analysis such as coordinate based altitude retrievals 
    (:func:`get_altitude`), change of the grid resolution (cf. 
    :func:`increase_grid_resolution`) or for customised analysis by directly
    using the numpy data array containing the topographic altitudes 
    (:attr:`data`), together with the corresponding dimension coordinates 
    which can be accessed via attributes :attr:`latitude`and :attr:`longitude`.
    
    To create an instance of this class, you may use the class 
    :class:`TopoDataAccess`.
    
    Note
    ----
    This object does not provide topography data access. For data access, 
    please use :class:`TopoDataAccess` which returns :class:`TopoData` 
    objects. This object is intended for storage and post processing of
    topographic data 
        
    Parameters
    ----------
    lats : ndarray 
        numpy array with latitude coordinates of the topographic dataset.
        Accessible via :attr:`latitude`.
    lons : ndarray 
        numpy array with longitude coordinates of the topographic dataset.
        Accessible via :attr:`longitude`.
    data : ndarray
        2D numpy array containing elevation values.
    data_id : str
        ID of this data set.
    repl_nan_minval : bool 
        coordinates containing NaN values are replaced with the minimum 
        altitude in the range.
    """
    def __init__(self, lats, lons, data, data_id="", repl_nan_minval=False):
        self.data_id = data_id #: ID of topodata file

        self.lats = lats #
        self.lons = lons #asarray(lons)
        
        if repl_nan_minval:
            self.replace_nans()
        
        self.data = data
        
    @property
    def latitude(self):
        """Wrapper for :attr:`lats`"""
        return self.lats
    
    @property
    def longitude(self):
        """Wrapper for :attr:`lats`"""
        return self.lats
    
    def replace_nans(self, fillval=None):
        """Replace NaNs in topographic data with a fill value
        
        Parameters
        ----------
        fillval : float
            value that is assigned to all coordinates that include NaNs
        """
        if fillval is None:
            fillval = np.nanmin(self.data)
        self.data[np.isnan(self.data)] = fillval
        
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
        return np.nanmax(self.data)
        
    @property
    def min(self):
        """Returns minimum altitude of topo data"""
        return np.nanmin(self.data)
    
    @property
    def shape(self):
        """Shape of 2D numpy array that contains topography"""
        return self.data.shape
    
    @property
    def alt_range(self):
        """Returns covered altitude range"""
        return self.max - self.min
        
    @property
    def resolution(self):
        """Returns tuple (lat, lon) of the current grid resolution in km
        
        Note
        ----
        The resolution is determined at the center of this grid
        """
        if not LATLON_AVAILABLE:
            raise ImportError('Feature disabled: Neither LatLon nor '
                                      'LatLon23 are installed')
        x_lon, x_lat = int(len(self.lons) / 2), int(len(self.lats) / 2)
        p0 = LatLon(self.lats[x_lat], self.lons[x_lon])
        r_lon = (p0 - LatLon(self.lats[x_lat], self.lons[x_lon + 1])).magnitude
        r_lat = (p0 - LatLon(self.lats[x_lat + 1], self.lons[x_lon])).magnitude
        return (r_lat, r_lon)
    
    def mean(self):
        """Mean elevation (ignores NaNs in data)"""
        return np.nanmean(self.data)  
    
    def std(self):
        """Standard deviation of topographic dataset"""
        return np.nanstd(self.data)
        
    def increase_grid_resolution(self, res=0.2, polyorder=2):
        """Gaussian pyramide based upscaling of topographic grid
        
        This function checks the current topographic resolution in the center
        of the grid. Based on this, an upsacling factor is determined and the
        :func:`cv2.pyrUp` is used to increase the resolution using
        interpolation. Note, that this does not increase the actual resolution
        of the topographic data grid.
        
        Note
        ----
        Requires that opencv library is installed.
        
        Parameters
        ----------
        res : int or float, optional
            desired grid resolution in km (default: 0.2)
        polyorder : int, optional 
            order of polynomial used for interpolation (default: 2)
        
        Returns
        -------
        TopoData
            new object with desired grid resolution
        
        Raises
        ------
        ImportError
            if opencv is not installed.
        """
        if not CV2_AVAILABLE or not LATLON_AVAILABLE:
            raise ImportError('Feature disabled: Require opencv and '
                                      'LatLon (or LatLon23) library to change '
                                      'grid resolution ')
        from cv2 import pyrUp
        lons = self.lons
        lats = self.lats
        vals = self.data
        if not all(np.mod(x, 2) == 0 for x in [len(lons), len(lats)]):
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
        fac = int(np.ceil(np.log2(res_fac))) #the corresponding up-factor for the scale space
        print(("Increasing spatial topography resolution by factor %s" %(2**fac)))
        for k in range(fac):
            vals = pyrUp(vals)
        p_lons = np.poly1d(np.polyfit(np.arange(len(lons)), lons, polyorder))
        p_lats = np.poly1d(np.polyfit(np.arange(len(lats)), lats, polyorder))
        lons_new = p_lons(np.linspace(0, len(lons) - 1, vals.shape[1]))
        lats_new = p_lats(np.linspace(0, len(lats) - 1, vals.shape[0]))
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
    
    def plot(self, plot3d=True, draw_coastlines=False, 
             draw_mapscale=False, **kwargs):
        from geonum.mapping import Map
        if not "projection" in kwargs:
            kwargs["projection"] = "lcc"
        if not "llcrnrlon" in kwargs:
            kwargs["llcrnrlat"] = self.lat0
            kwargs["llcrnrlon"] = self.lon0
            kwargs["urcrnrlat"] = self.lat1
            kwargs["urcrnrlon"] = self.lon1
        
            kwargs["lat_0"] = self.lat0 + (self.lat1-self.lat0)/2
            kwargs["lon_0"] = self.lon0 + (self.lon1-self.lon0)/2
        
        m = Map(**kwargs)
        m.set_topo_data(self)
        
        if plot3d:
            if draw_coastlines:
                raise NotImplementedError('3D plotting and coastline drawing is '
                                          'not supported yet')
            elif draw_mapscale:
                raise NotImplementedError('3D plotting and mapscale drawing is '
                                          'not supported yet')
            m.draw_topo_3d()
        else:  
            if draw_coastlines:
                m.draw_coastlines()
            m.draw_coordinates()
            if draw_mapscale:
                m.draw_mapscale_auto()
        return m
    
    def plot_2d(self, ax=None):
        """Plot 2D basemap of topodata"""
        from geonum.mapping import Map
        latc, lonc = self.center_coordinates
        m = Map(self.lon0, self.lat0, self.lon1, self.lat1,
                projection="merc", lat_0=latc, lon_0=lonc, ax=ax)
        m.set_topo_data(self)
        m.draw_topo(insert_colorbar=True)
        m.draw_topo_contour() 
        m.draw_coordinates()
        m.drawcoastlines()
        return m
        
    def plot_3d(self, ax=None):
        """Creates 3D surface plot of data
        
        Parameters
        ----------
        ax
            instance of matplotlib Axes3D object
        
        Returns
        -------
        geonum.mapping.Map
            plotted map object.
            
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
    
    def closest_index(self, lat, lon):
        """Finds closest index to input coordinate
        
        Parameters
        ----------
        lat : float
            latitude coordinate in decimal degrees
        lon : float 
            longitude coordinate in decimal degrees
        
        Returns
        -------
        tuple
            2-element tuple containing closest index of lat and lon arrays to
            to input index
        
        Raises
        ------
        ValueError
            if input coordinate is not included in this dataset
        """
        if not self.includes_coordinate(lat, lon):
            raise ValueError("Input values out of range...")
        
        idx_lat = np.argmin(abs(self.lats - lat))
        idx_lon = np.argmin(abs(self.lons - lon))
        return (idx_lat, idx_lon)
        
    def get_altitude(self, lat, lon):
        """Get altitude value at input coordinate
        
        Parameters
        ----------
        lat : int or float 
            latitude of coordinate
        lon : int or float 
            longitude of coordinate
            
        Returns
        -------
        float
            altitude at input coordinate
        
        Raises
        ------
        ValueError
            if retrieved altitude value is NaN
        """
        idx_lat, idx_lon = self.closest_index(lat, lon)
        dat = self.data[idx_lat, idx_lon]
        if np.isnan(dat):
            raise ValueError("Invalid value encountered in topodata...")
        return dat
        
    def __call__(self, lat, lon):
        """Get altitude value at input position"""
        try:
            return self.get_altitude(lat, lon)
        except ValueError as e:
            raise e