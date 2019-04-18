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

from LatLon23 import LatLon, GeoVector
    
from numpy import (radians, cos, sin, degrees, sqrt,tan, isnan, arctan2,
                   asarray, nanmin, nanmax)

from warnings import warn

class GeoPoint(LatLon):
    """The Geopoint object represents a location in the atmosphere 
    
    This object is in 3D and includes elevation information.
    
    Attributes
    ----------
    altitude : float
        altitude in m
    altitude_err : float
        uncertainty in altitude
    local_topo_path : str
        directory that is used to search for local topography if topodata is
        requested. Irrelevant for default topo
    Parameters
    ----------
    lat : float 
        latitude of point (decimal degrees)
    lon : float 
        longitude of point (decimal degrees)
    altitude : float
        elevation of point in m
    name : str
        name (ID) of this point
    topo_access_mode : str
        string specifying the current access mode for topographic data (in v1, 
        choose between srtm or etopo1)
    topo_path : str
        specify path where etopo1 data files are stored
    topo_data : TopoData 
        existing topographic dataset that can be assigned to this point (can 
        save time, e.g. for altitude access)
    """
    _ALTERR_DEFAULT = 99999999999.9
    def __init__(self, lat=0.0, lon=0.0, altitude=None, name="n/d",
                 topo_access_mode="srtm", topo_path=None, topo_data=None,
                 auto_topo_access=True):
        
        LatLon.__init__(self, float(lat), float(lon), name)
        self.altitude = altitude #altitude in m
        self.altitude_err = 0.0
        
        self.topo_access_mode = topo_access_mode
        self.local_topo_path = topo_path
        
        self.topo_data = None 
        
        if topo_data is not None:
            self.set_topo_data(topo_data)
        
        if self.altitude is None or isnan(self.altitude):
            if auto_topo_access:
                self.get_altitude()
            else:
                self.altitude = 0.0
                self.altitude_err = self._ALTERR_DEFAULT
    
    @property
    def longitude(self):
        """Longitude coordinate in decimal degrees"""
        return self.lon.decimal_degree

    @property
    def latitude(self):
        """Latitude coordinate in decimal degrees"""
        return self.lat.decimal_degree
    
    def offset(self, azimuth, dist_hor, dist_vert=0.0, ellipse="WGS84", 
               **kwargs):
        """Returns new GeoPoint at offset position
        
        Parameters
        -----------
        azimuth : float
            azimuth direction angle (in decimal degrees, e.g. 90 for E 
            direction, -90 or 270 for west)
        dist_hor : float 
            horizontal offset to this point in km 
        dist_vert : float 
            vertical offset to this point in m (positive if point is higher, 
            negative, if it is lower)
        ellipse : float 
            geodetic system (ellipsoid), default is "WGS84", i.e. 
            `World Geodetic System <https://confluence.qps.nl/pages/view page.
            action?pageId=29855173>`__
        **kwargs : 
            additional keyword arguments passed to init of new 
            :class:`GeoPoint` object
        
        Returns
        -------
        GeoPoint
            new location at offset position
        """
        p = LatLon(self.lat.decimal_degree, self.lon.decimal_degree)
        p1 = p.offset(azimuth, dist_hor, ellipse)
        return GeoPoint(p1.lat.decimal_degree, p1.lon.decimal_degree,
                        self.altitude + dist_vert, **kwargs)
    
    def set_topo_data(self, topo_data):
        """Assign topographic data to this point
        
        If input is valid, the topographic data set will be stored in 
        class attribute topo_data.
        
        Parameters
        ----------
        topo_data : TopoData
            topography dataset
        """
        from geonum.topodata import TopoData
        if not isinstance(topo_data, TopoData):
            return
        elif not topo_data.includes_coordinate(self.latitude, self.longitude):
            raise TypeError("Could not assign topodata to GeoPoint %s: "
             "topodata does not cover the coordinates" %self.name)
        self.topo_data = topo_data
    
    def range_borders(self, *points):
        """Get geographical borders (lower left, top right) of this point
        
        Additional points can be included as non keyword arguments. The
        range borders are then determined under consideration of all
        specified points. If no additional points are specified,
        the range corners will be set corresponding to 1km diagonal 
        distance to this point.
        
        Parameters
        ----------
        *points 
            additional :class:`GeoPoint` objects to be  considered
        
        Returns
        -------
        2-element tuple, containing:
            
            - :class:`GeoPoint`, lower left corner of regime
            - :class:`GeoPoint`, top right corner of regime
        
        """
        lats, lons = [self.latitude], [self.longitude]
        #retrieve all latitudes and longitudes
        for p in points:
            if isinstance(p, GeoPoint):
                lats.append(p.latitude)
                lons.append(p.longitude)
        lats, lons = asarray(lats), asarray(lons)
        if not len(lats) > 0:
            #print "Borders could not be initiated, no objects found..."
            return False
        lat_ll, lon_ll, lat_tr , lon_tr = (nanmin(lats), nanmin(lons), 
                                           nanmax(lats), nanmax(lons))
        pll, ptr = GeoPoint(lat_ll, lon_ll, 0.0), GeoPoint(lat_tr, lon_tr, 0.0)
        extend_km = (pll - ptr).magnitude * 0.1
        if extend_km == 0:
            extend_km = 1
        ll = pll.offset(azimuth=-135, dist_hor=float(extend_km), name="ll")
        tr = ptr.offset(azimuth=45, dist_hor=float(extend_km), name="tr")
        return  ll, tr        
        
    def get_topo_data(self, geo_point=None, azimuth=None, dist_hor=None,
                      lon1=None, lat1=None):
        """Retrieve topographic data grid spanned by this and another point
        
        The second coordinate can be specified in several different ways 
        dependent on which of the input arguments are used. First checks if the
        currently loaded topo data (``self.topo_data``) covers the range 
        spanned by the 2 points. If this is the case, the loaded data will be 
        used and nothing reloaded. If not, data will attempted to be loaded 
        from the current :class:`TopoDataAccess` object within the lon lat 
        range covered by the 2 points (set using :func:`range_borders`)
        
        Note
        ----
        The following input combinations work (and are preferentially 
        processed in the specified list order if multiple input is given):
        
            1. specify endpoint using geo_point
            #. specify endpoint using azimuth and dist_hor
            #. specify endpoint using lon1 and lat1
        
        Parameters
        ----------
        geo_point : GeoPoint 
            another geo_point,topo data will be retrieved between this point 
            and destination point
        azimuth : float 
            azimuth angle of heading in decimal degrees (to be used with 
            dist_hor)
        dist_hor : float 
            horizontal distance from this point in km (to be used with azimuth)
        lon1 : float
            longitude of destination point, topo data will be retrieved between 
            this point and destination point (to be used with lat1)
        lat1 : float 
            latitude of destination point, topo data will be retrieved between 
            this point and destination point (to be used with lon1)
        
        Returns
        -------
        2-element tuple, containing:
            
            - :class:`TopoData`, the topographic data
            - :class:`GeoPoint`, the second coordinate
        """
        pf = self
        if isinstance(geo_point, GeoPoint):
            pf = geo_point
        elif all(isinstance(x, (int, float)) for x in [dist_hor, azimuth]):
            pf = self.offset(azimuth, dist_hor)
        elif all(isinstance(x, (int, float)) for x in [lon1, lat1]):
            pf = GeoPoint(lat1, lon1)
        
        ll, tr = self.range_borders(self, pf)
        lat0, lon0 = ll.latitude, ll.longitude
        lat1, lon1 = tr.latitude, tr.longitude
        if not (self.check_topo(lat1, lon1) and self.check_topo(lat0, lon0)):
            print(("Topo data not avaialable in geo point %s, loading topodata"
                        %self.name))
            self._load_topo_data(lat0, lon0, lat1, lon1)
        else:
            print ("Topo data avaialable in geo point %s")
        return self.topo_data, pf
    
    def check_topo(self, lat1=None, lon1=None):
        """Check if topography is available between this point and another
        
        :param float lat1: latitude of end point, if None, only the 
            coordinates of this object will be considered
        :param float lon1: longitude of end point, if None, only the 
            coordinates of this object will be considered
        :returns: - bool, True if topo data is loaded, False if not
        """
        if any([x == None for x in [lat1, lon1]]):
            lat1, lon1 = self.latitude, self.longitude
        if self.topo_data is None:
            return False
        if not (self.topo_data.includes_coordinate(self.latitude, 
                                                   self.longitude) 
            and self.topo_data.includes_coordinate(lat1, lon1)):
            return False
        return True
        
    def get_elevation_profile(self, geo_point=None, azimuth=None, 
                              dist_hor=None, lon1=None, lat1=None, 
                              resolution=5., **mapping_opts):
        """Estimates the elevation profile for a given viewing direction up
        to a given distance from the coordinates of this point. For input 
        possibilities see docs of :func:`get_topo_data`.
        
        :param GeoPoint geo_point: another geo_point,topo data will be 
            retrieved between this point and destination point
        :param float azimuth: azimuth angle of heading in decimal degrees 
            (to be used with dist_hor)
        :param float dist_hor: horizontal distance from this point in km 
            (to be used with azimuth)
        :param float lon1: longitude of destination point, topo data will 
            be retrieved between this point and destination point (to be 
            used with lat1)
        :param float lat1: latitude of destination point, topo data will be 
            retrieved between this point and destination point 
            (to be used with lon1)
        :param float resolution: desired topo grid resolution in m. 
            1D interpolation of the elevation profile is performed is 
            applicable.
        :returns: 
            - :class:`ElevationProfile`, the profile object
          
        The following input combinations work (and are preferentially processed
        in the specified list order if multiple input is given):
        
            1. specify endpoint using geo_point
            #. specify endpoint using azimuth and dist_hor
            #. specify endpoint using lon1 and lat1
            
        
        """
        from geonum.elevationprofile import ElevationProfile
        data, pf = self.get_topo_data(geo_point, azimuth, dist_hor, lon1, 
                                      lat1)
        return ElevationProfile(data, self, pf, resolution=resolution,
                                **mapping_opts)
        
    def get_altitude(self):
        """Estimate the altitude (+/- uncertainty) of this point 
        
        The estimatation is done by retrieving a 2x2 grid of topography
        data enclosing this point. The altitude value is estimated using the 
        topo data tile closest to the coordinates of this point. The uncertainty
        is estimated conservatively using min/max difference of data in the 2x2
        grid.
        
        :returns: 
            - altitude (0 if retrieval failed)
            - altitude error (99999999999 if retrieval failes)
        :rtype:
            - int, float
        
        """
#==============================================================================
#         print ("Try retrieving altitude for point %s (lat, lon): %s, %s" 
#                                 %(self.name, self.latitude, self.longitude))
#==============================================================================
        try:
            #the data access class already catches access failure and in case,
            #tries to  access etopo1 data
            data = self._topo_access.get_data(self.latitude,
                                              self.longitude)
            z = data(self.latitude, self.longitude)# / 1000.
            z_err = (data.data.max() - data.data.min()) #/ 1000.'
#            print "Retrieved altitude: %s +/- %s m" %(z, z_err)
        except Exception as e:
            warn('Altitude of {} could not be retrieved...\n'
                 'Setting altitude to 0.0 m\n'
                 'Error: {}'.format(repr(self), repr(e)))
            z, z_err = 0.0, self._ALTERR_DEFAULT
    
        self.altitude, self.altitude_err = z, z_err
        
        return (z, z_err)
    
    """HELPERS, IO stuff, etc"""
    def update_topo_access(self, mode, local_path=None):
        """Update topo access mode
        
        :param str mode: valid access mode (e.g. "srtm", "etopo1")        
        :param str local_path (None): path for etopo1 access (if None, uses 
            the one which is currently set in ``self._topo_access``)
        """
        if local_path is not None:
            self.local_topo_path = local_path
        self.topo_access_mode = mode
        
        
    """Plotting etc..."""
    def plot_2d(self, map, add_name=False, dist_text=0.5, angle_text=-45, 
                **kwargs):
        """Plot this point into existing 2D basemap
                
        :param map basemap: Basemap object (drawn in an Axes3D object)
        :param bool add_name: add the name of this GeoPoint in the map
        :param float dist_text: distance of text annotation from point
        :param float angle_text: angle of text displacement
        :param kwargs: additional keyword arguments passed to matplotlib 
            plot function
        """
        if not "marker" in kwargs:
            kwargs["marker"] = "^"
        if not any([x in kwargs for x in ["c", "color"]]):
            kwargs["c"] = "lime"
        map.draw_geo_point_2d(self, **kwargs)
        if add_name:
            map.write_point_name_2d(self, dist_text, angle_text)
                                                    
    def plot_3d(self, map, add_name=False, dz_text=0.0, **kwargs):
        """Plot this point into existing 3D basemap
                
        :param map basemap: Basemap object (drawn in an Axes3D object)
        :param bool add_name: add the name of this GeoPoint in the map
        :param float dz_text: altitude offset of text(to point)
        :param kwargs: additional keyword arguments passed to matplotlib 
            plot function
        """
        if not "marker" in kwargs:
            kwargs["marker"] = "^"
        if not any([x in kwargs for x in ["c", "color"]]):
            kwargs["c"] = "lime"
        map.draw_geo_point_3d(self, **kwargs)
        x0, y0 = list(map(self.lon.decimal_degree, self.lat.decimal_degree))
        if add_name and self.name is not None:
            try:
                fs = kwargs["fontSize"]
            except:
                fs = 12
            #zt=self.altitude*1000. + dz_text
            zt = self.altitude + dz_text
            map.draw_text_3d(self.longitude, self.latitude, zt, self.name,
                             color="k", fontsize=fs)
            
    @property
    def _topo_access(self):
        """Topography data access class"""
        from geonum.topodataaccess import TopoDataAccess
        return TopoDataAccess(mode=self.topo_access_mode,
                              local_path=self.local_topo_path)
        
    def _load_topo_data(self, lat0, lon0, lat1, lon1):
        """Load and update topo data
        
        :param float lon0: start longitude
        :param float lat0: start latitude
        :param float lon1: stop longitude 
        :param float lat1: stop latitude
        :return:
            - :class:`TopoData` object with loaded data (it is also stored in
                ``self.topo_data``)
        """
        print("Trying to load topography data for {}".format(self))
        self.topo_data = self._topo_access.get_data(lat0, lon0, lat1, lon1)
        return self.topo_data
        
    def _sub_geo_vector_2d(self, other):
        """Subtract another geo vector
        
        Called when subtracting a GeoVector object from self (adapted from
        `LatLon` object,  only return type was changed from LatLon 
        object to GeoPoint object)
        
        :param (GeoVector, GeoVector3D) other: vector to be subtracted
        :return: :class:`GeoPoint` at difference position
        """
        heading, distance = other()
        heading = (heading + 180) % 360 # Flip heading
        p = GeoPoint(self.lat, self.lon, self.altitude)
        
        return p.offset(heading, distance)
    
    def _add_geo_vector_2d(self, other):
        """Subtract another geo vector
        
        Add a geo vector (adapted from `LatLon` object,  only return 
        type was changed from LatLon object to GeoPoint object)
        
        :param GeoVector other: Geovector added to this 
        
        """
        azimuth, dist_hor = other()
        p = GeoPoint(self.lat, self.lon, self.altitude)
        return p.offset(azimuth, dist_hor, 0.0)
        
    def _sub_geo_vector_3d(self, other):
        """Called when subtracting a GeoVector3D object from self"""
        azimuth, dist_hor, dz = other()
        azimuth = (azimuth + 180) % 360 # Flip heading
        p = GeoPoint(self.lat, self.lon, self.altitude) # Copy current position
        return p.offset(azimuth, dist_hor, -dz) # Offset position by GeoVector
    
    def _add_geo_vector_3d(self, other):
        """Add a 3d geo vector object
        
        :param GeoVector3D other: Geovector which is added to self
        """
        azimuth, dist_hor, dz = other()
        p = GeoPoint(self.lat, self.lon, self.altitude)
        return p.offset(azimuth, dist_hor, dz)
        
    def _sub_latlon(self, other):
        """Called when subtracting a LatLon object from this point"""
        from geonum import GeoVector3D
        inv = self._pyproj_inv(other)
        heading = inv['heading_reverse']
        distance = inv['distance']
        name=str(other.name + "->" + self.name)
        return GeoVector3D(dz=0.0, azimuth=heading, dist_hor=distance,
                           anchor=other, name=name)
    
    def _sub_geo_point(self, other):
        """Called when subtracting a LatLon object from self"""
        from geonum import GeoVector3D
        inv = self._pyproj_inv(other)
        heading = inv['heading_reverse']
        distance = inv['distance']
        name = str(other.name + "->" + self.name)
        dz = self.altitude - other.altitude
        return GeoVector3D(dz=dz, azimuth=heading, dist_hor=distance,
                           anchor=other, name=name)
        
    def __eq__(self, other):
        """Checks equality of this point with another
        
        :param other: another :class:`GeoPoint` object
        """
        if (self.lat == other.lat and self.lon == other.lon 
            and self.altitude == other.altitude):
            return True
        return False
        
    def __add__(self, other):
        """Add a geo vector to this point
        
        :param (GeoVector, GeoVector3D) other: a vector
        """
        object_operator = {'GeoVector'      :   self._add_geo_vector_2d,
                           'GeoVector3D'    :   self._add_geo_vector_3d}
        return object_operator[other.type()](other)

    def __sub__(self, other):
        """Subtraction"""
        object_operator = {'GeoVector'      :   self._sub_geo_vector_2d,
                           'GeoVector3D'    :   self._sub_geo_vector_3d,
                           'LatLon'         :   self._sub_latlon,
                           'GeoPoint'       :   self._sub_geo_point}
        return object_operator[other.type()](other)
    
    def __str__(self):
        """String formatting"""
        return ("GeoPoint %s\nLat: %s, Lon: %s, Alt: %s m\n" 
                %(self.name, self.latitude, self.longitude, self.altitude))
    
    def __repr__(self):
        """Obj. representation"""
        return ("%s, %s, Alt. %s m" 
                %(self.lat.__repr__(), self.lon.__repr__(), self.altitude))
    
    def __complex__(self):
        """Complex representation of lat and lon"""
        return self.complex()
    
    def type(self):
        """Object type identifier
        
        :returns: str, the string identifier        
        """
        return 'GeoPoint'
