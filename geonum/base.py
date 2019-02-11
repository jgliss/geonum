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
try:
    from LatLon23 import LatLon, GeoVector
except:
    from LatLon import LatLon, GeoVector
    
from numpy import (radians, cos, sin, degrees, sqrt,tan, isnan, arctan2,
                   asarray, nanmin, nanmax)

from warnings import warn
from geonum.topodata import TopoDataAccess, TopoData

class GeoPoint(LatLon):
    """The Geopoint object represents a location in the atmosphere 
    
    This object is in 3D and includes elevation information.
    
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
        
        self._topo_access = TopoDataAccess(mode=topo_access_mode,
                                           local_path=topo_path)
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
    def local_topo_path(self):
        """Returns current etopo1 data access path (str)"""
        return self._topo_access.local_path    
    
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
                              resolution=5.):
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
        from geonum.processing import ElevationProfile
        data, pf = self.get_topo_data(geo_point, azimuth, dist_hor, lon1, 
                                      lat1)
        return ElevationProfile(data, self, pf, resolution=resolution)
        
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
        except:
            warn("Altitude could not be retrieved...setting altitude to 0.0 m")
            z, z_err = 0.0, self._ALTERR_DEFAULT
    
        self.altitude, self.altitude_err = z, z_err
        
        return z, z_err
    
    """HELPERS, IO stuff, etc"""
    def update_topo_access(self, mode, local_path=None):
        """Update topo access mode
        
        :param str mode: valid access mode (e.g. "srtm", "etopo1")        
        :param str local_path (None): path for etopo1 access (if None, uses 
            the one which is currently set in ``self._topo_access``)
        """
        if local_path is None:
            local_path = self._topo_access.local_path
        self._topo_access = TopoDataAccess(mode=mode, 
                                           local_path=local_path)
        
    @property
    def topo_access_mode(self):
        """Return current topography access mode"""
        return self._topo_access.mode
        
    @property
    def longitude(self):
        """Returns longitude coordinate as decimal number"""
        return self.lon.decimal_degree

    @property
    def latitude(self):
        """Returns longitude coordinate as decimal number"""
        return self.lat.decimal_degree
        
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
            
    """PRIVATE METHODS"""  
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
        print("Trying to load topography data in GeoPoint obj..")
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
        inv = self._pyproj_inv(other)
        heading = inv['heading_reverse']
        distance = inv['distance']
        name=str(other.name + "->" + self.name)
        return GeoVector3D(dz=0.0, azimuth=heading, dist_hor=distance,
                           anchor=other, name=name)
    
    def _sub_geo_point(self, other):
        """Called when subtracting a LatLon object from self"""
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


class GeoVector3D(GeoVector):
    """A 3 dimensional geo vector object 
    
    3D vector representation for geo numerical calculations, intuitive 
    usage, i.e.::
    
        from geonum import GeoPoint, GeoVector3D
        p = GeoPoint(10.0, 15.0, name = "random_point") #lat, lon
        v = GeoVector3D(dx = 15, dy = 100, dz = 300) #dx, dy in km, dz in m
        
        new_point = p + v # GeoPoint object
        
    
    Note
    ----
    
        This class inherits and makes use of the functionality of 
        `GeoVector` objects of the `LatLon` library. Methods and 
        attributes are partly the same and partly overwritten
        (note that it is not initiated as :class:`GeoVector`)
 
    """
            
    def __init__(self, dx=None, dy=None, dz=None, azimuth=None, dist_hor=None, 
                 elevation=None, anchor=None, name="n/d"):
        """Class initialisation
        
        :param float dx: longitudinal length in km
        :param float dy: latitudinal length in km
        :param float dz: altitude length in m
        :param float azimuth: azimuth angle to North direction (corresponds to 
            heading_initial variable in :class:`GeoVector` object)
        :param float dist_hor: horizontal displacement of object in km
        :param float elevation: elevation angle in decimal degrees 
            (0 corresponds to horizon, 90, to zenith)
        
        For initiation use one of the following to input combinations:
        
            1. dx, dy, dz
            #. dz, azimuth, dist_hor
            #. azimuth, dist_hor, elevation
            
        Multiple input combinations will be processed in the preferred order
        as given in the list
        
        """
        self.name = name
        if any(x == None for x in [dx, dy]): # If only initial_heading and distance are given
            theta_rad = radians(self._angle_or_heading(azimuth))
            self.dx = dist_hor * cos(theta_rad)
            self.dy = dist_hor * sin(theta_rad)
        elif azimuth == None and dist_hor == None: # If only dx and dy are given
            self.dx = dx
            self.dy = dy
        else:
            raise NameError
        #Check input for altitude difference dz
        if dz is None or isnan(dz): #invalid for dz directly
            if elevation is not None and -90 <= elevation <= 90: #check if instead elevation is valid, then set dz
                #tan elev = dz/dist_hor
                dz = tan(radians(elevation)) * sqrt(self.dx**2 + self.dy**2)\
                                                                        * 1000
            else: #both dz input and elevation are invalid, set dz=0
                dz = 0.0
        self.dz = dz
        # dictionary with private attributes
        self._priv_attr= {"anchor" : None}
        # call setter for private attribute anchor (this ensures that input 
        # attr anchor is of right type
        try:
            self.set_anchor(anchor)
        except:
            pass
    
    @property
    def azimuth(self):
        """Azimuth angle in dec degrees 
        
        Horizontal orientation measured from N direction
        """
        return self._geom_hor()[0]

    @property
    def elevation(self):
        """Elevation angle in decimal degrees
        
        Measured from horizon (i.e. 0 -> horizon, 90 -> zenith)
        """
        return 90 - self.polar_angle    
    
    @property
    def polar_angle(self):
        """Polar angle in decimal degrees
        
        Returns the polar angle in decimal degrees (measured from zenith), 
        i.e.::
        
            90 - self.elevation
            
        """
        return degrees(arctan2(self.dist_hor, self.dz / 1000.))
    
    @property    
    def magnitude(self):
        """Return magnitude of vector (length)"""
        return sqrt(self.dist_hor**2 + (self.dz / 1000.)**2)

    @property
    def norm(self):
        """Return magnitude of vector (length)"""
        return self.magnitude
        
    @property
    def dist_hor(self):
        """Horizontal distance spanned by this vector"""
        return self._geom_hor()[1]
    
    @property
    def anchor(self):
        """Getter method for private attribute anchor"""
        return self._priv_attr["anchor"]
        
    @anchor.setter
    def anchor(self, value):
        """Setter method of private property anchor
        
        See :func:`set_anchor` for input / output specs
        """     
        self.set_anchor(value)
    
    def set_anchor(self, geo_point):
        """Set anchor of this vector
        
        :param GeoPoint geo_point: anchor point of this vector
        :raises TypeError: if input is not :class:`GeoPoint` object        
        """
        if not isinstance(geo_point, GeoPoint):
            raise TypeError("Could not set anchor: Invalid input type")
        self._priv_attr["anchor"] = geo_point
           
    def intersect_hor(self, other):
        """Determine the horizontal intersection of this vector with other 
        input vector
        
        :param GeoVector3D: other vector
        
        .. note::
        
            Only works if anchor is defined
        """
        if not all([isinstance(x,GeoPoint) for x in [self.anchor, other.anchor]]):
            raise ValueError("Intersection can not be determined, anchor "
                " of one of the vectors not set..")
        v = other.anchor - self.anchor
        v.dz = 0.0
        other_az = radians((other.azimuth + 360) % 360)
        self_az = radians((self.azimuth + 360) % 360)
        dy = (v.dx - tan(other_az) * v.dy) / (tan(self_az) - tan(other_az))
        dx = tan(self_az) * dy
        return GeoVector3D(dx, dy, dz=0.0)
        
    def _geom_hor(self):
        """Returns horizontal heading and horizontal magnitude"""
        return (degrees(arctan2(self.dx, self.dy)),
                sqrt(self.dx**2 + self.dy**2))
        
    """Plotting etc..."""
    def plot(self, map, add_anchor=False, **kwargs):
        """Plot this vector into existing basemap
        
        :param Map map: map object
        :param bool add_anchor: If true, the anchor point is plotted as well
        :param kwargs: additional keyword arguments
        
        .. note::
        
            1. The basemap needs to be set up with Axes3D
            #. ``self.anchor`` must be set with a :class:`Geopoint` instance
            
        """
        map.draw_geo_vector_3d(self, **kwargs)
        if add_anchor:
            self.anchor.plot(map, add_name=True, dz_text=self.dz*.1)
        
    """Redifining magic methods from base class :class:`GeoVector` object for
    3D calcs
    """
    def __call__(self):
        """Call function
        
        :returns:
            - float azimuth angle
            - float, horizontal length in km
            - float, vertical length in m
        """
        azimuth, dist_hor = self._geom_hor()
        return azimuth, dist_hor, self.dz
    
    
        
    def __neg__(self):
        """Returns negative of this vector"""
        return GeoVector3D(-self.dx, -self.dy, -self.dz)
    
    def __add__(self, other):
        """Add another geo vector
        
        
        :param GeoVector3D other: another geo vector  
        """
        return GeoVector3D(self.dx + other.dx, self.dy + other.dy,
                           self.dz + other.dz)

    def __str__(self):
        """String representation"""
        return ('GeoVector3D %s\nAzimuth: %s, Elevation: %s, Magnitude: %s m\n'
                %(self.name,self.azimuth, self.elevation, self.magnitude))
            
    
    def __repr__(self):
        return ('Az %s, Elev %s, Mag. %s' 
                %(self.azimuth, self.elevation, self.magnitude))
    
    def type(self):
        """Returns object type"""
        return 'GeoVector3D'
