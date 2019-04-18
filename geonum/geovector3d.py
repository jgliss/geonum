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

from LatLon23 import GeoVector
from geonum.geopoint import GeoPoint
    
from numpy import (radians, cos, sin, degrees, sqrt,tan, isnan, arctan2)

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
                dz = tan(radians(elevation))*sqrt(self.dx**2+self.dy**2)*1000
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
        return ('GeoVector3D {}\n'
                'Azimuth: {:.2f}°, Elevation: {:.4f}°, Magnitude: {:.2f} km '
                '(hor: {:.2f} km)'
                .format(self.name, self.azimuth, self.elevation, 
                        self.magnitude, self.dist_hor))
            
    
    def __repr__(self):
        return ('Az %s, Elev %s, Mag. %s' 
                %(self.azimuth, self.elevation, self.magnitude))
    
    def type(self):
        """Returns object type"""
        return 'GeoVector3D'
