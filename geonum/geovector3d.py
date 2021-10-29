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
import warnings

from LatLon23 import GeoVector
from numpy import (radians, cos, sin, degrees, sqrt, tan, isnan, arctan2)

from geonum.geopoint import GeoPoint


class GeoVector3D(GeoVector):
    """A 3-dimensional geo vector object
    
    3D vector representation for geodesic connections and arithmetics of
    locations (:class:`geonum.GeoPoint`).

    Note
    ----
    To create an instnace of :class:`GeoVector3D` use one of the following to
    input combinations:

        1. dx, dy, dz
        #. dz, azimuth, dist_hor
        #. azimuth, dist_hor, elevation

    Multiple input combinations will be processed in the preferred order
    as given in the list.


    Note
    ----
    Horitontal displacement components dx (longitude) and dy (latitude) are
    in units of km, while vertical displacement component dz (altitude) is
    in units of m.

    Attributes
    ----------
    dx : float
        longitudinal component in units of km
    dy : float
        latitudinal component in units of km
    dz : float
        altitude component in units of m
    name : str
        name of vector


    Parameters
    ----------
    dx : float, optional
        longitudinal component in units of km
    dy : float, optional
        latitudinal component in units of km
    dz : float, optional
        altitude component in units of m
    azimuth : float, optional
        azimuth orientation angle in units of decimal degrees, relative to
        North direction.
    dist_hor : float, optional
        horizontal displacement length in input `azimuth` direction.
    elevation : float, optional
        elevation angle in decimal degrees (0 corresponds to horizon,
        90 to zenith)
    anchor : GeoPoint, optional
        anchor point of this vector.
    name : str, optional
        name of vector, defaults to "undefined".

    Example
    -------
    
    >>> from geonum import GeoPoint, GeoVector3D
    >>> p = GeoPoint(latitude=10.0, longitude=15.0, name="random_point")
    >>> v = GeoVector3D(dx=15, dy=100, dz=300) #dx, dy in km, dz in m
    >>> new_point = p + v # GeoPoint object
    >>> print(new_point)
    GeoPoint undefined
    Lat: 10.904041412793307,  Lon: 15.137202747069097, Alt: 300 m

    """
            
    def __init__(self, dx=None, dy=None, dz=None, azimuth=None, dist_hor=None, 
                 elevation=None, anchor=None, name=None):
        if name is None:
            name = 'undefined'

        self.name = name

        self.dx = None
        self.dy = None
        self.dz = None

        self._eval_input(dx,dy,dz,azimuth,dist_hor,elevation)

        # dictionary with private attributes
        self._priv_attr= {"anchor" : None}

        # call setter for private attribute anchor (this ensures that input 
        # attr anchor is of right type
        if anchor is not None:
            self.set_anchor(anchor)

    def _eval_input(self,dx,dy,dz,azimuth,dist_hor,elevation):
        if any(x is None for x in [dx, dy]): # If only initial_heading and
            # distance are given
            theta_rad = radians(self._angle_or_heading(azimuth))
            self.dx = dist_hor * cos(theta_rad)
            self.dy = dist_hor * sin(theta_rad)
        elif azimuth is None and dist_hor is None: # If only dx and dy are given
            self.dx = dx
            self.dy = dy
        else:
            raise ValueError('invalid input')

        # Check input for altitude component dz
        if dz is None or isnan(dz): #invalid for dz directly
            if elevation is not None and -90 <= elevation <= 90: #check if instead elevation is valid, then set dz
                #tan elev = dz/dist_hor
                dz = tan(radians(elevation))*sqrt(self.dx**2+self.dy**2)*1000
            else: #both dz input and elevation are invalid, set dz=0
                dz = 0.0
        self.dz = dz

    @property
    def dz_km(self) -> float:
        """:attr:`dz` converted to units of km

        E.g. used for internal arithmetics.
        """
        return self.dz/1000

    @property
    def azimuth(self) -> float:
        """Horizontal orientation angle relative to North direction"""
        return degrees(arctan2(self.dx, self.dy))

    @property
    def elevation(self) -> float:
        """Elevation angle in decimal degrees relative to horizon"""
        return (90-self.polar_angle)
    
    @property
    def polar_angle(self) -> float:
        """Polar angle in decimal degrees relative to zenith"""
        return degrees(arctan2(self.dist_hor, self.dz_km))
    
    @property    
    def magnitude(self) -> float:
        """Magnitude of vector (length) in units of km"""
        return sqrt(self.dist_hor**2 + self.dz_km**2)

    @property
    def norm(self) -> float:
        """Norm of vector, wrapper for :meth:`magnitude`"""
        return self.magnitude
        
    @property
    def dist_hor(self) -> float:
        """Horizontal distance spanned by this vector"""
        return sqrt(self.dx**2 + self.dy**2)
    
    @property
    def anchor(self) -> GeoPoint:
        """Anchor point of vector"""
        return self._priv_attr["anchor"]
        
    @anchor.setter
    def anchor(self, value):
        self.set_anchor(value)
    
    def set_anchor(self, geo_point):
        """Set anchor of this vector

        Parameters
        ----------
        geo_point : GeoPoint
            anchor location.

        Raises
        ------
        TypeError
            if input point is not instance of :class:`GeoPoint`.
        """
        if not isinstance(geo_point, GeoPoint):
            raise TypeError("Could not set anchor: Invalid input type")
        self._priv_attr["anchor"] = geo_point
           
    def intersect_hor(self, other) -> 'GeoVector3D':
        """Find horizontal intersection of this vector with another vector

        Note
        ----
        Only works if anchor point (:attr:`anchor`) is set in both vectors.

        Parameters
        ----------
        other : GeoVector3D
            Other vector for which the intersection is to be determined.

        Raises
        ------
        AttributeError
            if :attr:`anchor` is not set in this vector or input vector

        Returns
        -------
        GeoVector3D
            New vector pointing in direction of this vector but with correct
            horizontal magnitude to the intersection with the other vector.

        """
        if not all([isinstance(x,GeoPoint) for x in [self.anchor, other.anchor]]):
            raise AttributeError(
                "Intersection can not be determined, anchor "
                "of one of the vectors not set..")
        v = other.anchor - self.anchor
        v.dz = 0.0
        other_az = radians((other.azimuth + 360) % 360)
        self_az = radians((self.azimuth + 360) % 360)
        dy = (v.dx - tan(other_az) * v.dy) / (tan(self_az) - tan(other_az))
        dx = tan(self_az) * dy
        return GeoVector3D(dx, dy, dz=0.0)

    def plot(self, map, add_anchor=False, **kwargs): # pragma: no cover
        """Plot this vector into existing basemap"""
        warnings.warn(DeprecationWarning(
            'See https://github.com/jgliss/geonum/issues/4'))
        map.draw_geo_vector_3d(self, **kwargs)
        if add_anchor:
            self.anchor.plot(map, add_name=True, dz_text=self.dz*.1)

    # ToDo: check if __call__ is really needed, a little confusing...
    def __call__(self):
        """Call function

        Returns
        -------
        float
            azimuth angle in units of decimal degrees wrt to N direction
        float
            horizontal length in km
        float
            vertical length in m
        """
        return (self.azimuth, self.dist_hor, self.dz)
        
    def __neg__(self):
        """Returns negative of this vector"""
        return GeoVector3D(-self.dx, -self.dy, -self.dz)
    
    def __add__(self, other):
        """Add another geo vector"""
        return GeoVector3D(self.dx + other.dx, self.dy + other.dy,
                           self.dz + other.dz)

    def __str__(self) -> str:
        """String representation"""
        return ('GeoVector3D {}\n'
                'Azimuth: {:.2f}°, Elevation: {:.4f}°, Magnitude: {:.2f} km '
                '(hor: {:.2f} km)'
                .format(self.name, self.azimuth, self.elevation, 
                        self.magnitude, self.dist_hor))
            
    
    def __repr__(self) -> str:
        return ('Az %s, Elev %s, Mag. %s' 
                %(self.azimuth, self.elevation, self.magnitude))
    
    def type(self) -> str:
        """Returns object type"""
        return 'GeoVector3D'
