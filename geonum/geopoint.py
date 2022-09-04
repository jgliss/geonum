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
from LatLon23 import LatLon

from geonum.exceptions import OutOfDomain
from geonum.topodata import TopoData
from geonum.topodataaccess import TopoDataAccess


class GeoPoint(LatLon):
    """The Geopoint object represents a location in the atmosphere

    This object is in 3D and includes altitude information.

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
    latitude : float, optional
        latitude of point (decimal degrees), defaults to 0.
    longitude : float, optional
        longitude of point (decimal degrees), defaults to 0.
    altitude : float, optional
        elevation of point in m, defaults to 0.
    name : str, optional
        name (ID) of this point, defaults to undefined
    topo_access_mode : str, optional
        string specifying the current access mode for topographic data (in v1,
        choose between srtm or etopo1), defaults to None, in which case default
        of :class:`TopoDataAccess` is used if topographic data is accessed via
        this point (e.g. if `auto_topo_access` is True). As of v1.4.X, default
        access mode is srtm, which uses SRTM database, which does not have
        full global coverage!
    topo_path : str, optional
        local path where local topography data is stored (at the moment only
        relevant for topo_mode etopo1), defaults to None.
    topo_data : TopoData, optional
        existing topographic dataset that can be assigned to this point (can
        save time, e.g. for altitude access)
    auto_topo_access : bool
        if True, try set :attr:`altitude`and :attr:`altitude_err` based on
        local topography (note that defa)
    """
    _ALTERR_DEFAULT = 1e9

    def __init__(self, latitude=0, longitude=0, altitude=None, name=None,
                 topo_access_mode=None, topo_path=None, topo_data=None,
                 auto_topo_access=False):
        if name is None:
            name = 'undefined'
        LatLon.__init__(self, float(latitude), float(longitude), name)
        if altitude is None:
            altitude = 0
        elif auto_topo_access:  # altitude is not None
            raise ValueError(
                'either provide altitude or deactivate auto_topo_access'
            )
        self.altitude = altitude  # altitude in m
        self.altitude_err = self._ALTERR_DEFAULT

        self.topo_access_mode = topo_access_mode
        self.local_topo_path = topo_path

        self.topo_data = None

        if topo_data is not None:
            self.set_topo_data(topo_data)

        if auto_topo_access:
            self.set_topo_altitude()

    @property
    def longitude(self) -> float:
        """Longitude coordinate in decimal degrees"""
        return self.lon.decimal_degree

    @property
    def latitude(self) -> float:
        """Latitude coordinate in decimal degrees"""
        return self.lat.decimal_degree

    @property
    def _topo_access(self) -> TopoDataAccess:
        """Topography data access class"""
        return TopoDataAccess(mode=self.topo_access_mode,
                              local_path=self.local_topo_path)

    def offset(self, azimuth, dist_hor, dist_vert=0.0, ellipse=None,
               **kwargs) -> 'GeoPoint':
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
        if ellipse is None:
            ellipse = "WGS84"
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
            raise ValueError('need instance of geonum.TopoData')
        elif not topo_data.includes_coordinate(self.latitude, self.longitude):
            raise OutOfDomain(
                f'Cannot assign input TopoData {topo_data} to {self}. '
                f'The location of {self} are outside the borders of TopoData'
            )
        self.topo_data = topo_data

    def range_borders(self, *points, extend_fac=0.1):
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
        extend_fac : float
            domain extension factor wrt to lat / lon span covered by input
            points, defaults to 0.1 (ca 10% extension)

        Returns
        -------
        GeoPoint
            lower left corner of regime
        GeoPoint
            top right corner of regime

        """
        lats, lons = [self.latitude], [self.longitude]
        # retrieve all latitudes and longitudes
        for p in points:
            if isinstance(p, GeoPoint):
                lats.append(p.latitude)
                lons.append(p.longitude)
        lats, lons = np.asarray(lats), np.asarray(lons)

        (lat_ll, lon_ll,
         lat_tr, lon_tr) = (np.nanmin(lats), np.nanmin(lons),
                            np.nanmax(lats), np.nanmax(lons))
        pll, ptr = GeoPoint(lat_ll, lon_ll, 0.0), GeoPoint(lat_tr, lon_tr, 0.0)
        extend_km = (pll - ptr).magnitude * extend_fac
        if extend_km == 0:
            extend_km = 1
        ll = pll.offset(azimuth=-135, dist_hor=float(extend_km), name="ll")
        tr = ptr.offset(azimuth=45, dist_hor=float(extend_km), name="tr")
        return (ll, tr)

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
        processed in the specified order if multiple input is given):

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
        TopoData
            topographic data object.
        GeoPoint
            the second coordinate which deterimes the domain wrt to this
            location.
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
            self._load_topo_data(lat0, lon0, lat1, lon1)
        return self.topo_data, pf

    def check_topo(self, lat1=None, lon1=None):
        """Check if topography is available between this point and another

        Parameters
        ----------
        lat1 : float
            latitude of end point, if None, only the coordinates of this
            object will be considered
        lon1 : float
            longitude of end point, if None, only the coordinates of this
            object will be considered

        Returns
        -------
        bool
            True if topo data is loaded, False if not
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
        """Estimates the elevation profile for a given viewing direction

        The profile is retrieved up to a given distance from the coordinate of
        this point. For input possibilities see docs of :func:`get_topo_data`.

        Note
        ----
        The following input combinations work (and are preferentially processed
        in the specified list order if multiple input is given):

            1. specify endpoint using geo_point
            #. specify endpoint using azimuth and dist_hor
            #. specify endpoint using lon1 and lat1

        Parameters
        ----------
        geo_point : GeoPoint
            end location of elevation profile (if this is provided, the
            following other input options to determine the end point will be
            ignored).
        azimuth : float
            azimuth angle (direction) of elevation profile in decimal degrees
            (to be used with input arg dist_hor)
        dist_hor : float
            horizontal distance from this point in km (to be used with azimuth)
        lon1 : float
            longitude of destination point, topo data will
            be retrieved between this point and destination point (to be
            used with lat1)
        lat1 : float
            latitude of destination point, topo data will be
            retrieved between this point and destination point
            (to be used with lon1)
        resolution : float
            desired topo grid resolution in m.
            1D interpolation of the elevation profile is performed if
            applicable.

        Returns
        -------
        ElevationProfile
            the profile object
        """
        from geonum.elevationprofile import ElevationProfile
        topo, pf = self.get_topo_data(geo_point, azimuth, dist_hor, lon1, lat1)
        if pf == self:
            raise ValueError('please specify endpoint location')
        return ElevationProfile(observer=self,
                                endpoint=pf,
                                topo_data=topo,
                                calc_on_init=True,
                                resolution=resolution,
                                **mapping_opts)

    def set_topo_altitude(self):
        """Set altitude using topographic terrain height at lat / lon position

        The estimation is done by retrieving a 2x2 grid of topography
        data enclosing this point. The altitude value is estimated using the
        topo data tile closest to the coordinates of this point. The
        uncertainty is estimated conservatively using min/max difference of
        data in the 2x2 grid.

        Returns
        -------
        float
            altitude
        float
            uncertainty in altitude value
        """
        if not isinstance(self.topo_data, TopoData):
            data = self._topo_access.get_data(self.latitude,
                                              self.longitude)
            # store for later usage
            self.set_topo_data(data)
        else:
            data = self.topo_data
        z = data(self.latitude, self.longitude)
        z_err = (data.data.max() - data.data.min())
        self.altitude, self.altitude_err = z, z_err
        return (z, z_err)

    def plot_2d(self, map, add_name=False, dist_text=0.5, angle_text=-45,
                **kwargs):  # pragma: no cover
        """Plot this point into existing 2D basemap

        Parameters
        ----------
        map : basemap
            Basemap object (drawn in an Axes3D object)
        add_name : bool
            add the name of this GeoPoint in the map
        dist_text : float
            distance of text annotation from point
        angle_text : float
            angle of text displacement
        **kwargs
            additional keyword arguments passed to matplotlib plot function
        """
        if not "marker" in kwargs:
            kwargs["marker"] = "^"
        if not any([x in kwargs for x in ["c", "color"]]):
            kwargs["c"] = "lime"
        map.draw_geo_point_2d(self, **kwargs)
        if add_name:
            map.write_point_name_2d(self, dist_text, angle_text)

    def plot_3d(self, map, add_name=False, dz_text=0.0,
                **kwargs):  # pragma: no cover
        """Plot this point into existing 3D basemap

        map : basemap
            Basemap object (drawn in an Axes3D object)
        add_name : bool
            add the name of this GeoPoint in the map
        dz_text : float
            altitude offset of text (to point location in map)
        **kwargs
            additional keyword arguments passed to matplotlib plot function
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
            # zt=self.altitude*1000. + dz_text
            zt = self.altitude + dz_text
            map.draw_text_3d(self.longitude, self.latitude, zt, self.name,
                             color="k", fontsize=fs)

    def _load_topo_data(self, lat0, lon0, lat1, lon1) -> TopoData:
        """Load topo data

        Parameters
        ----------
        lon0 : float
            start longitude
        lat0 : float
            start latitude
        lon1 : float
            stop longitude
        lat1 : float
            stop latitude

        Return
        ------
        TopoData
            loaded topographic data
        """
        self.topo_data = self._topo_access.get_data(lat0, lon0, lat1, lon1)
        return self.topo_data

    def _sub_geo_vector_2d(self, other):
        """Subtract another geo vector

        Called when subtracting a GeoVector object from self (adapted from
        `LatLon` object,  only return type was changed from LatLon
        object to GeoPoint object)

        Parameters
        ----------
        other : GeoVector
            vector to be subtracted from this location

        Returns
        -------
        GeoPoint
            endpoint location
        """
        heading, distance = other()
        heading = (heading + 180) % 360  # Flip heading
        p = GeoPoint(self.lat, self.lon, self.altitude)

        return p.offset(heading, distance)

    def _add_geo_vector_2d(self, other):
        """Add another geo vector

        Add a geo vector (adapted from `LatLon` object,  only return
        type was changed from LatLon object to GeoPoint object)

        Parameters
        ----------
        other : GeoVector
            vector to be added to this location

        Returns
        -------
        GeoPoint
            endpoint location

        """
        azimuth, dist_hor = other()
        p = GeoPoint(self.lat, self.lon, self.altitude)
        return p.offset(azimuth, dist_hor, 0.0)

    def _sub_geo_vector_3d(self, other):
        """Subtract a GeoVector3D object from self

        Parameters
        ----------
        other : GeoVector3D
            vector to be subtracted from this location

        Returns
        -------
        GeoPoint
            endpoint location
        """
        azimuth, dist_hor, dz = other()
        azimuth = (azimuth + 180) % 360  # Flip heading
        p = GeoPoint(self.lat, self.lon, self.altitude)  # Copy current position
        return p.offset(azimuth, dist_hor, -dz)  # Offset position by GeoVector

    def _add_geo_vector_3d(self, other):
        """Add a 3d geo vector object

        Parameters
        ----------
        other : GeoVector3D
            vector to be added to this location

        Returns
        -------
        GeoPoint
            endpoint location
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
        on = other.name
        if on is None:
            on = 'undefined'
        name = f'{on}->{self.name}'
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

        Parameters
        ----------
        other : GeoPoint
            other location

        Returns
        -------
        bool
            True if other is equal, else False
        """
        if (self.latitude == other.latitude and
                self.longitude == other.longitude and
                self.altitude == other.altitude):
            return True
        return False

    def __add__(self, other):
        """Add a geo vector or geo point to this point"""
        object_operator = {'GeoVector': self._add_geo_vector_2d,
                           'GeoVector3D': self._add_geo_vector_3d}
        try:
            tp = other.type()
            if not tp in object_operator.keys():
                raise AttributeError
        except AttributeError:
            raise ValueError(f'invalid input type {type(other)}, choose from '
                             f'{list(object_operator)}')
        return object_operator[other.type()](other)

    def __sub__(self, other):
        """Subtraction"""
        object_operator = {'GeoVector': self._sub_geo_vector_2d,
                           'GeoVector3D': self._sub_geo_vector_3d,
                           'LatLon': self._sub_latlon,
                           'GeoPoint': self._sub_geo_point}
        try:
            tp = other.type()
            if not tp in object_operator.keys():
                raise AttributeError
        except AttributeError:
            raise ValueError(f'invalid input type {type(other)}, choose from '
                             f'{list(object_operator)}')
        return object_operator[other.type()](other)

    def __str__(self):
        """String formatting"""
        return (
            f"GeoPoint {self.name}\nLat: {self.latitude},  Lon:"
            f" {self.longitude}, Alt: {self.altitude} m\n")

    def __repr__(self):
        """Obj. representation"""
        return (
            f"Lat: {self.latitude}, Lon: {self.longitude}, Alt:"
            f" {self.altitude} m")

    def type(self):
        """Object type identifier

        Returns
        -------
        str
            value='GeoPoint'
        """
        return 'GeoPoint'

    @staticmethod
    def from_LatLon(coord):
        if not isinstance(coord, LatLon):
            raise ValueError('Need LatLon...')
        return GeoPoint(latitude=coord.lon.decimal_degree,
                        longitude=coord.lat.decimal_degree,
                        altitude=0,
                        auto_topo_access=False)
