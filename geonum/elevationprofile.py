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
from scipy.interpolate import interp1d

from geonum.exceptions import IntersectNotFound
from geonum.geopoint import GeoPoint
from geonum.processing import LineOnGrid
from geonum.topodata import TopoData

class ElevationProfile(object):
    """Class for calculating elevation profiles

    This class can be used to compute elevation profiles for arbitrary
    locations on earth. This is done by using topographic data which is
    provided (instance of :class:`geonum.TopoData`). Profiles are computed with respect to
    input location (:attr:`observer`) and the azimuthal direction of the
    profile is specified via the provided enpoint location (:attr:`endpoint`).
    The profile is calculated between two :class:`GeoPoint` objects using a
    provided topographic data grid which must cover the range spanned by
    both points.

    Parameters
    ----------
    observer : GeoPoint
        start point of profile
    endpoint : GeoPoint
        end point of profile
    topo_data : TopoData
        topographic dataset
    calc_on_init : bool, optional
        if True, then the profile is calculated on class inialisation.
        Default is True.
    **kwargs
        additional keyword args passed to :func:`det_profile` (only relevant
        if `calc_on_init` is True)
    """
    def __init__(self, observer, endpoint, topo_data=None,
                 calc_on_init=True, **kwargs):

        self._topo_data = None
        self._observer = None #: coordinate of observer (start of profile)
        self._endpoint = None #: coordinate of endpoint of profile

        self._check_set_observer(observer)
        self._check_set_endpoint(endpoint)

        if topo_data is not None:
            self.topo_data = topo_data

        self._init_attrs()

        if calc_on_init:
            self.det_profile(**kwargs)

    @property
    def observer(self):
        """
        GeoPoint: start coordinate of elevation profile
        """
        return self._observer

    @observer.setter
    def observer(self, value):
        self._check_set_observer(value)
        self._init_attrs()

    @property
    def endpoint(self):
        """
        GeoPoint: end coordinate of elevation profile
        """
        return self._endpoint

    @endpoint.setter
    def endpoint(self, value):
        self._check_set_endpoint(value)
        self._init_attrs()

    @property
    def topo_data(self):
        """
        TopoData: topographic data used to extract elevation profile
        """
        if self._topo_data is None:
            topo, _ = self.observer.get_topo_data(geo_point=self.endpoint)
            self._topo_data = topo
        return self._topo_data

    @topo_data.setter
    def topo_data(self, value):
        self._check_set_topo_data(value)
        self._init_attrs()

    @property
    def line(self):
        """
        LineOnGrid: line along which the profile is calculated

        Note
        ----
        private attribute that cannot be set by the user but is calculated
        automatically based on attributes :attr:`observer`, :attr:`endpoint`
        and :attr:`topo_data`.
        """
        if not isinstance(self._line, LineOnGrid):
            if self._coords_topo is None:
                self._set_latlon_indices_topogrid()
            self._line = LineOnGrid(**self._coords_topo)
        return self._line

    @property
    def observer_topogrid(self):
        """
        GeoPoint: Location of endpoint on topogrid (depends on topo resolution)

        Note
        ----
        private attribute that cannot be set by the user but is calculated
        automatically based on attributes :attr:`observer` and
        :attr:`topo_data`.
        """
        if self._observer_topogrid is None:
            if self._coords_topo is None:
                self._set_latlon_indices_topogrid()
            self._observer_topogrid = GeoPoint(
                self.topo_data.lats[self._coords_topo['y0']],
                self.topo_data.lons[self._coords_topo['x0']],
                topo_data=self.topo_data)
        return self._observer_topogrid

    @property
    def endpoint_topogrid(self):
        """
        GeoPoint: Location of endpoint on topogrid (depends on topo resolution)

        Note
        ----
        private attribute that cannot be set by the user but is calculated
        automatically based on attributes :attr:`endpoint` and
        :attr:`topo_data`.
        """
        if self._endpoint_topogrid is None:
            if self._coords_topo is None:
                self._set_latlon_indices_topogrid()
            self._endpoint_topogrid = GeoPoint(
                self.topo_data.lats[self._coords_topo['y1']],
                self.topo_data.lons[self._coords_topo['x1']],
                topo_data=self.topo_data)
        return self._endpoint_topogrid

    @property
    def dist_hor(self):
        """
        float: Horizontal distance in km between start and endpoint
        """
        return (self.endpoint_topogrid - self.observer_topogrid).dist_hor


    @property
    def azimuth(self):
        """
        float: Azimuth angle of profile (wrt to North direction)
        """
        return (self.endpoint_topogrid - self.observer_topogrid).azimuth

    @property
    def profile(self):
        """
        Retrived altitude levels in m along :attr:`line`

        Raises
        ------
        AttributeError
            if profile has not been calculated yet

        Returns
        -------
        ndarray
        """
        if self._profile is None:
            raise AttributeError(
                'Profile information not available, call det_profile first')
        return self._profile

    @property
    def dists(self):
        """
        Distances from observer in km along :attr:`line`

        Raises
        ------
        AttributeError
            if profile has not been calculated yet

        Returns
        -------
        ndarray
        """
        if self._dists is None:
            raise AttributeError(
                'Distance information not available, call det_profile first')
        return self._dists

    @property
    def resolution(self):
        """Get profile x (distance) resolution (averaged from distance array)

        Note
        ----
        Only works if profile was already determined

        Returns
        -------
        float
            average resolution along horizontal dimension

        """
        return abs((self.dists[1:] - self.dists[:-1]).mean())

    @property
    def gradient(self):
        """Return gradient of profile

        Uses numpy function ``gradient`` on :attr:`profile`

        Returns
        -------
        ndarray

        """
        return np.gradient(self.profile)

    @property
    def slope(self):
        """Returns slope of profile

        Determines dx and dy arrays and returns slope vector::

            slope = dy / dx = gradient(self.profile) /
                (gradient(self.dists) * 1000.0)

        """
        return np.gradient(self.profile)/(np.gradient(self.dists)*1000)

    @property
    def start_point(self):
        """Return position of observer"""
        return self.observer

    @property
    def min(self):
        """Minimum altitude in profile"""
        return np.nanmin(self.profile)

    @property
    def max(self):
        """Maximum altitude in profile"""
        return np.nanmax(self.profile)

    @property
    def alt_range(self):
        """Altitude range of profile (dh = max - min)
        """
        return self.max - self.min

    def slope_angles(self, decimal_degrees=True):
        """
        Get slope angles of profile (in each sample point)

        Parameters
        ----------
        decimal_degrees : bool
            if True, angles are converted from radians to degrees

        Returns
        -------
        ndarray
            slopes at each profile point
        """
        a = np.tan(self.slope)
        if decimal_degrees:
            a = np.rad2deg(a)
        return a

    def slope_angle(self, dist):
        """
        Get slope angle of profile at input distance from observer

        Parameters
        ----------
        dist : float
            distance from observer in km for which slope is to be retrieved

        Returns
        -------
        float
            retrieved slope at input distance.
        """
        if not dist >= 0 or dist > np.nanmax(self.dists):
            raise ValueError('invalid input: dist must be positive')
        idx = np.argmin(abs(self.dists - dist))
        return self.slope_angles()[idx]

    def det_profile(self, interpolate=True, resolution=5.0, itp_type="linear",
                    **mapping_opts):
        """
        Determine the elevation profile

        Searches the closest tiles in the topo data grid for both observer and
        endpoint and based on these two points the elevation profile is
        extracted using a :class:`LineOnGrid` object (which extracts altitudes
        from the topo data along the connection vector of the 2 points using
        2D spline intepolation)

        Parameters
        ----------
        interpolate : bool
            if True, the profile is interpolated to a certain horizontal
            resolution
        resolution : float
            desired grid resolution in m for interpolation. Interpolation is
            performed if :attr:`interpolate` is True and if the actual
            resolution of the topo data is smaller than input, else, nothing
            is done
        itp_type : str
            interpolation type (e.g. "linear", "cubic")
        **mapping_opts
            additional keyword arguments that are passed to
            :func:`LineOnGrid.get_line_profile`

        Returns
        -------
        ndarray
            the array containing the retrieved altitude levels along the
            retrieval direction (from the observer)
        """
        line = self.line
        topo = self.topo_data
        altitudes = line.get_line_profile(topo.data, **mapping_opts)

        dists = np.linspace(0, self.dist_hor, line.length + 1)
        if interpolate:
            dists, altitudes = self._interp_helper(dists, altitudes,
                                                   resolution, itp_type)
        self._dists = dists
        self._profile = altitudes
        return altitudes

    def _interp_helper(self, dists, altitudes, resolution, itp_type):
        """
        Helper to interpolate profile

        Parameters
        ----------
        dists : ndarray
            Array with distances from observer
        altitudes : ndarray
            Array with retrieved altitudes at input distances
        resolution : int or float
            desired horizontal resolution in m
        itp_type : str
            interpolation type (passed to :func:`scipy.interpolate.interp1d`)

        Returns
        -------
        dists : ndarray
            interpolated distances from observer.
        altitudes :  ndarray
            interpolated altitude levels.
        """
        res0 = (dists[1] - dists[0]) * 1000
        fac = int(np.ceil(res0 / resolution))
        if fac > 1:
            fz = interp1d(dists, altitudes, kind=itp_type)
            dists = np.linspace(0, self.dist_hor, self.line.length * fac)
            altitudes = fz(dists)
        else:
            print('No interpolation of elevation profile needed, data already '
                  'has sufficient resolution')
        return dists, altitudes

    def get_altitudes_view_dir(self, elev_angle, view_above_topo_m=0):
        """
        Get vector containing altitudes for a viewing direction

        The viewing direction is specified by the azimuth angle of the
        connection vector between observer and endpoint, the elevation angle
        needs to be specified via the input parameters, and the start point
        (first elevation value) corresponds to the altitude at the observer
        position plus an offset in m which can be specified.

        Parameters
        ----------
        elev_angle : float
            elevation angle of viewing direction
        view_above_topo_m : float
            altitude offset of start point in m, defaults to 0.

        Returns
        -------
        ndarray
            vector with altitude values (same length as ``self.profile``)

        """
        dh = 1000 * np.tan(np.radians(elev_angle))*self.dists
        result = dh+self.profile[0] + view_above_topo_m
        return result

    def get_first_intersection(self, elev_angle,view_above_topo_m=1.5,
                               min_dist=None,local_tolerance=3,
                               max_diff=None, plot=False):

        """Finds first intersection of a viewing direction with topography

        Start point of the viewing vector is the observer position (or more
        accurately, the center position of the closest topography tile in
        the topography data set). The relative altitude of the start point
        with respect to the topography altitude at the observer position
        can be controlled on input (arg ``view_above_topo_m``) as well as
        the elevation angle of the viewing direction. The azimuth angle
        corresponds to the profile azimuth (i.e. azimuth between
        :attr:`observer` and :attr:`endpoint` of this profile).

        The signal analysed for finding the intersection is a vector
        containing relative elevation positions of the viewing direction
        with respect to the profile altitude::

            diff_signal = self.get_altitudes_view_dir - self.profile

        The first intersection of the viewing direction with the topography
        is identified by the first zero crossing of this ``diff_signal``.

        Parameters
        ----------
        elev_angle : float
            elevation angle of viewing direction (in decimal degrees)
        view_above_topo_m : float
            altitude offset in m relative to altitude of start point (
            :attr:`observer`) in m. Defaults to 1.5 m.
        min_dist : float
            minimum distance (in km) of first considered intersection with
            topography from observer. If None, use 1% of distance between
            :attr:`observer` and :attr:`endpoint`.
        local_tolerance : int
            tolerance factor to estimate distance uncertainty based on topo
            grid resolution. Defaults to 3.
        max_diff : float, optional
            maximum allowed altitude difference in m between vector specifying
            viewing direction and the topography. Default to None, in which
            case the value assigned to :attr:`resolution` x 1000.
        plot : bool
            If true, the profile is plotted, including the intersection
            result.

        Returns
        -------
        float
            retrieved horizontal distance of intersection point.
        float
            error of retrieved horizontal distance of intersection point.
        GeoPoint
            location of intersect
        ndarray
            elevations along viewing direction vector arising from input
            elevation and altitude offset as well as the azimuth of the
            profile
        Axes or None
            matplotlib axes instance if input arg `plot` is True, else None.
        """

        if min_dist is None:
            min_dist = self.dist_hor*.01
        if max_diff is None:
            max_diff = self.resolution * 1000
        view_elevations = self.get_altitudes_view_dir(elev_angle,
                                                      view_above_topo_m)

        #: determine the difference signal
        diff_signal = view_elevations - self.profile
        #: First condition: consider only points in certain distance from observer
        cond1 = self.dists > min_dist
        #: Second condition: consider only "close to zero" points (this might be
        #: redundant, I'll leave it for now because it works)
        cond2 = abs(diff_signal) < max_diff

        #: create array with all distances matching the 2 conditions
        dists_0 = self.dists[cond1 * cond2]
        #: relax condition 2 if nothing was found
        if len(dists_0) == 0:
            raise IntersectNotFound(
                'could not establish initial array of candidate distances '
                'for retrieval of intersection point, you might succeed by '
                'setting or increasing the value of input parameter max_diff')

        try:
            #: get all diff vals matching the 2 conditions
            diff_vals = diff_signal[cond1 * cond2]
            #: get the index of the first zero crossing
            first_idx = np.where(np.diff(np.sign(diff_vals)))[0][0]

            #: get distance and local tolerance value
            dist, tol = dists_0[first_idx], self.resolution * local_tolerance
            #: make mask to access all distances within tolerance range
            cond4 = np.logical_and(dist - tol <= dists_0, dist + tol >= dists_0)
            #: estimate distance error from standard deviation of all
            #: distances within tolerance range
            dist_err = dists_0[cond4].std()

            #: create geopoint corresponding to intersection
            intersect = self._observer_topogrid.offset(self.azimuth, dist)

            #: set the altitude of this point (retrieved from profile)
            intersect.altitude = self.get_altitude_at_distance(dist) # / 1000.

        except IndexError:
            from traceback import format_exc
            raise IntersectNotFound(
                f'No intersections could be detected, traceback:\n'
                f'{format_exc()}')

        ax = None
        if plot:
            ax = self._plot_intersect_search_result(view_elevations, dist)
            ax.set_title("Azim: %.1f, Elev: %.1f, Intersect @ "
                "dist =%.1f km" %(self.azimuth, elev_angle, dist))

        return dist, dist_err, intersect, view_elevations, ax

    def find_horizon_elev(self, elev_start=0.0, elev_stop=60.0, step_deg=0.1,
                          **kwargs):
        """Find first elevation angle which does not intersect with topo

        Scans towards zenith elevation within input elevation range
        in elevation steps provided via input arg `step_deg`, starting at
        `elev_start`. For each elevation, an attempt is made to find an
        intersection of the resulting viewing direction vector with the local
        topography using :func:`get_first_intersection`. The horizon
        elevation is identified as the first elevation angle for which no
        intersection can be found.

        Parameters
        ----------
        elev_start : float
            start search elevation angle
        elev_stop : float
            stop search elevation angle
        step_deg : float
            angle step for search (coarser is faster)
        **kwargs:
            additional keyword agruments passed to
            :func:`get_first_intersection`

        Raises
        ------
        IntersectNotFound
            if no intersection with topography can be found for input
            elevation range.

        Returns
        -------
        float
            detected elevation angle of horizon
        list
            list of elevation angles for which intersections with topography
            could be identified.
        list
            corresponding horizontal distances of :attr:`observer` to the
            terrain features.
        """
        elevs = np.arange(elev_start, elev_stop + step_deg, step_deg)
        elev_sects = []
        dists_sects = []

        prev_elev = None
        for elev in elevs:
            try:
                (dist,dist_err,intersect,view_elevations,_) = \
                    self.get_first_intersection(elev, plot=False, **kwargs)
            except IntersectNotFound:
                if prev_elev is None:
                    raise IntersectNotFound(
                        f'failed to determine elevation angle of horizon '
                        f'consider starting at a lower elevation than '
                        f'{elev_start}')
                else:
                    return elev, elev_sects, dists_sects
            dists_sects.append(dist), elev_sects.append(elev)
            prev_elev = elev

    def get_altitude_at_distance(self, dist):
        """Returns altitude at a ceratain distance from observer

        :param float dist: horizontal distance from obsever along profile
        """
        idx = np.argmin(abs(self.dists - dist))
        return self.profile[idx]

    def plot(self, ax = None):
        """Plot the profile into an axis object

        :param ax: matplotlib axis object
        """
        if ax is None:
            from matplotlib.pyplot import subplots
            fig, ax = subplots(1,1)
        ax.fill_between(self.dists, self.profile, facecolor="#994d00",
                        alpha=0.20)
        ax.set_xlabel("Distance [km]")
        ax.set_ylabel("Altitude [m]")
        ax.set_ylim([self.min - .1 * self.alt_range,
                     self.max + .1 * self.alt_range])
        return ax

    def __call__(self, dist):
        """Returns altitude at a certain distance

        :param float dist: distance in km
        """
        return self.get_altitude_at_distance(dist)

    def _init_attrs(self):
        """
        Initiate private attributes required in this class
        """
        self._observer_topogrid = None
        self._endpoint_topogrid = None
        self._coords_topo = None
        self._line = None
        self._profile=None
        self._dists=None

    def _check_set_observer(self, value):
        """Setter helper for :attr:`observer`

        Raises
        ------
        ValueError
            if input is not instance of :class:`GeoPoint`
        """
        if not isinstance(value, GeoPoint):
            raise ValueError('Need instance of geonum.GeoPoint')
        self._observer = value

    def _check_set_endpoint(self, value):
        """Setter helper for :attr:`endpoint`

        Raises
        ------
        ValueError
            if input is not instance of :class:`GeoPoint`
        """
        if not isinstance(value, GeoPoint):
            raise ValueError('Need instance of geonum.GeoPoint')
        self._endpoint = value

    def _check_set_topo_data(self, value):
        """Setter helper for :attr:`topo_data`

        Raises
        ------
        ValueError
            if input is not instance of :class:`TopoData`
        """
        if not isinstance(value, TopoData):
            raise ValueError('Need instance of geonum.TopoData')
        self._topo_data = value

    def _set_latlon_indices_topogrid(self):
        """Calculate and set indeices of start and stop point on topogrid"""
        topo = self.topo_data
        self._coords_topo = {
         'x0' : np.argmin(abs(topo.lons - self.observer.lon.decimal_degree)),
         'y0' : np.argmin(abs(topo.lats - self.observer.lat.decimal_degree)),
         'x1' : np.argmin(abs(topo.lons - self.endpoint.lon.decimal_degree)),
         'y1' : np.argmin(abs(topo.lats - self.endpoint.lat.decimal_degree))
        }

    def _plot_intersect_search_result(self, view_elevations, dist=None):
        ax = self.plot()
        ax.plot(self.dists, view_elevations, label = "Viewing direction")
        try:
            ax.axvline(dist, ls = "--", label = "Intersection")
        except: # pragma: no cover
            pass
        ax.legend(loc="best", fancybox=True, framealpha=0.4)
        return ax
