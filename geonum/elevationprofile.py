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
from warnings import warn
from scipy.interpolate import interp1d

from geonum.geopoint import GeoPoint

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
    topo_data : TopoData 
        topography data object
    observer : :obj:`GeoPoint`
        starting point of profile 
    endpoint : :obj:`GeoPoint`
        stop point of profile 
    interpolate : bool
        if True, the profile is interpolated to a certain horizontal resolution
    resolution : float
        desired grid resolution in m for interpolation. Interpolation is 
        performed if :attr:`interpolate` is True and if the actual 
        resolution of the topo data is smaller than input, else, nothing is done
    itp_type : str
        interpolation type (e.g. "linear", "cubic")
    mapping_opts
        additional keyword args that are passed to :func:`det_profile`
    """        
    def __init__(self, topo_data, observer, endpoint, interpolate=True, 
                 resolution=5.0, itp_type="linear", **mapping_opts):
     
        self.topo_data = topo_data
        self.observer = observer #: coordinate of observer (start of profile)
        self.endpoint = endpoint #: coordinate of endpoint of profile
        
        self._observer_topogrid = None #: closest to observer on topo data grid
        self._endpoint_topogrid = None #: closest to endpoint in topo data grid
        
        # In the following parameters the results will be stored
        self.line = None
        self.profile = None
        self.dists = None

        try:
            self.det_profile(interpolate, resolution, itp_type, **mapping_opts)
        except Exception as e:
            warn("Failed to compute elevation profile. Error msg: %s"
                 %repr(e))
        
    @property
    def dist_hor(self):
        """Returns the horizontal distance between the 2 points"""
        return (self._endpoint_topogrid - self._observer_topogrid).dist_hor
        
    
    @property
    def azimuth(self):
        """Returns the azimuth angle of the profile"""
        return (self._endpoint_topogrid - self._observer_topogrid).azimuth
        
    @property
    def resolution(self):
        """Get profile x (distance) resolution (averaged from distance array)
        
        Note
        ----
        Only works if profile was already determined
            
        """
        return abs((self.dists[1:] - self.dists[:-1]).mean())
    
    @property
    def gradient(self):
        """Return gradient of profile
        
        Uses numpy function ``gradient``
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
        """Return minimum altitude in profile"""
        return np.nanmin(self.profile)
        
    @property
    def max(self):
        """Return maximum altitude in profile"""
        return np.nanmax(self.profile)
    
    @property
    def alt_range(self):
        """Return altitude range of profile
        """
        return self.max - self.min
    
    def slope_angles(self, decimal_degrees=True):
        """Returns slope angle of profile (in each sample point)
        
        :param bool decimal_degrees: rad or degrees (default True)
        """
        a = np.tan(self.slope)
        if decimal_degrees:
            a = np.rad2deg(a)
        return a
        
    def slope_angle(self, dist):
        """Returns slope angle of profile at input distance
        
        :param float dist: distance in km
        """
        idx = np.argmin(abs(self.dists - dist))
        return self.slope_angles()[idx]
        
    def det_profile(self, interpolate=True, resolution=5.0, itp_type="linear",
                    **mapping_opts):
        """Determines the elevation profile
        
        Searches the closest tiles in the topo data grid for both observer and
        endpoint and based on these two points the elevation profile is
        extracted using a :class:`LineOnGrid` object (which extracts altitudes 
        from the topo data along the connection vector of the 2 points using 
        2D spline intepolation)
        
        Parameters
        ----------
        interpolate : bool
            if True, the profile is interpolated to a certain horizontal resolution
        resolution : float
            desired grid resolution in m for interpolation. Interpolation is 
            performed if :attr:`interpolate` is True and if the actual 
            resolution of the topo data is smaller than input, else, nothing is done
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
        from geonum.processing import LineOnGrid
        data = self.topo_data
        idx_lon_0 = np.argmin(abs(data.lons - self.observer.lon.decimal_degree))
        idx_lat_0 = np.argmin(abs(data.lats - self.observer.lat.decimal_degree))
        idx_lon_1 = np.argmin(abs(data.lons - self.endpoint.lon.decimal_degree))
        idx_lat_1 = np.argmin(abs(data.lats - self.endpoint.lat.decimal_degree))
        
        self.line = l = LineOnGrid(idx_lon_0, idx_lat_0, 
                                   idx_lon_1, idx_lat_1)
        
        z = l.get_line_profile(data.data, **mapping_opts)
        self._observer_topogrid = GeoPoint(data.lats[idx_lat_0],
                                           data.lons[idx_lon_0], 
                                           topo_data=data)
        self._endpoint_topogrid = GeoPoint(data.lats[idx_lat_1],
                                           data.lons[idx_lon_1], 
                                           topo_data=data)
        dists = np.linspace(0, self.dist_hor, l.length + 1)
        if interpolate:
            try:
                res0 = (dists[1] - dists[0]) * 1000
                fac = int(np.ceil(res0 / resolution))
                if fac > 1:
                    fz = interp1d(dists, z, kind=itp_type)
                    dists = np.linspace(0, self.dist_hor, l.length * fac)
                    z = fz(dists)
            except Exception as e:
                warn("Failed to perform interpolation of retrieved elevation "
                     "profile. Error msg: %s" %repr(e))
        self.dists = dists
        self.profile = z
        return z
     
    def get_altitudes_view_dir(self, elev_angle, view_above_topo_m=1.5):
        """Get vector containing altitudes for a viewing direction 
        
        The viewing direction is specified by the azimuth angle of the 
        connection vector between observer and endpoint, the elevation angle
        needs to be specified via the input parameters, and the start point
        (first elevation value) corresponds to the altitude at the observer
        position plus an offset in m which can be specified.
        
        :param float elev_angle: elevation angle of viewing direction
        :param float view_above_topo_m (1.5): altitude offset of start point 
            in m
        :return:
            - vector with altitude values (same length as ``self.profile``)
            
        """
        return (1000 * np.tan(np.radians(elev_angle)) * self.dists + 
                self.profile[0] + view_above_topo_m)
    
    def find_horizon_elev(self, elev_start=0.0, elev_stop=60.0, step_deg=0.1,
                          **kwargs):
        """Find first elevation angle which does not intersect with topo
        
        :param float elev_start: start search elevation angle
        :param float elev_stop: stop search elevation angle
        :param float step_deg: angle step for search (coarser is faster)
        :param **kwargs: additional keyword agruments passed to 
            :func:`get_first_intersection`
        """
        elevs = np.arange(elev_start, elev_stop + step_deg, step_deg)
        elev_sects = []
        dists_sects = []
        for elev in elevs:
            (dist, 
             dist_err, 
             intersect, 
             view_elevations, 
             _) = self.get_first_intersection(elev, **kwargs)
            if dist is None:
                return elev, elev_sects, dists_sects
            else:
                dists_sects.append(dist), elev_sects.append(elev)
        raise Exception("Unexpected exception..")
        
    def get_first_intersection(self, elev_angle, view_above_topo_m=1.5,
                               min_dist=None, local_tolerance=3, plot=False):
                                
        """Finds first intersection of a viewing direction with topography
        
        Start point of the viewing vector is the observer position (or more
        accurately, the center position of the closest topography tile in 
        the topography data set). The relative altitude of the start point 
        with respect to the topography altitude at the observer position 
        can be controlled on input (arg ``view_above_topo_m``) as well as 
        the elevation angle of the viewing direction. The azimuth angle 
        corresponds to the profile azimuth (i.e. azimuth between observer 
        position and endpoint of this profile).
        
        The signal analysed for finding the intersection is a vector 
        containing relative elevation positions of the viewing direction 
        with respect to the profile altitude::
        
            diff_signal = self.get_altitudes_view_dir - self.profile
            
        The first intersection of the viewing direction with the topography
        is identified by the first zero crossing of this ``diff_signal``.
        
        :param float elev_angle: elevation angle of viewing direction (in 
            decimal degrees)
        :param int view_above_topo_m: altitude offset of start point 
            (``observer``) in m (default: 1.5)
        :param float min_dist: minimum distance (in km) of first 
            considered intersection with topography from observer. If None, 
            use 1% of distance between ``observer`` and ``endpoint``.
        :param float local_tolerance: tolerance factor to estimate 
            distance uncertainty based on topo grid resolution  (default: 3)
        :param bool plot: creates a plot of the profile including the
            intersection information

        
        """
         
        if min_dist == None:
            min_dist = self.dist_hor*.01
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
        dist, dist_err, intersect = None, None, None
        #: relax condition 2 if nothing was found 
        if not len(dists_0) > 0:
            max_diff = self.resolution * 1000 * local_tolerance
            cond2 = abs(diff_signal) < max_diff
            dists_0 = self.dists[cond1 * cond2]
                   
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
           
        except IndexError as e:
            print(("No intersections could be detected, err: %s" %repr(e)))
        
        ax = None
        if plot:
            ax = self._plot_intersect_search_result(view_elevations, dist)
            if dist == None:
                dist = np.nan
            ax.set_title("Azim: %.1f, Elev: %.1f, Intersect @ "
                "dist =%.1f km" %(self.azimuth, elev_angle, dist))
            
        return dist, dist_err, intersect, view_elevations, ax
    
    def _plot_intersect_search_result(self, view_elevations, dist=None):
        ax = self.plot()
        ax.plot(self.dists, view_elevations, label = "Viewing direction")
        try:
            ax.axvline(dist, ls = "--", label = "Intersection")
        except:
            pass
        ax.legend(loc="best", fancybox=True, framealpha=0.4)
        return ax
        
        
    def get_altitude_at_distance(self, dist):
        """Returns altitude at a ceratain distance from observer
        
        :param float dist: horizontal distance from obsever along profile        
        """
        idx = np.argmin(abs(self.dists - dist))
        return self.profile[idx]
        
    def __call__(self, dist):
        """Returns altitude at a certain distance
        
        :param float dist: distance in km
        """
        return self.get_altitude_at_distance(dist)
        
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