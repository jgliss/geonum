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

from os.path import exists
from warnings import warn

from matplotlib.pyplot import get_cmap
from numpy import asarray, nanmin, nanmax

from geonum import BASEMAP_AVAILABLE
from geonum.exceptions import OutOfDomain
from geonum.geopoint import GeoPoint
from geonum.geovector3d import GeoVector3D
from geonum.topodata import TopoData
from geonum.topodataaccess import TopoDataAccess

class GeoSetup(object):
    """The GeoSetup class represents a collection of GeoPoints and vectors

    Attributes
    ----------
    id : str
        name of this setup
    points : dict
       contains :class:`GeoPoint` objects assigned to this setup
    vectors : dict
        contains :class:`GeoVector3D` objects assigned to this setup

    Parameters
    ----------
    points : list
        list of :class:`GeoPoint` objects to be included in this setup
    vectors : list
        list of :class:`GeoVector3D` objects to be included in this setup
    lat_ll : float, optional
        lower left latitude of regime
    lon_ll : float, optional
        lower left longitude of regime
    lat_tr : float, optional
        top right latitude of regime
    lon_tr : float, optional
        top right longitude of regime
    id : str
        identification string of this setup
    topo_access_mode : str
        topo data mode, default is SRTM
        (see :class:`TopoDataAccess` for details)
    local_topo_path : str
        local path were topography data (e.g. ETOPO1 data) is stored
    cmap_vecs : str
        String specifying a valid matplotlib colormap supposed to be
        used for drawing :class:`GeoVector3D` objects into overview
        maps
    init_borders : bool
        Whether or not lower left (ll) and top right (tr) border points are
        supposed to be created.

    """

    def __init__(self, points=None, vectors=None, lat_ll=None, lon_ll=None,
                 lat_tr=None, lon_tr=None, id=None,
                 topo_access_mode=None, local_topo_path=None,
                 cmap_vecs=None, init_borders=True):
        if id is None:
            id = 'MyGeoSetup'
        if topo_access_mode is None:
            topo_access_mode = "srtm"
        if cmap_vecs is None:
            cmap_vecs = 'Greens'

        self._cmap_vecs = cmap_vecs
        self.id = id
        self.points = {}
        self.vectors = {}

        self.topo_access_mode = topo_access_mode
        self.local_topo_path = local_topo_path

        self.topo_data = None
        if points is None:
            points = []
        elif isinstance(points, GeoPoint):
            points = [points]
        if not isinstance(points, list):
            raise ValueError('invalid input for points. need GeoPoint or list')

        for pt in points:
            if isinstance(pt, GeoPoint):
                self.add_geo_point(pt, assert_in_domain=False)

        if vectors is None:
            vectors = []
        elif isinstance(vectors, GeoVector3D):
            vectors = [vectors]
        if not isinstance(vectors, list):
            raise ValueError(
                'invalid input for points. need GeoVector3D or list')

        for vec in vectors:
            if isinstance(vec, GeoVector3D):
                self.add_geo_vector(vec)

        if lat_ll is not None:
            try:
                self.new_geo_point(lat_ll, lon_ll, name="ll")
                self.new_geo_point(lat_tr, lon_tr, name="tr")
            except (TypeError, ValueError):
                pass
        if init_borders and not self.borders_set and self.has_points():
            self.set_borders_from_points()

    @property
    def ll(self) -> GeoPoint:
        """Lower left coordinate of domain"""
        if not 'll' in self.points:
            raise AttributeError(
                'Lower left corner of GeoSetup domain is not defined'
            )
        return self.points["ll"]

    @ll.setter
    def ll(self, value):
        if not isinstance(value, GeoPoint):
            raise ValueError("Could not set lower left coordinate in "
                             "GeoSetup: need GeoPoint object")
        self.points["ll"] = value
        if self.has_topo_data():
            print('updated lower-left coordinate, consider reloading '
                  'topographic dataset using method load_topo_data for '
                  'updated domain.')

    @property
    def tr(self) -> GeoPoint:
        """Top right coordinate of domain"""
        if not 'tr' in self.points:
            raise AttributeError(
                'Top right corner of GeoSetup domain is not defined'
            )
        return self.points["tr"]

    @tr.setter
    def tr(self, value):
        if not isinstance(value, GeoPoint):
            raise ValueError("Could not set top right coordinate in "
                             "GeoSetup: need GeoPoint object")
        self.points["tr"] = value
        if self.has_topo_data():
            print('updated top-right coordinate, consider reloading '
                  'topographic dataset using method load_topo_data for '
                  'updated domain.')

    @property
    def lon_ll(self) -> float:
        """Longitude in decimal degrees of lower left coordinate of domain"""
        return self.ll.longitude

    @property
    def lat_ll(self) -> float:
        """Latitude in decimal degrees of lower left coordinate of domain"""
        return self.ll.latitude

    @property
    def lon_tr(self) -> float:
        """Longitude in decimal degrees of top right coordinate of domain"""
        return self.tr.longitude

    @property
    def lat_tr(self) -> float:
        """Latitude in decimal degrees of top right coordinate of domain"""
        return self.tr.latitude

    @property
    def delta_lon(self) -> float:
        """Longitude range of domain (in decimal degrees)"""
        return abs(self.lon_tr - self.lon_ll)

    @property
    def delta_lat(self) -> float:
        """Latitude range of domain (in decimal degrees)"""
        return abs(self.lat_tr - self.lat_ll)

    @property
    def borders_set(self) -> bool:
        """
        Boolean specifying whether domain borders are set or not
        """
        return True if self.has_point('ll') and self.has_point('tr') else False

    @property
    def center_coordinates(self) -> tuple:
        """Lat / Lon coordinates of center of domain"""
        return (self.lat_ll + self.delta_lat / 2.,
                self.lon_ll + self.delta_lon / 2.)

    @property
    def magnitude(self):
        """Returns dimension (in km) of area covered by this setup"""
        return (self.tr - self.ll).norm

    @property
    def topo_access(self):
        """Topography data access class"""
        return TopoDataAccess(self.topo_access_mode,
                              self.local_topo_path)

    @property
    def cmap_vecs(self):
        """Default colormap used for drawing vectors"""
        return get_cmap(self._cmap_vecs)

    def has_points(self):
        """
        Determine whether any points are set in this GeoSetup

        Returns
        -------
        bool
            True if one or more GeoPoints are registered, else False
        """
        if len(self.points) > 0:
            return True
        return False

    def has_topo_data(self):
        """
        Check whether topographic data is assigned or not

        Returns
        -------
        bool
            True, if topo data is available, else not

        """
        return True if isinstance(self.topo_data, TopoData) else False

    def set_local_topo_path(self, p):
        """Sets local path for Etopo1 data files can be found

        Note
        ----
        The default topo mode is "srtm" which provides online access, so
        it is not mandatory to provide topography data locally. However,
        the SRTM model has no global coverage, so there might be need to
        use another of the provided topo modes and provide the respective
        files locally.

        Parameters
        ----------
        p : str
            new search path for topography data
        """
        if not exists(p):
            raise FileExistsError("Input path does not exist")
        self.local_topo_path = p

    def change_topo_mode(self, new_mode=None, local_path=None):
        """Change the current mode for topography data access

        Parameters
        ----------
        new_mode : str
            new topo access mode, default to "srtm"
        local_path : str, optional
            if not None and valid, update local topo access

        """
        if new_mode is None:
            new_mode = "srtm"
        if local_path is not None:
            self.set_local_topo_path(local_path)
        self.topo_access_mode = new_mode

    def get_topo(self) -> TopoData:
        """Get current topo data"""
        if not isinstance(self.topo_data, TopoData):
            self.load_topo_data()
        return self.topo_data

    def load_topo_data(self):
        """Load topography data

        Note
        ----

        The loaded :class:`TopoData` object will also be set in all
        :class:`GeoPoint` objects belonging to this setup

        Parameters
        ----------
        topo_access_mode : str, optional
            topo dataset that is supposed to be used. If None, then
            :attr:`topo_access_mode` is used.

        Returns
        -------
        TopoData
            topographic data (is also assigned to :attr:`topo_data`).
        """
        if "ll" not in self.points:
            self.set_borders_from_points()
        self.topo_data = self.topo_access.get_data(self.ll.latitude,
                                                   self.ll.longitude,
                                                   self.tr.latitude,
                                                self.tr.longitude)
        for p in list(self.points.values()):
            p.set_topo_data(self.topo_data)

        return self.topo_data


    def has_point(self, name):
        """Checks if point with input name exists

        Parameters
        ----------
        name : str
            name of GeoPoint to be checked

        Returns
        -------
        bool
            Whether or not point exists
        """
        if name in self.points:
            return True
        return False

    def has_vector(self, name):
        """Checks if vector with input name exists

        Parameters
        ----------
        name : str
            name of vector

        Returns
        -------
        bool
            whether or not such a vector with input name exists
        """
        if name in self.vectors:
            return True
        return False

    def contains_coordinate(self, lat, lon):
        """
        Check if input coordinate is within domain of this GeoSetup

        Parameters
        ----------
        lat : float
            latitude of coordinate
        lon : float
            longitude of coordinate

        Returns
        -------
        bool

        """
        latok = self.lat_ll <= lat <= self.lat_tr
        lonok = self.lon_ll <= lon <= self.lon_tr
        return True if latok and lonok else False

    def add_geo_point(self, pt, assert_in_domain=True,
                      overwrite_existing=False) -> None:
        """Add a GeoPoint to this collection

        Parameters
        ----------
        pt : GeoPoint
            point to be added
        assert_in_domain : bool
            if True, a check is performed whether the input point is within
            domain or not (does not apply if input point is "ll" or "tr", that
            is, one of the points defining the domain itself). Defaults to
            True.
        overwrite_existing : bool
            if True and a point with the same name already exists in this
            setup, then the existing point will be overwritten with input
            point, else an exception is raises. Defaults to False.

        Raises
        ------
        OutOfDomain
            if paramerter `assert_in_domain` is True and input point is not within
            current domain.
        ValueError
            if parameter `overwrite_existing` is False and a point with the input n
            name already exists in this setup.

        Returns
        -------
        None
        """
        if not isinstance(pt, GeoPoint):
            raise ValueError(f'invalid input: {pt}. Need GeoPoint')
        if pt.name == 'll':
            self.ll = pt
        elif pt.name == 'tr':
            self.tr = pt
        elif pt.name in self.points and not overwrite_existing:
            raise ValueError(
                f'GeoPoint with name {pt.name} already exists in GeoSetup. '
            )
        elif self.borders_set and assert_in_domain and not \
                self.contains_coordinate(pt.latitude, pt.longitude):

            raise OutOfDomain(f'{pt} is not within domain of GeoSetup')

        else:
            self.points[pt.name] = pt
            if (isinstance(self.topo_data, TopoData) and
                    not isinstance(pt.topo_data, TopoData)):
                try:
                    pt.set_topo_data(self.topo_data)
                except OutOfDomain:
                    pass

    def add_geo_points(self, *pts, assert_in_domain=False) -> None:
        """Add multiple GeoPoints to the collection

        Parameters
        ----------
        *pts
            points to add
        assert_in_domain : bool
            if True, check assert that each point is within domain

        Returns
        -------
        None
        """
        for pt in pts:
            self.add_geo_point(pt, assert_in_domain)

    def add_geo_vector(self, vec, overwrite_existing=True) -> None:
        """Add :class:`GeoVector3D` to this collection

        Parameters
        ----------
        vec : GeoVector3D
            vector to be added
        overwrite_existing : bool
            If True and if vector with same name already exists, the former one

        Raises
        ------
        ValueError
            if input is not GeoVector3D

        """
        if not isinstance(vec, GeoVector3D):
            raise ValueError(f'invalid input: {vec}. Need GeoVector3D')
        if vec.name in self.vectors and not overwrite_existing:
            raise ValueError(f'Vector with name {vec.name} already exists.')
        self.vectors[vec.name] = vec

    def add_geo_vectors(self, *vecs) -> None:
        """Add multiple GeoVector3D objects to the collection

        Parameters
        ----------
        *vecs
            Instances of :class:`GeoVector3D` to be added
        """
        for vec in vecs:
            self.add_geo_vector(vec)

    def delete_geo_point(self, name) -> None:
        """Remove one of the geo points from the collection

        Parameters
        ----------
        name : str
            name of the point

        Raises
        ------
        ValueError
            if no point with with input name exists in this setup

        Returns
        -------
        None
        """
        if not name in self.points:
            raise ValueError(f'no such GeoPoint ({name}) in GeoSetup')
        del self.points[name]

    def delete_geo_vector(self, name):
        """Remove one of the vectors from the collection

        Parameters
        ----------
        name : str
            name of the vector

        Raises
        ------
        ValueError
            if no vector with with input name exists in this setup

        Returns
        -------
        None
        """
        if not name in self.vectors:
            raise ValueError(f'no such GeoVector3D ({name}) in GeoSetup')
        del self.vectors[name]

    def new_geo_point(self, *args, **kwargs) -> None:
        """Create new geo_point and add to collection

        Create GeoPoint using input args and kwargs and then parse that
        GeoPoint to :func:`add_geo_point`.

        Parameters
        ----------
        *args
            input args passed to init function of GeoPoint
        **kwargs
            input keyword args passed to init function of GeoPoint

        Returns
        -------
        None
        """
        self.add_geo_point(GeoPoint(*args, **kwargs))

    def _all_lats_lons(self) -> tuple:
        """Get list of all latitude and longitude coordinates of all points

        Returns
        ----
        ndarray
            latitude coordinates
        ndarray
            longitude coordinates
        """
        lats, lons = [], []
        for p in self.points.values():
            lats.append(p.latitude)
            lons.append(p.longitude)
        return (asarray(lats), asarray(lons))

    def set_borders_from_points(self, extend_km=1, to_square=True) -> None:
        """Set lower left (ll) and top right (tr) corners of domain

        The domain is inferred from all points associated with this
        setup.

        Parameters
        ----------
        extend_km : float
            extend range from the outermost points by this number in km
        to_square : bool
            extend the shorter base side to the size of the longer one (
            quadratic range)

        Raises
        ------
        AttributeError
            if no points are available in the setup.
        """
        lats, lons = self._all_lats_lons()
        if not len(lats) > 0:
            raise AttributeError('Cannot initiate range coordinates in empty '
                                 'GeoSetup, please add at least one point')

        lat_ll, lon_ll, lat_tr, lon_tr = (nanmin(lats), nanmin(lons),
                                          nanmax(lats), nanmax(lons))

        pll, ptr = GeoPoint(lat_ll, lon_ll, 0.0), GeoPoint(lat_tr, lon_tr, 0.0)

        if to_square:
            v = ptr - pll
            add = (abs(v.dx) - abs(v.dy)) / 2
            if add > 0:  # E/W extend (dx) is smaller than N/S extend (dy)
                pll = pll.offset(azimuth=180, dist_hor=add)
                ptr = ptr.offset(azimuth=0, dist_hor=add)
            else:
                pll = pll.offset(azimuth=270, dist_hor=-add)
                ptr = ptr.offset(azimuth=90, dist_hor=-add)

        ll = pll.offset(azimuth=-135,
                        dist_hor=float(extend_km),
                        name='ll')
        self.add_geo_point(ll, assert_in_domain=False,
                           overwrite_existing=True)

        tr = ptr.offset(azimuth=45,
                        dist_hor=float(extend_km),
                        name='tr')
        self.add_geo_point(tr, assert_in_domain=False,
                           overwrite_existing=True)

    def points_close(self, p, radius=None):
        """
        Find all GeoPoints within a certain radius around input point

        The search is performed against all points defined in this GeoSetup.

        Parameters
        ----------
        p : GeoPoint
            location for which the search is performed
        radius : float, optional
            radius (in km) specifying considered range around point (is set
            to 10% of the magnitude of this setup if unspecified)

        Returns
        -------
        list
            names of all points that are in vicinity around input point
        """
        if not isinstance(p, GeoPoint):
            raise ValueError('invalid input, need GeoPoint')
        if radius == None:
            radius = self.magnitude * .1
        names = []
        for pt in list(self.points.values()):
            if not pt is p and (pt - p).magnitude < radius:
                names.append(pt.name)
        print(f"Found {len(names):d} points within radius of {radius:.1f} km "
              f"of point {p.name}")
        return names

    @staticmethod
    def create_test_setup() -> 'GeoSetup':
        """Initiate example test data set

        Returns
        -------
        GeoSetup
            example setup
        """
        gs = GeoSetup()
        source = GeoPoint(latitude=37.751005, longitude=14.993435,
                          name="Etna", auto_topo_access=True)
        instrument = GeoPoint(latitude=37.765755, longitude=15.016696,
                              name="Observatory", auto_topo_access=True)
        gs.add_geo_points(source, instrument)
        gs.set_borders_from_points()
        plume = GeoVector3D(azimuth=83, dist_hor=gs.magnitude,
                            elevation=0, anchor=source, name="plume")
        view_dir = GeoVector3D(azimuth=160, dist_hor=gs.magnitude,
                               elevation=8, anchor=instrument, name="cfov")

        gs.add_geo_vectors(plume, view_dir)
        return gs

    def create_map(self, *args, **kwargs):  # pragma: no cover
        """Create a Basemap object for this regime"""
        if not BASEMAP_AVAILABLE:
            raise ImportError("Cannot create map: Basemap library is not available")

        from geonum.mapping import Map
        if not isinstance(self.topo_data, TopoData):
            self.load_topo_data()
        if not "projection" in kwargs and self.magnitude < 150:
            kwargs["projection"] = "lcc"
        if "llcrnrlon" not in kwargs:
            kwargs["llcrnrlat"] = self.lat_ll
            kwargs["llcrnrlon"] = self.lon_ll
            kwargs["urcrnrlat"] = self.lat_tr
            kwargs["urcrnrlon"] = self.lon_tr
            kwargs["lat_0"], kwargs["lon_0"] = self.center_coordinates

        m = Map(*args, **kwargs)
        m.set_topo_data(self.topo_data)
        return m

    def plot_2d(self, draw_all_points=True, draw_all_vectors=True,
                draw_topo=True, draw_coastline=True, draw_mapscale=True,
                draw_legend=True, *args, **kwargs):  # pragma: no cover
        """
        To be updated
        """
        if not BASEMAP_AVAILABLE:
            raise ImportError("Cannot create overview map: Basemap module "
                              "is not available")
        if "ax" not in kwargs:
            from matplotlib.pyplot import figure
            fig = figure(figsize=(10, 8))
            ax = fig.add_axes([0.12, 0.15, 0.8, 0.8])
            kwargs["ax"] = ax
        m = self.create_map(*args, **kwargs)
        if draw_coastline:
            try:
                m.drawcoastlines()
            except ValueError:
                pass  # see https://github.com/jgliss/geonum/issues/5
        if draw_topo:
            m.draw_topo(insert_colorbar=True)
            m.draw_topo_contour()
        m.draw_coordinates()
        if draw_mapscale:
            m.draw_mapscale_auto()
        p_close_count = 0
        if draw_all_points:
            dist = self.magnitude * .05
            for pt in list(self.points.values()):
                if not any([pt.name == x for x in ["ll", "tr"]]):
                    m.draw_geo_point_2d(pt)
                    ang = -45
                    num_close = len(self.points_close(pt))
                    if num_close > 0:
                        step = 360. / (4 * num_close)
                        ang = ang - step * p_close_count
                        p_close_count += 1
                    m.write_point_name_2d(pt, dist, ang)

        # create some color indices for colormap
        nums = [int(255.0 / k) for k in range(1, len(self.vectors) + 3)]
        if draw_all_vectors:
            for i, vec in enumerate(self.vectors.values()):
                m.draw_geo_vector_2d(vec,
                                     ls="-",
                                     c=self.cmap_vecs(nums[i]),
                                     label=vec.name)
        if draw_legend:
            try:
                m.legend()
            except Exception as e:
                warn(f"Failed to draw legend in GeoSetup...: {e}")

        return m

    def plot_3d(self, draw_all_points=True, draw_all_vectors=True,
                cmap_topo="Oranges", contour_color="#708090",
                contour_lw=0.2, contour_antialiased=True,
                *args, **kwargs):  # pragma: no cover
        """Create a 3D overview map of the current setup

        Parameters
        ----------
        draw_all_points : bool
            if True, all current GeoPoint objects are plotted, defaults to
            True
        draw_all_vectors : bool
            if True, all current GeoVector3D objects are plotted, defaults to
            True
        cmap : str
            string ID of the colormap used to plot the topographic data,
            defaults to "Oranges"
        contour_color : str
            string specifying color of contour lines colors of contour lines
            (default: "#708090")
        contour_lw :
            width of drawn contour lines, defaults to 0.5, use 0 if you do not
            want contour lines inserted
        contour_antialiased : bool
            apply antialiasing to surface plot of topography, defaults to False
        *args :
            additional non-keyword parameters (passed to `basemap
            <http://matplotlib.org/basemap/api/basemap_api.html#mpl
            _toolkits.basemap.Basemap>`_)
        **kwargs :
            additional keyword parameters (passed to `basemap
            <http://matplotlib.org/basemap/api/basemap_api.html#mpl
            _toolkits.basemap.Basemap>`_)

        Returns
        -------
        Map
            plotted 3D basemap
        """
        if not BASEMAP_AVAILABLE:
            raise ImportError("Cannot create overview map: Basemap module "
                              "is not available.")
        m = self.create_map(*args, **kwargs)
        m.draw_topo_3d(cmap=cmap_topo, contour_color=contour_color,
                       contour_lw=contour_lw,
                       contour_antialiased=contour_antialiased)
        if draw_all_points:
            zr = self.topo_data.alt_range * 0.05
            alts = []
            for name, pt in self.points.items():
                if not any([name == x for x in ["ll", "tr"]]):
                    try:
                        add_alt = 0  # in m
                        for alt in alts:
                            if abs(alt - pt.altitude) < zr:
                                add_alt = 3 * zr
                                print("Add " + str(add_alt))

                        pt.plot_3d(m, add_name=True, dz_text=zr + add_alt)
                        alts.append(pt.altitude)

                    except Exception as e:
                        warn("Point %s could not be drawn: %s"
                             % (pt.name, repr(e)))

        if draw_all_vectors:
            nums = [int(255.0 / k) for k in range(1, len(self.vectors) + 3)]
            for i, vec in enumerate(self.vectors.values()):
                try:
                    m.draw_geo_vector_3d(vec, label=vec.name,
                                         c=self.cmap_vecs(nums[i]),
                                         ls="-",
                                         **kwargs)
                except Exception as e:
                    warn("Vector %s could not be drawn: %s"
                         % (vec.name, repr(e)))

        return m
