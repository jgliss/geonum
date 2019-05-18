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
This module contains the GeoSetup class, a high level object for managing 
muliple geo point and geo vector objects.
"""
from geonum import BASEMAP_AVAILABLE
if BASEMAP_AVAILABLE:
    from geonum.mapping import Map
from numpy import asarray, nanmin, nanmax
from os.path import exists
from traceback import print_exc
from warnings import warn

from geonum.geopoint import GeoPoint
from geonum.geovector3d import GeoVector3D
from geonum.topodataaccess import TopoDataAccess
from geonum.topodata import TopoData

class GeoSetup(object):
    """The GeoSetup class represents a collection of GeoPoints and vectors
    
    Attributes
    ----------
    id : str
        name of this setup
    points : list
        list of :class:`GeoPoint` objects assigned to this setup
    vectors : list
        list of :class:`GeoVector3D` objects assigned to this setup
        
    Parameters
    ----------
    points : list
        list of :class:`GeoPoint` objects to be included in this setup
    vectors : list
        list of :class:`GeoVector3D` objects to be included in this setup
    lat_ll : :obj:`float`, optional
        lower left latitude of regime
    lon_ll : :obj:`float`, optional
        lower left longitude of regime
    lat_tr : :obj:`float`, optional
        top right latitude of regime
    lon_tr : :obj:`float`, optional
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
    
        
    """
    def __init__(self, points=[], vectors=[], lat_ll=None, lon_ll=None, 
                 lat_tr=None, lon_tr=None, id="MyGeoSetup", 
                 topo_access_mode="srtm", local_topo_path=None, 
                 cmap_vecs="Greens"):

        self.id = id
        self.points = {}
        self.vectors = {}
        from matplotlib.pyplot import get_cmap
        try:
            cmap = get_cmap(cmap_vecs)
        except:
            cmap = get_cmap("Greens")
        self.cmap = cmap
        
        self.topo_access_mode = topo_access_mode
        self.local_topo_path = local_topo_path
        
        self.topo_data = None
        try:
            iter(points)
        except TypeError:
            if isinstance(points, GeoPoint):
                points = [points]
            else:
                print("Invalid input points: ", points) 
                points = []
        try:
            iter(vectors)
        except:
            if isinstance(points, GeoPoint):
                points = [points]
            else:
                print("invalid input vectors: ", vectors)
                vectors = []
                
        for pt in points:
            if isinstance(pt, GeoPoint):
                self.add_geo_point(pt)
        for vec in vectors:
            if isinstance(vec, GeoVector3D):
                self.add_geo_vector(vec)
                
        self.set_borders_from_points()
        #temporary creation of points ll and tr in case there is some valid input
        try:
            self.new_geo_point(lat_ll, lon_ll, name="ll")
            self.new_geo_point(lat_tr, lon_tr, name="tr")
        except (TypeError, ValueError):
            pass
      
    @property
    def topo_access(self):
        """Topograph data access class"""
        return TopoDataAccess(self.topo_access_mode, 
                              self.local_topo_path) 
    def has_points(self):
        """Returns True, if this setup includes GeoPoints, False if not"""
        if not bool(self.points):
            return False
        return True
    
    def create_test_data(self):
        """Create exemplary test data set"""
        source = GeoPoint(37.751005,  14.993435, name="Etna")
        instrument = GeoPoint(37.765755,  15.016696, name="Observatory")
        self.add_geo_points(source, instrument)
        self.set_borders_from_points()
        plume = GeoVector3D(azimuth=83, dist_hor = self.magnitude,
                            elevation=0, anchor=source, name="plume")
        view_dir = GeoVector3D(azimuth=160, dist_hor=self.magnitude,
                               elevation=8, anchor=instrument, name="cfov")
                                
        self.add_geo_vectors(plume, view_dir)
        
    def set_local_topo_path(self, p):
        """Sets local path for Etopo1 data files can be found
        
        Note
        ----
        The default topomode is "srtm" which provides online access, so
        it is not mandatory to provide topography data locally. However, 
        the SRTM model has no global coverage, so there might be need to 
        use another of the provided topomodes and provide the respective
        files locally.
        
        Parameters
        ----------
        p : str
            new search path for topography data
        """
        if not exists(p):
            raise IOError("Input path does not exist")
        self.topo_access.local_path = p
    
    def change_topo_mode(self, new_mode="srtm", local_path=None):
        """Change the current mode for topography data access
        
        Parameters
        ----------
        new_mode : str
            new topo access mode
        local_path : :obj:`str`, optional 
            if not None and valid, update local topo access
            
        """
        if local_path is not None and exists(local_path):
            self.load_topo_path = local_path
        self.topo_access_mode = new_mode
    
    def get_topo(self):
        """Get current topo data"""
        if not isinstance(self.topo_data, TopoData):
            self.load_topo_data()
        return self.topo_data
        
    def load_topo_data(self):
        """Load topography data 
        
        .. note:: 
        
            The loaded :class:`TopoData` object will also be set in all 
            :class:`GeoPoint` objects belonging to this setup
            
        """
        if "ll" not in self.points:
            self.set_borders_from_points()
        self.topo_data = self.topo_access.get_data(self.ll.latitude,
                                                   self.ll.longitude,
                                                   self.tr.latitude,
                                                   self.tr.longitude)
        for p in list(self.points.values()):
            p.set_topo_data(self.topo_data)
        
    @property
    def ll(self):
        """Return lower left point of topo data range """
        try:
            return self.points["ll"]
        except AttributeError:
            print("Lower left corner (GeoPoint) not yet defined in GeoSetup")
    
    @ll.setter
    def ll(self, value):
        if not isinstance(value, GeoPoint):
            raise TypeError("Could not set lower left coordinate in "
                "GeoSetup: need GeoPoint object")
        self.points["ll"] = value
        
    @property
    def tr(self):
        """Return lower left point of topo data range"""
        try:
            return self.points["tr"]
        except AttributeError:
            print("Top right corner (GeoPoint) not yet defined in GeoSetup")
            pass
    
    @tr.setter
    def tr(self, value):
        if not isinstance(value, GeoPoint):
            raise TypeError("Could not set top right coordinate in "
                            "GeoSetup: need GeoPoint object")
        self.points["tr"] = value
        
    @property
    def lon_ll(self):
        """Lower left corner of object regime"""
        return self.ll.longitude
    
    @property
    def lat_ll(self):
        """Lower left corner of object regime"""
        return self.ll.latitude
        
    @property
    def lon_tr(self):
        """Lower left corner of object regime"""
        return self.tr.longitude
    
    @property
    def lat_tr(self):
        """Lower left corner of object regime"""
        return self.tr.latitude
    
    @property
    def delta_lon(self):
        """Returns longitude range"""
        return abs(self.lon_tr - self.lon_ll)
        
    @property    
    def delta_lat(self):
        """Returns latitude range"""
        return abs(self.lat_tr - self.lat_ll)
        
    @property
    def center_coordinates(self):
        """Lat / Lon coordinates of center of data"""
        return (self.lat_ll + self.delta_lat / 2., 
                self.lon_ll + self.delta_lon / 2.)
        
    def add_geo_points(self, *args):
        """Add multiple GeoPoints to the collection
        
        :param *args: arbitrary amount of new geo points        
        """
        for arg in args:
            self.add_geo_point(arg)
    
    def has_point(self, name):
        """Checks if point with input name exists
        
        :param str key: name of point
        :return: bool
        """
        if name in self.points:
            return True
        return False
    
    def has_vector(self, name):
        """Checks if vector with input name exists"""
        if name in self.vectors:
            return True
        return False
        
    def add_geo_point(self, pt):
        """Add :class:`GeoPoint` to this collection
        
        :param GeoPoint pt: the new point
        """
        try:
            if pt.name in self.points:
                print(("Point ID %s already exists in GeoSetup" %(pt.name)))
                pt2 = self.points[pt.name]
                if pt.almost_equal(pt2) and pt.altitude == pt2.altitude:
                    print("Point is unchanged")
                    return
                print("Updating name of existing GeoPoint to: %s_old" %pt.name)
                pt2.name = pt2.name + "_old"
                self.points[pt2.name] = pt2
            self.points[pt.name] = pt
            if not isinstance(pt.topo_data, TopoData):
                pt.set_topo_data(self.topo_data)
        
        except Exception as e:
            print("Geopoint could not be added: " + repr(e))
            
    def set_geo_point(self, p_id, pt):
        """Update an existing GeoPoint in the collection
        
        :param str id: id of existing point in ``self.points``
        :param GeoPoint pt: a new geo_point
        """
        if not isinstance(pt, GeoPoint):
            raise TypeError("Wrong input: " + type(pt))
        self.points[p_id] = pt
        
    def add_geo_vectors(self, *args):
        """Add multiple GeoPoints to the collection"""
        for arg in args:
            self.add_geo_vector(arg)
            
    def add_geo_vector(self, vec):
        """Add :class:`GeoVector3D` to this collection
        
        :param GeoVector3D vec: should be clear ;)
        """
        if not isinstance(vec, GeoVector3D):
            print(("Error adding GeoVector3D, wrong input type, need "
                  ":class:`GeoVector3D` object, input type: %s" %type(vec)))
            return
        if vec.name in self.vectors:
            print(("Vector ID %s already exists in %s" %(vec.name, self)))
            vec2 = self.vectors[vec.name]
            if (vec2.almost_equals(vec) and vec2.dz == vec.dz 
                and vec.anchor == vec2.anchor):
                print("Vector is unchanged")
                return 
            print("Updating name of existing vector to: %s_old" %vec.name)
            vec2.name = vec2.name + "_old"
            self.vectors[vec2.name] = vec2
        self.vectors[vec.name] = vec
    
    def delete_geo_point(self, name):
        """Remove one of the geo_points from the collection
        
        :param str name: name of geo point        
        """
        del self.points[name]
        
    
    def delete_geo_vector(self, name):
        """Remove one of the vectors from the collection
        
        :param str name: name of geo vector
        """
        del self.vectors[name]
            
    def new_geo_point(self, *args, **kwargs):
        """Create new geo_point and add to collection
        
        :param **kwargs: see :class:`GeoPoint` for initiation info
        """
        try:
            self.add_geo_point(GeoPoint(*args, **kwargs))
        except (TypeError, ValueError):
            return
        except:
            raise Exception(print_exc())
    
    def _all_lats_lons(self):
        """Get 2 arrays including all latitudes and all longitudes of all 
        points included in this collection
        
        .. note::
        
            Existing points specifying the regime (i.e. lower left / top 
            right corner) are not considered here
            
        """
        lats, lons = [], []
        for id, p in self.points.items():
            if not any([id == x for x in ["ll","tr"]]):
                lats.append(p.latitude)
                lons.append(p.longitude)
        return asarray(lats), asarray(lons)
    
    @property
    def magnitude(self):
        """Returns dimension (in km) of area covered by this setup"""
        return (self.tr - self.ll).norm
        
    def set_borders_from_points(self, extend_km=1, to_square=True):
        """Set range of setup (lower left and upper right coordinates) 
        considering all points in this collection
        
        :param float extend_km: extend range from the outermost points by 
            this number in km
        :param float to_square (True): extend the shorter base side to the 
            size of the longer one (quadratic range)
        """
        lats, lons= self._all_lats_lons()
        if not len(lats) > 0:
            #print "Borders could not be initiated, no objects found..."
            return False
        
        lat_ll, lon_ll, lat_tr , lon_tr = (nanmin(lats), nanmin(lons),
                                           nanmax(lats), nanmax(lons))
                                           
        pll, ptr = GeoPoint(lat_ll, lon_ll, 0.0), GeoPoint(lat_tr, lon_tr, 0.0)
        
        if to_square:
            v = ptr - pll
            add = (abs(v.dx) - abs(v.dy)) / 2
            if add > 0: #E/W extend (dx) is smaller than N/S extend (dy)
                pll = pll.offset(azimuth = 180, dist_hor = add)
                ptr = ptr.offset(azimuth = 0, dist_hor = add)
            else:
                pll = pll.offset(azimuth = 270, dist_hor = -add)
                ptr = ptr.offset(azimuth = 90, dist_hor = -add)
                 
        self.set_geo_point("ll", pll.offset(azimuth=-135, 
                                            dist_hor=float(extend_km),
                                            name="ll"))
        
        self.set_geo_point("tr", ptr.offset(azimuth=45,
                                            dist_hor=float(extend_km),
                                            name="tr"))
        return True
   
    def create_map(self, *args, **kwargs):     
        """Create a Basemap object for this regime"""
        if not BASEMAP_AVAILABLE:
            raise ImportError("Cannot create map: "
                              "Basemap library is not available")
        
        if not isinstance(self.topo_data, TopoData):
            self.load_topo_data()
        if not "projection" in kwargs and self.magnitude < 150:
            kwargs["projection"] = "lcc"
        if not "llcrnrlon" in kwargs:
            kwargs["llcrnrlat"] = self.lat_ll
            kwargs["llcrnrlon"] = self.lon_ll
            kwargs["urcrnrlat"] = self.lat_tr
            kwargs["urcrnrlon"] = self.lon_tr
            kwargs["lat_0"], kwargs["lon_0"] = self.center_coordinates
       
        m = Map(*args, **kwargs)
        m.set_topo_data(self.topo_data)
        return m
    
    def points_close(self, p, radius=None):
        """Finds all GeoPoints which are within a certain radius around another
        point
        
        :param GeoPoint p: the actual point for which the search is performed
        :param float radius (None): radius (in km) specifying considered range
            around point (is set to 10% of the magnitude of this setup if 
            unspecified)
        :returns:
            - list of point string IDs which are within the specified radius
                around input point
        """
        if radius == None:
            radius = self.magnitude * .1
        ids = []
        for pt in list(self.points.values()):
            if not pt is p and (pt - p).magnitude < radius:
                ids.append(pt.name)
        print(("Found %d points within radius of %.1f km of point %s" 
                                                %(len(ids), radius, p.name)))
        return ids
        
    def plot_2d(self, draw_all_points=True, draw_all_vectors=True,
                draw_topo=True, draw_coastline=True, draw_mapscale=True,
                draw_legend=True, *args, **kwargs):
        """Draw overview map of the current setup
        
        :param bool draw_all_points (True): if true, all points are included
        :param bool draw_all_vectors (True): if true, all vectors (with anchor) 
            are included
        :param bool draw_topo (True): include topography into map
        :param bool draw_coastline (True): include coastline into map
        :param bool draw_mapscale (True): insert map scale
        :param bool draw_legend (True): insert a (draggable) legend
        :param *args: additional non-keyword parameters (passed to `basemap 
            <http://matplotlib.org/basemap/api/basemap_api.html#mpl
            _toolkits.basemap.Basemap>`_)
        :param **kwargs: additional keyword parameters (passed to `basemap 
            <http://matplotlib.org/basemap/api/basemap_api.html#mpl
            _toolkits.basemap.Basemap>`_)
        :return: 
            - :class:`geonum.mapping.Map` object
            
        """
        if not BASEMAP_AVAILABLE:
            raise ImportError("Cannot create overview map: Basemap module "
                              "is not available")
        if not "ax" in kwargs:
            #fig, ax = subplots(1,1)
            from matplotlib.pyplot import figure
            fig = figure(figsize=(10,8))
            ax = fig.add_axes([0.12,0.15,0.8,0.8])
            kwargs["ax"] = ax
        m = self.create_map(*args, **kwargs)
        if draw_coastline:
            m.drawcoastlines()
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
        
        #create some color indices for colormap
        nums = [int(255.0 / k) for k in range(1, len(self.vectors)+3)]
        if draw_all_vectors:
            for i, vec in enumerate(self.vectors.values()):
                m.draw_geo_vector_2d(vec, 
                                     ls="-",
                                     c=self.cmap(nums[i]),
                                     label=vec.name)
        if draw_legend:
            try:
                m.legend()
            except:
                warn("Failed to draw legend in GeoSetup...")
                
        return m
                
    def plot_3d(self, draw_all_points=True, draw_all_vectors=True, 
                cmap_topo="Oranges", contour_color="#708090", 
                contour_lw=0.2, contour_antialiased=True, *args, **kwargs):
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
                        add_alt = 0 #in m
                        for alt in alts:
                            if abs(alt - pt.altitude) < zr:
                                add_alt = 3 * zr
                                print("Add " + str(add_alt))
                                
                        pt.plot_3d(m, add_name = True, dz_text = zr + add_alt)
                        alts.append(pt.altitude)
    
                    except Exception as e:
                        warn("Point %s could not be drawn: %s"
                             %(pt.name, repr(e)))
                        pass
        if draw_all_vectors:
            nums = [int(255.0 / k) for k in range(1, len(self.vectors)+3)]
            for i, vec in enumerate(self.vectors.values()):
                try:
                    m.draw_geo_vector_3d(vec, label=vec.name,
                                         c=self.cmap(nums[i]),
                                         ls="-",
                                         **kwargs)
                except Exception as e:
                    warn("Vector %s could not be drawn: %s"
                         %(vec.name, repr(e)))
                    pass
                
        return m
    
def show_coordinate(geo_point=None, lat_pt=None, lon_pt=None, extend_km=10.0, 
                    *args, **kwargs):
    """Draw overview map for a given point
    
    Parameters
    ----------
    geo_point : GeoPoint
        Geographical location around which overview map is drawn
    lat_pt : float
        Latitude of geographical location around which overview map is 
        drawn (only considered if :attr:`geo_point` is invalid)
    lon_pt : float
        Longitude of geographical location around which overview map is 
        drawn (only considered if :attr:`geo_point` is invalid)
    extend_km : float
        map extend in km around considered geolocation
    *args :
        non-keyword arguments passed to :func:`plot_2d` of the 
        :class:`GeoSetup` instance that is created in order to draw the map
    
    Returns
    -------
    Map
        instance of :class:`geonum.Map`
        
    """
    if not isinstance(geo_point, GeoPoint):
        try: 
            geo_point = GeoPoint(lat=lat_pt, lon=lon_pt)
        except:
            raise TypeError("Invalid input, please provide information "
                            "about location of GeoPoint")
    stp = GeoSetup(points=[geo_point])
    stp.set_borders_from_points(extend_km=extend_km)
    m = stp.plot_2d(*args, **kwargs)
    return m 