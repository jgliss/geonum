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
    from mpl_toolkits.basemap import Basemap
except:
    print('Plotting of maps etc. is deactivated, please install Basemap')

import numpy as np    

from geonum.topodataaccess import TopoDataAccess
from geonum.topodata import TopoData
from geonum.exceptions import TopoAccessError
from geonum.helpers import haversine_formula, shifted_color_map, exponent

from geonum import CV2_AVAILABLE

class Map(Basemap):
    """Basemap object for drawing and plotting (on) a map
    
    This object is initiated as `Basemap <http://matplotlib.org/basemap/
    users/examples.html>`_ and can therefore be used as such. 
    Added functionality mainly includes:
    
    1. Including topography data
    #. 2D and 3D plotting
    #. Handle of :class:`GeoVector3D` and :class:`GeoPoint` objects
    #. Some more high level functionality (i.e. draw text, points, polygons \
       or plot data onto map)

    Note
    ----
    
    The mapping functionality was initally developed for geographical 
    setups on local scales (i.e. several 10 km grids) including 
    handling of high resolution topography data and was not tested to 
    create large maps on global scales. In principle, this should work,
    though.
       
    """ 
    def __init__(self, *args, **kwargs):
        """Initialisation of the map object
        
        :param *args: additional non-keyword parameters (passed to `basemap 
            <http://matplotlib.org/basemap/api/basemap_api.html#mpl
            _toolkits.basemap.Basemap>`_)
        :param **kwargs: additional keyword parameters (passed to `basemap 
            <http://matplotlib.org/basemap/api/basemap_api.html#mpl
            _toolkits.basemap.Basemap>`_)
            
        .. note::
        
            Additional possible input in **kwargs wit to ``Basemap`` 
            objects is "topo_data" which will be set in case it is a valid 
            input (i.e. :class:`TopoData`) and will be used for plotting 
            topography
        """
        super(Map, self).__init__(*args, **kwargs)

        self.topo_data = None
        #IMPORTANT HANDLES
        self.contour_lines = None
        self.contour_filled = None
        
        self.colorbars = {}
        self.points = {}
        self.lines = {}
        self.polygons = {}
        self.texts = {}
        self.plotted_data = {}
        
        self.map_scale = None
        self.meridians = None
        self.parallels = None
        
        self.default_colors = {"contour_lines"  :   "#708090",
                               "land"           :   "#FFE0B2",
                               "water"          :   "#e6e6fa"}
                               
    @property
    def color_land(self):
        """Default color of land tiles"""
        return self.default_colors["land"]

    @property
    def color_water(self):
        """Default color of ocean / sea tiles"""
        return self.default_colors["water"]
                
    def set_topo_data(self, topo):
        """Update topo data
        
        :class TopoData topo: object containing topography data
        """
        if isinstance(topo, TopoData):
            self.topo_data = topo
        
    def load_topo_data(self, mode="srtm", local_path=None):
        """Geonum wrapper for topography access"""
        ta = TopoDataAccess(mode, local_path)
        try:
            self.topo_data = ta.get_data(self.lat_ll, self.lon_ll,
                                         self.lat_tr, self.lon_tr)
            return self.topo_data
        except Exception as e:
            raise TopoAccessError(repr(e))
            
    @property
    def lon_ll(self):
        """Lower left longitude of plotted range"""
        return self.llcrnrlon
    @property
    def lat_ll(self):
        """Lower left latitude of plotted range"""
        return self.llcrnrlat
    @property
    def lon_tr(self):
        """Top right longitude of plotted range"""        
        return self.urcrnrlon
    @property
    def lat_tr(self):
        """Top right latitude of plotted range"""        
        return self.urcrnrlat
    
    @property
    def delta_lon(self):
        """Returns longitude range"""
        return abs(self.lon_tr - self.lon_ll)
        
    @property    
    def delta_lat(self):
        """Returns latitude range"""
        return abs(self.lat_tr - self.lat_ll)
        
    def fill_map(self):
        """Fill the map with default colors"""
        self.drawmapboundary(fill_color=self.color_water)
        self.fillcontinents(color=self.color_land, lake_color=self.color_water)
    
    def _prep_topo_data(self, grid_points=100):
        """Prepare topography data for map
        
        This function prepares high resolution topography for plotting, i.e.
        it reduces the resolution based on the map size
        
        :param int grid_points: number of plotted grid points (default: 100)
        
        """
        
            
        if not isinstance(self.topo_data, TopoData):
            try:
                self.load_topo_data()
            except:
                raise
    
        topo = self.topo_data
        # determine the nu
        pyr_steps = int(np.ceil(np.log2(float(len(topo.lons)) / grid_points)))
        z_max = float(topo.max)
        z_min = float(topo.min)
        z_order = np.floor(np.log10(topo.alt_range))
        X,Y = np.meshgrid(topo.lons, topo.lats)
        x, y = self(X, Y)
        z = topo.data
        if not CV2_AVAILABLE:
            import warnings
            warnings.warn("Could not reduce resolution of topographic data, "
                          "opencv library is not available")
        else:
            from cv2 import pyrDown
            if pyr_steps > 0:
                for k in range(pyr_steps):
                    x = pyrDown(x)
                    y = pyrDown(y)
                    z = pyrDown(z)            
        return (x, y, z, z_min, z_max, z_order)
            
    def draw_topo_contour(self, include_seabed=1, separation_levels=500):
        """Draw topography contour lines
        
         :param bool include_seabed: include seabed topography 
             (default: True)
         :param int separation_levels: separation in m of contour lines
             (default: 500)
        """

        x, y, z, z_min, z_max, z_order = self._prep_topo_data()
        if z_min > 0:
            include_seabed=1
        min_c = int(round(z_min / 10**(z_order - 1)) * 10**(z_order - 1))
        max_c = int(round(z_max / 10**(z_order - 1)) * 10**(z_order - 1))
        if include_seabed:
            levels_contour = np.arange(min_c, max_c, separation_levels)
        else:
            levels_contour = np.arange(0, max_c, separation_levels)
        
        CS1 = self.contour(x, y, z, levels_contour, linewidths=0.5,
                           colors=self.default_colors["contour_lines"])
        CS1.levels = [self._convert_float(val) for val in CS1.levels]
             
        self.ax.clabel(CS1, CS1.levels, inline=True)
        self.contour_lines = CS1   
            
    def _convert_float(self, val):
        """Custom float to string conversion for topo contour plotting
        
        :param float val: the float value to be converted into string        
        """
        s = "%.1f" %val
        if s[-1] == "0":
            return "%.0f" %val
        return "%.1f" %val
    
    def _check_ax3d(self, ax):
        """Check if input is :class:`Axes3D`"""

        from mpl_toolkits import mplot3d
        if isinstance(ax, mplot3d.Axes3D):
            return True
        return False
    
    def set_ticks_topo_colorbar(self, start, stop, step):
        """Update ticks of topo colorbar"""
        cb = self.colorbars["topo"]
        ticks = np.arange(start, stop, step)
        cb.set_ticks(ticks)
        
    def draw_topo(self, insert_colorbar=False, include_seabed=True,
                    max_grid_points=500, cmap_div=None,
                    cmap_seq=None, alpha=0.5, ax=None):
        """Draw topography into map
        
        :param bool insert_colorbar: draws a colorbar for altitude
            range (default: False)
        :param bool include_seabed: include seabed topography 
            (default: True)
        :param int max_grid_points: resolution of displayed topo data 
            points (makes it faster in interactive mode, default: 500)
        :param str cmap_div: name of a diverging colormap (this one is 
            used if :arg:`include_seabed` is True, and the cmap is shifted 
            such , that white colors correspond to sea level altitude, 
            default: "coolwarm")
        :param str cmap_seq: name of a sequential colormap (this one is 
            used if :arg:`include_seabed` is False, default: "Oranges")
        :param float alpha: Alpha value (transparency) of plotted 
            topography
        :param ax: matplotlib axes object
        """
        import matplotlib.cm as colormaps
        if cmap_div is None:
            cmap_div = colormaps.coolwarm
        if cmap_seq is None:
            cmap_seq = colormaps.Oranges
        
        try:  
            if ax is None:
                ax = self.ax
            if ax is None:
                from matplotlib.pyplot import subplots
                fig, ax = subplots(1, 1, figsize=(16,10))
                self.ax = ax
    
            (x, y, z, z_min, z_max, 
             z_order) = self._prep_topo_data(grid_points=max_grid_points)
                
            
            delz = z_max - z_min #in m
            exp = exponent(delz) - 1
            start_alt = round((z_min - delz*0.1)*10**(-exp))*10**exp
            stop_alt = round((z_max + delz*0.1)*10**(-exp))*10**exp
            if z_min > 0:
                include_seabed = 1
                
            if not include_seabed:
                start_alt = 0
            
            z_step = (stop_alt - start_alt) / 1000. 
            levels_filled = np.arange(start_alt, stop_alt, z_step)
            
            if levels_filled[0] < 0:          
                shifted_cmap = shifted_color_map(z_min, z_max, cmap_div)
                
                cs2 = ax.contourf(x, y, z, levels_filled, cmap=shifted_cmap,
                                  extend="both", alpha=alpha)
                                  
            elif levels_filled[0] >= 0:        
                cs2 = ax.contourf(x, y, z, levels_filled, cmap=cmap_seq,
                                  alpha=1.0, extend="min")
                self.contour_filled = cs2
                    
            if insert_colorbar:              
                self.insert_colorbar("topo", cs2, label="Altitude [m]")
                exp = exponent(delz)
                step = round(delz*10**(-exp))*10**(exp-1)
                self.set_ticks_topo_colorbar(start_alt, stop_alt, step)
        except Exception as e:
            raise
            msg=("Could not draw topography in high res, using default "
                 "etopo() instead...")
            print(msg + repr(e))
            self.etopo()
                 
    def draw_topo_3d(self, num_ticks=4, cmap="Oranges", alpha=0.5,
                     contour_color="#708090", contour_antialiased=True, 
                     contour_lw=0.2, ax=None, figsize=(16,8)):
        """Draw topography into 3D axis
        
        Parameters
        ----------        
        num_ticks : int 
            number of lon / lat axis ticks (default: 4)
        cmap : str
            string ID of the colormap used (default: "Oranges")
        alpha : flaot
            Alpha value (transparency) of plotted topography (default: 0.5)
        contour_color : str
            string specifying color of contour lines colors of contour lines 
            (default: "#708090")
        contour_antialiased : bool
            apply antialiasing to surface plot of topography, defaults to False
        contour_lw : 
            width of drawn contour lines, defaults to 0.5, use 0 if you do not 
            want contour lines inserted
            
        ax : 
            axes object
        """
        if not self._check_ax3d(ax):
            if self._check_ax3d(self.ax):
                ax = self.ax
            else:
                from mpl_toolkits import mplot3d
                from matplotlib.pyplot import figure
                fig = figure(figsize=figsize)
                ax = mplot3d.Axes3D(fig)
    
        x, y, z, z_min, z_max, z_order = self._prep_topo_data()
        tickformatter = "{:.2f}"    
        ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cmap, alpha=alpha, 
                        linewidth=contour_lw, antialiased=contour_antialiased, 
                        edgecolors=contour_color, vmin=z_min, vmax=z_max, 
                        zorder=1)
        #xlabels=ax.get_xticklabels()
        ax.set_xlabel("Lons")
        ax.set_ylabel("Lats")
        ax.set_zlabel("Altitude [m]")
        lon_step, lat_step = self.delta_lon * .33, self.delta_lat * .33
        lon_tick_array, lat_tick_array=[], []
        for k in range(num_ticks):
            lon_tick_array.append(self.lon_ll + lon_step * k)
            lat_tick_array.append(self.lat_ll + lat_step * k)
        #lon_tick_array, lat_tick_array = self._prep_coord_ticks()
        lon_coords=[]
        lon_labels=[]
        for lon in lon_tick_array:
            x,_=self(lon,self.lat_ll)
            lon_coords.append(x)
            lon_labels.append(tickformatter.format(lon))
        
        lat_coords=[]
        lat_labels=[]
        for lat in lat_tick_array:
            _, y = self(self.lon_ll, lat)
            lat_coords.append(y)
            lat_labels.append(tickformatter.format(lat))
            
        ax.set_yticks(lat_coords)
        ax.set_yticklabels(lat_labels)
        ax.set_xticks(lon_coords)
        ax.set_xticklabels(lon_labels)
        self.ax = ax
        
    def draw_mapscale_auto(self, **kwargs):
        """Insert a map scale 
        
        Determines missing input parameters automatically and calls 
        :func:`drawmapscale` 
        
        :param kwargs: key word arguments for :func:`drawmapscale` (missing
            ones will be set automatically)
        """
        l = self._len_diag()  
        exp = exponent(l)
        l = round(l * 10**(-exp))*10**(exp)/5.0
        lat_center, lon_center = self.get_map_center()
        str_format = '%d'
        if not "length" in kwargs:
            kwargs["length"] = l
        if l < 40:
            str_format = '%.1f'
# =============================================================================
#         if not "fontsize" in kwargs:
#             kwargs["fontsize"] = 8
# =============================================================================
        if not "lon0" in kwargs:
            kwargs["lon0"] = lon_center
            kwargs["lat0"] = lat_center
        if not "lon" in kwargs:
            kwargs["lon"] = lon_center - self.delta_lon * 0.3
            kwargs["lat"] = lat_center - self.delta_lat * 0.4
        if not "barstyle" in kwargs:
            kwargs["barstyle"] = "fancy"
        if not "units" in kwargs:
            kwargs["units"] = "km"
        if not "format" in kwargs:
            kwargs["format"] = str_format
        self.map_scale = self.drawmapscale(**kwargs)
        
    def remove_map_scale(self):
        """Remove scale"""
        if self.map_scale is not None:
            for i in self.map_scale:
                try:
                    i.remove()
                except:
                    print(("Warning in remove_map_scale: item %s could not be "
                        "removed" %i))
        self.fig.canvas.draw()
    
    def _prep_coord_ticks(self, lon_tick=None, lat_tick=None):
        """Prepare coordinate ticks for parallels and meridians
        
        :param float lon_tick: longitude separation (in decimal degrees).
            Will be set automatically if unspecified
        :param float lat_tick: latitude separation (in decimal degrees)
            Will be set automatically if unspecified
        """
        if lon_tick is None:
            pot_lon = np.floor(np.log10(self.delta_lon))
            lon_tick = np.floor(self.delta_lon / 10**pot_lon) * 10**pot_lon / 4
        if lat_tick is None:
            pot_lat = np.floor(np.log10(self.delta_lat))
            lat_tick = np.floor(self.delta_lat / 10**pot_lat) * 10**pot_lat / 3
            
        lon_tick_array = np.arange(lon_tick * int((self.lon_ll - self.delta_lon 
                                * 0.3) / lon_tick), lon_tick * int((self.lon_tr 
                                + self.delta_lon * 0.3) / lon_tick), lon_tick)
                                
        lat_tick_array = np.arange(lat_tick * int((self.lat_ll - self.delta_lat 
                                * 0.3) / lat_tick), lat_tick * int((self.lat_tr 
                                + self.delta_lat * 0.3) / lat_tick), lat_tick)
                                
        return lon_tick_array, lat_tick_array
    
    def _len_diag(self):
        """Returns the lenght from LL point to TR point in km"""
        return self.haversine(self.lon_ll,self.lat_ll,self.lon_tr,self.lat_tr)
    
        
    def draw_coordinates(self, lat_tick=None, lon_tick=None,
                         labelslon=[0,0,0,1], labelslat=[1,0,0,0], 
                         **kwargs):
        """Draws meridians and parallels 
        
        :param float lat_tick: latitude separation (in decimal degrees)
        :param float lon_tick: longitude separation (in decimal degrees)
        :param labelslon: see basemap docs (default: [0,0,0,1])
        :param labelslat: see basemap docs (default: [1,0,0,0])
        :param **kwargs: additional keyword arguments (passed to 
            :func:`drawmeridians` and :func:`drawparallels`)
        """
        if not "color" in kwargs:
            kwargs["color"] = "gray"
        if not "fmt" in kwargs:
            digs = '%d' %(4 - int(np.floor(np.log10(self._len_diag()))))
            kwargs["fmt"] = "%." + digs + "f"
        lon_tick_array, lat_tick_array = self._prep_coord_ticks(lon_tick,
                                                                lat_tick)
        meridians = self.drawmeridians(lon_tick_array, labels=labelslon,
                                       **kwargs)
        for m in meridians:
            try:
                meridians[m][1][0].set_rotation(25)
            except:
                pass
        self.meridians = meridians
        self.parallels = ps =self.drawparallels(lat_tick_array, 
                                                labels=labelslat, 
                                                **kwargs)
        for p in ps:
            try:
                ps[p][1][0].set_rotation(25)
            except:
                print("Could not rotate parallel...")
            
    def remove_coordinates(self):
        """Remove drawn parallels and meridians"""
        if self.meridians is not None:
            for Dict in [self.meridians, self.parallels]:
                for key in Dict:
                    for item in Dict[key]:
                        try:
                            item[0].remove()
                        except:
                            print("Object %s could not be removed" %item)
        self.fig.canvas.draw()
    
    def legend(self, ax=None, **kwargs):
        """Insert a legend"""
        if ax is None:
            ax = self.ax
# =============================================================================
#         if not "fontsize" in kwargs:
#             kwargs["fontsize"] = 10
# =============================================================================
        ax.legend(loc="best", fancybox=True, framealpha=0.7, **kwargs)
                                                                
    def draw_geo_point_2d(self, p, add_name=False, ax=None, **kwargs):
        """Draw a GeoPoint into 2D basemap
        
        :param GeoPoint p: the actual point        
        :param **kwargs: passed to matplotlib plotting function
            
        """
        if ax is None:
            ax = self.ax
        if not "marker" in kwargs:
            kwargs["marker"] = "^"
        if not "c" in kwargs:
            kwargs["c"]="lime"
        x, y = self(p.longitude, p.latitude) #maps the geo coordinates to figure coordinates
                
        handle = ax.plot(x, y, **kwargs) 
        self.points[p.name] = handle
        if add_name:
            self.write_point_name_2d(p)
        return handle
    
    def draw_geo_vector_2d(self, vec, **kwargs):
        """Draw a :class:`GeoVector3D` into 2D map
        
        :param GeoVector3D vec: the vector 
        
        .. note::
        
            Anchor must be set in the :class:`GeoVector3D` object  
        """
        if not vec.type() == "GeoVector3D":
            raise AttributeError("Wrong input, need :class:`GeoVector3D` "
                                                                "object")
        elif not vec.anchor.type() == "GeoPoint":
            raise AttributeError("Vector anchor not set or wrong type..")

        if not any([x in kwargs for x in["ls", "linestyle"]]):
            kwargs["ls"] = "--"
            
        a = vec.anchor
        pf = a + vec
        x0, y0 = self(a.longitude, a.latitude)
        x1, y1 = self(pf.longitude, pf.latitude)
        return self.draw_line_2d(vec.name, a.latitude, a.longitude,
                                 pf.latitude, pf.longitude, **kwargs)
    
    def add_geo_points_3d(self, pts, distances=None, marker="x", color="b", 
                          connect=True, connect_style="--", alt_offset_m=0):
        """Draws a list of :class:`GeoPoint` objects into the map
        
        :param list pts: geopoint objects
        :param color: color of points (passed to plot function, default: "b")
        :param marker: marker of points (default: "x")
        :param connect: if True, points the points are connected with each 
            other
        :param str connect_style: line style of connection
        """
        xs, ys, zs = [], [], []
        for p in pts:
            try:
                px, py= self(p.lon.decimal_degree,p.lat.decimal_degree)
                xs.append(px), ys.append(py), zs.append(p.altitude)
                self.draw_geo_point_3d(p, alt_offset_m=alt_offset_m, 
                                       marker=marker, s=20, c=color)
            except:
                print("Failed to add %s to map" %p)
        if distances is not None and len(distances) == len(pts):
            sc = self.ax.scatter(xs, ys, np.asarray(zs)+alt_offset_m, zdir='z', 
                                 s=20, c=distances, cmap="Oranges", 
                                 zorder=100000)
            zlabel = "Distance [%km]"
            self.insert_colorbar("dists", sc, label=zlabel)
        elif connect:
            self.ax.plot(xs, ys, np.asarray(zs)+alt_offset_m, 
                         ls=connect_style, 
                         c=color, lw=2, zorder=100000)   
                                                
    def draw_geo_point_3d(self, p, ax=None, alt_offset_m=0.0, **kwargs):
        """Draw a GeoPoint into 3D basemap
        
        :param GeoPoint p: the actual point        
        :param ax: matplotlib axes object (3d axis)
        :param float alt_offset_m: plot point with altitude offset, default 0
        :param **kwargs: passed to matplotlib plotting function
        
        .. note::
        
            The basemap needs to be set up with Axes3d object
            
        """
        if not self._check_ax3d(ax):
            if not self._check_ax3d(self.ax):    
                raise TypeError("No 3D Axes object available...")
            ax=self.ax
        if not "marker" in kwargs:
            kwargs["marker"] = "^"
        if not any([x in kwargs for x in ["c", "color"]]):
            kwargs["c"] = "lime"
        from mpl_toolkits import mplot3d
        if not isinstance(self.ax, mplot3d.Axes3D):
            raise ValueError("Need :class:`Axes3D` object as input...")
        x0, y0 = self(p.longitude, p.latitude)
        z0 = p.altitude + alt_offset_m#*1000
        handle = ax.scatter(x0, y0, z0, zorder=100000, **kwargs)
        self.points[p.name] = handle
        return handle
                                 
    def draw_geo_vector_3d(self, vec, ax=None,  **kwargs):
        """Draw a :class:`GeoVector3D` into 3D map
        
        :param GeoVector3D vec: the vector 
        
        .. note::
        
            Anchor must be set in the :class:`GeoVector3D` object  
            
        """ 
        if ax is None:
            ax = self.ax
        try:
            from mpl_toolkits import mplot3d
            if not isinstance(ax, mplot3d.Axes3D):
                raise ValueError("Need :class:`Axes3D` object as input...")
            elif not vec.type() == "GeoVector3D":
                raise AttributeError("Wrong input, need :class:`GeoVector3D` "
                                                                    "object")
            elif not vec.anchor.type() == "GeoPoint":
                raise AttributeError("Vector anchor not set or wrong type..")
    
            if not any([x in kwargs for x in["ls", "linestyle"]]):
                kwargs["ls"] = "--"
            
            a = vec.anchor
            pf = a + vec
            x0, y0 = self(a.longitude, a.latitude)
            z0 = a.altitude#*1000
            x1, y1 = self(pf.longitude, pf.latitude)
            z1 = pf.altitude#*1000
    
            l = mplot3d.art3d.Line3D((x0,x1), (y0,y1), (z0,z1), **kwargs)
            handle = ax.add_line(l)
            self.lines[vec.name] = handle
            return pf

        except:
            print(vec.anchor.type())
            raise
    
    def add_polygon_2d(self, points=[], poly_id="undefined", ax=None, **kwargs):
        """Add a polygon specified by list of input points
        
        :param list points: list with :class:`GeoPoint` objects
        :param str poly_id: string ID of this object (e.g. for 
            deletion, default: "undefined")
        """
        if ax is None:
            ax = self.ax
        if not "label" in kwargs:
            kwargs["label"] = poly_id
        coords=[]
        for p in points:
            try:
                coords.append(self(p.longitude, p.latitude))
            except Exception as e:
                print("Failed to add one point to poly: " + repr(e))
        from matplotlib.pyplot import Polygon
        polygon = Polygon(coords, **kwargs)
        ax.add_patch(polygon)
        
    def add_polygon_3d(self, points=[], poly_id="undefined", ax=None, **kwargs):
        """Add a polygon specified by list of input points 
        
        :param list points: list with :class:`GeoPoint` objects
        :param str poly_id: string ID of this object (e.g. for 
            deletion, default: "undefined")
        """
        from mpl_toolkits import mplot3d
        if ax is None:
            ax = self.ax
        if not "label" in kwargs:
            kwargs["label"] = poly_id
        xs, ys, zs = [], [], []
        for p in points:
            x,y=self(p.longitude, p.latitude)
            xs.append(x)
            ys.append(y)
            zs.append(p.altitude)#*1000)
        coords = [list(zip(xs, ys, zs))]
        poly_coll = mplot3d.art3d.Poly3DCollection(coords, **kwargs)
        ax.add_collection3d(poly_coll)
           
    def draw_line_2d(self,line_id, lat0, lon0, lat1, lon1, **kwargs):
        """Draw a line between 2 geo coordinates
        
        :param str line_id: ID of the line artist (for deletion management)
        :param float lon0: start longitude 
        :param float lat0: start latitude 
        :param float lon1: stop longitude 
        :param float lat1: stop latitude
        :returns: line object
        """
        line = self.drawgreatcircle(lon0, lat0, lon1, lat1, **kwargs)
        self.lines[line_id] = line
        return line
    
    def write_point_name_2d(self, p, dist=0, angle=0, ax=None, **kwargs):
        """Annotate name of point to point in map
        
        :param GeoPoint p: the point
        :param float dist: distance of annotation in km (default: 0)
        :param float angle: angular direction of  diplacement (default: 0)
        :param ax: matplotlib axes instance
        """
        if dist == 0:
            self.draw_text_2d(p.longitude, p.latitude, p.name, ax, **kwargs) 
        else:
            pt = p.offset(azimuth=angle, dist_hor=dist)
            self.annotate_text_2d(pt.longitude, pt.latitude, p.name,
                                  p.longitude, p.latitude, ax, **kwargs)
    
    def annotate_text_2d(self, lon_text, lat_text, text, lon_point, lat_point,
                         ax=None, **kwargs):
        """Annotate text to a certain position in a 2D basemap
        
        :param float lon_text: longitude coordinate of text position
        :param float lat_text: latitude coordinate of text position
        :param str text: the actual text
        :param float lon_point: longitude of text annotation
        :param float lat_point: latitude of text annotation
        :param **kwargs: additional keyword arguments passed to 
            :func:`annotate`
        """
        if ax is None:
            ax = self.ax
        x_text, y_text = self(lon_text, lat_text)  
        x_point, y_point = self(lon_point, lat_point)
        t = ax.annotate(text, xy=(x_point, y_point), xytext=(x_text, y_text), 
                        arrowprops=dict(arrowstyle="->", connectionstyle=
                        "arc,angleA=10,armA=20,rad=6", shrinkA=2.0, 
                        shrinkB=10), **kwargs)
                        
        self.texts[text] = t
        return t
            
    def draw_text_2d(self, lon, lat, text, ax=None, **kwargs):
        """Draw a text into a 2D map
                
        :param float lon: longitude of text
        :param float lat: latitude of text
        :param str text: the actual text
        :param **kwargs: draw parameters 
        """
        if ax is None:
            ax = self.ax
        x, y = self(lon, lat)  
        t = ax.text(x, y, text, **kwargs)
        self.texts[text] = t
        return t
    
    def draw_text_3d(self, lon, lat, alt, text, ax=None, **kwargs):
        """Draw a text into a 3D map
        
        :param float lon: longitude of text
        :param float lat: latitude of text
        :param float alt: altitude of text (in m)
        :param str text: the actual text
        :param **kwargs: draw parameters 
        """
        if not self._check_ax3d(ax):
            if not self._check_ax3d(self.ax):    
                raise TypeError("No 3D Axes object available...")
            ax = self.ax
        x, y = self(lon, lat)  
        t = ax.text(x, y, alt, text, zorder=1e6, **kwargs)
        self.texts[text] = t
        return t
    
    def draw_data_scatter_2d(self, data=None, lons=None, lats=None,
                             unit=None, data_id="n/d", ax=None, **kwargs):
        """Draw data into map
        
        If input is unspecified, a random test data set will be created and
        plotted onto the map.
        
        :param float data: the actual data
        :param float lons: array with longitude coordinates
        :param float lats: array with latitude coordinates
        :param str unit: Unit of Data for labelling
        :param str data_id: string ID under which data will be stored
        :param ax: axes object
        :param **kwargs: keyword arguments passed to :func:`scatter`
        
        """
        if ax is None:
            ax = self.ax
        if not "marker" in kwargs:
            kwargs["marker"] = "^"
        if not "s" in kwargs:
            kwargs["s"] = 30
        if not any([x in kwargs for x in ["edgecolor", "edgecolors"]]):
            kwargs["edgecolor"] = "none"
        
        if data is None:
            print("No input data, create random data")
            lons, lats, data = self.make_random_data()
            data_id = "TestData"
            unit = 'n/d'
    
        x, y = self(lons, lats)
        z = np.array(data)
        sc = ax.scatter(x, y, c=z, **kwargs)
        zlabel = "%s [%s]" %(data_id, unit)
        self.insert_colorbar(data_id, sc, label=zlabel)
        self.plotted_data[data_id] = sc

    def insert_colorbar(self, obj_id, obj, label=None, ax=None, **kwargs):
        """Insert a colorbar into the map
        
        :param obj_id: save ID
        :param obj: the actual artist to which the colorbar belongs
        :param str label: colorbar label
        :param **kwargs: further arguments for :func:`colorbar` call 
            (see matplotlib docs)
        """
        if ax is None:
            ax = self.ax
        cb = self.fig.colorbar(obj, shrink=0.9, ax=ax, **kwargs)
        if isinstance(label, str):
            cb.set_label(label, fontsize=16)
        self.colorbars[obj_id] = cb
        ax.figure.canvas.draw()
        
    def remove_colorbar(self, cb_id):
        """Remove a drawn colorbar
        
        :param str cb_id: string ID of colorbar        
        """
        if not cb_id in self.colorbars:
            raise KeyError("No colorbar found with ID: " + cb_id)
        cb = self.colorbars[cb_id]
        cb.remove()
        del self.colorbars[cb_id]
        self.fig.canvas.draw()
    
    @property
    def set_ax(self):
        """Check current axes and write into ``self.ax``"""
        self.ax = self._check_ax()
        return self.ax  
        
    @property
    def fig(self):
        """The current figure canvas"""
        return self.set_ax.figure
        
    def remove_all_colorbars(self):
        """Remove all colorbars from map"""
        cbs = self.colorbars
        for key, cb in cbs.items():
            cb.remove()
            del self.colorbars[key]
        self.colorbars = {}
        self.fig.canvas.draw()
        
    """Calculations etc
    """    
    
    def haversine(self, lon0, lat0, lon1, lat1):
        """Haversine formula (Distance between two geo coordinates)
        
        Approximate horizontal distance between 2 points assuming a spherical 
        earth
        
        :param float lon0: longitude of first point in decimal degrees
        :param float lat0: latitude of first point in decimal degrees
        :param float lon1: longitude of second point in decimal degrees
        :param float lat1: latitude of second point in decimal degrees
        """
        return haversine_formula(lon0, lat0, lon1, lat1)
            
    def get_map_center(self):
        """Returns center coordinates of map (lat, lon)"""
        lon_center = self.lon_ll + self.delta_lon / 2.0
        lat_center = self.lat_ll + self.delta_lat / 2.0
        return (lat_center, lon_center)
        
    """Other stuff"""
    def make_random_data(self,total_number=200, value_range=[-100, 100, 1]):
        """Create random data in the current regime"""
        from random import randrange
        lons, lats = np.zeros(total_number), np.zeros(total_number)
        data = np.zeros(total_number)
        for k in range(total_number):
            lons[k] = 0.1 * randrange(self.lon_ll * 10, self.lon_tr * 10, 1)
            lats[k] = 0.1 * randrange(self.lat_ll * 10, self.lat_tr * 10, 1)
            data[k] = randrange(value_range[0], value_range[1], value_range[2])
            
        return lons, lats, data
    
    """I/O, printing, etc...
    """
    def print_map_info(self):
        """Print some basic info about this map"""
        lat_center,lon_center=self.get_map_center()
        print()
        print("-------------------------------")
        print("------- MAP INFORMATION -------")
        print("-------------------------------")
        print()
        print("(Lon|Lat) lower left: (%s|%s)" %(self.lon_ll, self.lat_ll))
        print("(Lon|Lat) top right: (%s|%s)" %(self.lon_tr, self.lat_tr))
        print("(Lon|Lat) center: (%s|%s)" %(lon_center, lat_center))
        print("Projection: %s" %self.projection)
        
    def save_as(self, save_dir=None, save_name=None, ftype="png"):
        """Save current figure as image"""
        import os
        if save_dir is None:
            save_dir = os.getcwd()
        if save_name is None:
            save_name = "geonum_map"
        self.fig.savefig(os.path.join(save_dir, save_name + ftype))
        
if __name__ == "__main__":
    from matplotlib.pyplot import close
    from geonum import GeoSetup
    close("all")
    s = GeoSetup()
    s.create_test_data()
    s.plot_2d()
    