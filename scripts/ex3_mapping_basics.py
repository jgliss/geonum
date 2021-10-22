# -*- coding: utf-8 -*-
"""
geonum example script 3

Creating basemaps including topographic data and some further features
(i.e. draw stuff into a map)
"""

import geonum
from numpy import mean
from matplotlib.pyplot import close, subplots

def create_map_and_load_topodata():
    """Create and return map object around Guallatiri volcano"""
    # coordinate range of map
    lat0 = -18.48
    lon0 = -69.15
    lat1 = -18.39
    lon1 = -69.05

    # trivial ...
    lat_center = mean([lat0, lat1])
    lon_center = mean([lon0, lon1])

    # Create a map object - the Map class is initialised as Basemap object
    # and is extended by some features
    m = geonum.mapping.Map(projection="lcc", llcrnrlat=lat0,
                           llcrnrlon=lon0, urcrnrlat=lat1, urcrnrlon=lon1,
                           lat_0=lat_center, lon_0=lon_center)

    # Compared to a normal Basemap object, a geonum Map can directly access
    # and include SRTM topographic data
    m.load_topo_data("srtm")

    return m

def add_points_and_plot_map(basemap):
    """Add some points and plot in 2d and 3d"""
    basemap.draw_topo_3d()

    # create a geopoint for Guallatiri summit region
    # (here, altitude is set manually)
    summit = geonum.GeoPoint(latitude=-18.423672, longitude=-69.090369,
                             altitude=6071.0, name="Guallatiri")

    # draw this point and a text into the map
    basemap.draw_geo_point_3d(summit)
    basemap.draw_text_3d(lon=-69.090369, lat=-18.423672, alt=6100.0,
                         text="Guallatiri", color="k")
    basemap.ax.set_title("Overview map Guallatiri volcano")

    # create two more objects (without specifying altitude -> is retrieved
    # automatically from topo data)
    p1 = geonum.GeoPoint(-18.45, -69.12, name="Tourist (breathing hard)",
                         auto_topo_access=True)
    p2 = geonum.GeoPoint(-18.40, -69.12, name="Llamas",
                         auto_topo_access=True)


    # points can also be drawn directly into the map (here including the
    # name which is placed 50 m above the topographic altitude of Tourist)
    p1.plot_3d(basemap, add_name=True, dz_text=50)
    p2.plot_3d(basemap, add_name=True, dz_text=50)

    # You can include a polygon connecting different points
    basemap.add_polygon_3d([summit, p2, p1], color="lime")

    # just some little adjustment of the 3D viewing direction
    basemap.ax.view_init(20, 240)

    fig, ax = subplots(1,1)

    # update the axes object in the map
    basemap.ax = ax
    basemap.draw_topo() # draws topography into 2D map
    basemap.draw_topo_contour()

    basemap.draw_coordinates()
    basemap.draw_mapscale_auto()
    p1.plot_2d(basemap, add_name=True)
    p2.plot_2d(basemap, add_name=True)

    summit.plot_2d(basemap, add_name=True, dist_text=1.0, alpha=0.2,
                   angle_text=170, c="y")

    # You can include a polygon connecting different points (for whatever
    # reason)
    basemap.add_polygon_2d([summit, p2, p1], fc="lime", alpha=0.2)

if __name__ == "__main__":
    close('all')
    if not geonum.BASEMAP_AVAILABLE:
        print('Cannot run example script 3. Basemap is not installed')
    else:
        m = create_map_and_load_topodata()
        add_points_and_plot_map(m)