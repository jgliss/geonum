# -*- coding: utf-8 -*-
"""
geonum example script 3

Mapping and advanced functionality
"""

import geonum
from numpy import mean
from matplotlib.pyplot import close, subplots, show
from os.path import join
from os import getcwd

### Set save directory for figures
save_path = join(getcwd(), "scripts_out")

close("all")

def create_map_and_load_topodata():
    """Create and return map object around Guallatiri volcano"""
    lat0 = -18.48
    lon0 = -69.14
    lat1 = -18.39  
    lon1 = -69.06
    
    lat_center = mean([lat0, lat1])
    lon_center = mean([lon0, lon1])
    
    m = geonum.mapping.Map(projection = "lcc", llcrnrlat = lat0, llcrnrlon =\
            lon0, urcrnrlat = lat1, urcrnrlon = lon1, lat_0 = lat_center,\
            lon_0 = lon_center)
    
    # Compared to a normal Basemap object, a geonum Map can directly access 
    # SRTM topo data
    m.load_topo_data("srtm")
    
    return m

def add_points_and_plot_map(basemap):
    """Add some points and plot in 2d and 3d"""
    basemap.draw_topo_3d()
    
    # create a geopoint for the summit region (here, altitude is set manually)
    summit = geonum.GeoPoint(-18.423672, -69.090369, 6071.0, name = "summit")
    
    # draw this point and a text into the map
    basemap.draw_geo_point_3d(summit)
    basemap.draw_text_3d(-69.090369, -18.423672, 6100.0,\
                "Guallatiri summit", color = "k")
    basemap.ax.set_title("Guallatiri volcano")
    
    # create two geopoints, one for a tough guy aiming to save the world from 2 
    # threats, the other one one of the 2 threats, probably the worse one
    # (without altitude spec -> is retrieved automatically from topo data)
    p1 = geonum.GeoPoint(-18.47, -69.12, name = "Observer")
    p2 = geonum.GeoPoint(-18.40, -69.12, name = "Llamas")
    
    
    # points can also be drawn directly into the map (here including the name which
    # is placed 50 m above the topographic altitude of Scientist)
    p1.plot_3d(basemap, add_name = True, dz_text = 50)
    p2.plot_3d(basemap, add_name = True, dz_text = 50)
    
    # You can include a polygon connecting different points
    basemap.add_polygon_3d([summit, p2, p1])
    
    # just some little adjustment of the 3D viewing direction
    basemap.ax.view_init(20, 240)
    
    # save the 3D figure
    basemap.ax.figure.savefig(join(save_path, "ex3_out_1_map3D.png"))

    fig, ax = subplots(1,1)
    
    # update the axes object in the map
    basemap.ax = ax
    basemap.draw_topo() # draws 2D map
    basemap.draw_coordinates()
    p1.plot_2d(basemap, add_name = True)
    p2.plot_2d(basemap, add_name = True)
    
    summit.plot_2d(basemap, add_name = True, dist_text = 1.0,alpha = 0.2,\
                                                       angle_text = 170)
    
    # You can include a polygon connecting different points
    basemap.add_polygon_2d([summit, p2, p1], fc = "r", alpha = 0.2)
    basemap.ax.figure.savefig(join(save_path, "ex3_out_2_map2D.png"))

    show()

if __name__ == "__main__":
    m = create_map_and_load_topodata()
    add_points_and_plot_map(m)


