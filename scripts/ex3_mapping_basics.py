# -*- coding: utf-8 -*-
"""
geonum example script 3

Mapping and advanced functionality
"""

import geonum
from numpy import mean
from matplotlib.pyplot import close, subplots, show

close("all")

lat0 = -18.48
lon0 = -69.14
lat1 = -18.39  
lon1 = -69.06

lat_center = mean([lat0, lat1])
lon_center = mean([lon0, lon1])

m = geonum.mapping.Map(projection = "lcc", llcrnrlat = lat0, llcrnrlon = lon0,\
        urcrnrlat = lat1, urcrnrlon = lon1, lat_0 = lat_center,\
        lon_0 = lon_center)

# Compared to a normal Basemap object, a geonum Map can directly access 
# SRTM topo data
m.load_topo_data("srtm")

m.draw_topo_3d()

# create a geopoint for the summit region (here, altitude is set manually)
summit = geonum.GeoPoint(-18.423672, -69.090369, 6071.0, name = "Mordor (Threat 1)")

# draw this point and a text into the map
m.draw_geo_point_3d(summit)
m.draw_text_3d(-69.090369, -18.423672, 6100.0, "Mordor (Threat 1)", color="k")
m.ax.set_title("Guallatiri volcano")

# create two geopoints, one for a tough guy aiming to save the world from 2 
# threats, the other one one of the 2 threats, probably the worse one
# (without altitude spec -> is retrieved automatically from topo data)
frodo = geonum.GeoPoint(-18.47, -69.12, name = "Frodo (The cure)")
trump = geonum.GeoPoint(-18.40, -69.12, name = "D. T. (Threat 2)")


# points can also be drawn directly into the map (here including the name which
# is placed 50 m above the topographic altitude of Frodo)
frodo.plot_3d(m, add_name = True, dz_text = 50)
trump.plot_3d(m, add_name = True, dz_text = 50)

# You can include a polygon connecting different points
m.add_polygon_3d([summit, trump, frodo])

# just some little adjustment of the 3D viewing direction
m.ax.view_init(20, 240)

fig, ax = subplots(1,1)
# update the axes object in the map
m.ax = ax
m.draw_topo()
m.draw_coordinates()
frodo.plot_2d(m, add_name = True)
trump.plot_2d(m, add_name = True)

summit.plot_2d(m, add_name = True, dist_text = 1.0,alpha = 0.2,  angle_text = 170)
# You can include a polygon connecting different points
m.add_polygon_2d([summit, trump, frodo], fc = "r", alpha = 0.2)
##m.ax.text(m())

show()



