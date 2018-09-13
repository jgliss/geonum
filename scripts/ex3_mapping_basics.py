# -*- coding: utf-8 -*-
"""
geonum example script 3

Creating basemaps including topographic data and some further features 
(i.e. draw stuff into a map)
"""

import geonum
from numpy import mean
from matplotlib.pyplot import close, subplots, show
from os.path import join, exists
from os import getcwd
from SETTINGS import OPTPARSE
from numpy import testing as npt

### Set save directory for figures
save_path = join(getcwd(), "scripts_out")

close("all")
exceptions = []

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
    summit = geonum.GeoPoint(lat=-18.423672, lon=-69.090369,
                             altitude=6071.0, name="Guallatiri")
    
    # draw this point and a text into the map
    basemap.draw_geo_point_3d(summit)
    basemap.draw_text_3d(lon=-69.090369, lat=-18.423672, alt=6100.0,
                         text="Guallatiri", color="k")
    basemap.ax.set_title("Overview map Guallatiri volcano")
    
    # create two more objects (without specifying altitude -> is retrieved 
    # automatically from topo data)
    p1 = geonum.GeoPoint(-18.45, -69.12, name="Tourist (breathing hard)")
    p2 = geonum.GeoPoint(-18.40, -69.12, name="Llamas")
    
    
    # points can also be drawn directly into the map (here including the 
    # name which is placed 50 m above the topographic altitude of Tourist)
    p1.plot_3d(basemap, add_name=True, dz_text=50)
    p2.plot_3d(basemap, add_name=True, dz_text=50)
    
    # You can include a polygon connecting different points
    basemap.add_polygon_3d([summit, p2, p1], color="lime")
    
    # just some little adjustment of the 3D viewing direction
    basemap.ax.view_init(20, 240)
    
    
    # save the 3D figure
    print(("Save path: %s (exists: %s)" %(save_path, exists(save_path))))
    
    try:
        basemap.ax.figure.savefig(join(save_path, "ex3_out_1_map3D.png"))
    except Exception as e:
        exceptions.append("Failed to save first plot...%s" %repr(e))
        

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
    try:
        basemap.ax.figure.savefig(join(save_path, "ex3_out_2_map2D.png"))
    except:
        exceptions.append("Failed to save second plot... %s" %repr(e))
    for e in exceptions:
        print(e)

if __name__ == "__main__":
    m = create_map_and_load_topodata()
    add_points_and_plot_map(m)
    
    # Import script options
    (options, args) = OPTPARSE.parse_args()
    
    # If applicable, do some tests. This is done only if TESTMODE is active: 
    # testmode can be activated globally (see SETTINGS.py) or can also be 
    # activated from the command line when executing the script using the 
    # option --test 1
    if int(options.test):
        from os.path import basename
        npt.assert_array_equal([m.topo_data.data.shape],
                               [(110, 122)])
        
        
        actual = [m.topo_data.data.mean(),
                  m.delta_lat,
                  m.delta_lon]
        
        actual.extend(m.get_map_center())
        npt.assert_allclose(actual=actual,
                            desired=[4982.759463487,
                                     0.089999999999,
                                     0.100000000000,
                                     -18.43500000000,
                                     -69.10000000000],
                            rtol=1e-7)
        print(("All tests passed in script: %s" %basename(__file__))) 
    try:
        if int(options.show) == 1:
            show()
    except:
        print("Use option --show 1 if you want the plots to be displayed")


