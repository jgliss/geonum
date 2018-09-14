"""geonum example script 5

Low level access to topographic data (here SRTM) using the example of 
Oslo (Norway).

Note
----
Very short example, does not provide the deepest insights but gives an idea

"""
import geonum
from matplotlib.pyplot import show
from os.path import join
from os import getcwd
from SETTINGS import OPTPARSE
from numpy import testing as npt
### Set save directory for figures
save_path = join(getcwd(), "scripts_out")

def load_topo_oslo():
    dat = geonum.topodata.TopoDataAccess(mode='srtm')
    topodata = dat.get_data(59.8728, 10.375, 60.026, 10.9144)
    return topodata
    
    
if __name__ == "__main__":
    topo_data = load_topo_oslo()
    my_flat = geonum.base.GeoPoint(59.919386, 10.714970,
                                   name="Somewhere in Frogner (Oslo)")
    
    if geonum.BASEMAP_AVAILABLE:
        basemap = topo_data.plot_2d()
        basemap.ax.set_title("SRTM topo data Oslo")
        basemap.draw_mapscale_auto()
        
        
        my_flat.plot_2d(basemap, add_name=True)
        basemap.ax.figure.savefig(join(save_path, "ex5_out_1_oslo_map.png"))
        
    # Import script options
    (options, args) = OPTPARSE.parse_args()
    
    # If applicable, do some tests. This is done only if TESTMODE is active: 
    # testmode can be activated globally (see SETTINGS.py) or can also be 
    # activated from the command line when executing the script using the 
    # option --test 1
    if int(options.test):
        from os.path import basename
        npt.assert_array_equal([], [])
        
        actual = [my_flat.altitude]
        npt.assert_allclose(actual=actual,
                            desired=[44.0],
                            rtol=1e-6)
        print(("All tests passed in script: %s" %basename(__file__))) 
    try:
        if int(options.show) == 1:
            show()
    except:
        print("Use option --show 1 if you want the plots to be displayed")
    
    