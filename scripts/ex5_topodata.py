"""geonum example script 4

Elevation profiles and directional topographic analysis
"""
import geonum
from matplotlib.pyplot import show
from os.path import join
from os import getcwd
### Set save directory for figures
save_path = join(getcwd(), "scripts_out")

def load_topo_oslo():
    dat = geonum.topodata.TopoDataAccess()
    topodata = dat.get_data(59.8728, 10.375, 60.026, 10.9144)
    return topodata
    
    
if __name__ == "__main__":
    topo_data = load_topo_oslo()
    basemap = topo_data.plot_2d()
    basemap.ax.set_title("SRTM topo data Oslo")
    basemap.draw_mapscale_auto()
    
    my_flat = geonum.base.GeoPoint(59.919386, 10.714970, name =\
                                            "the author lives here")
    my_flat.plot_2d(basemap, add_name=True)
    basemap.ax.figure.savefig(join(save_path, "ex5_out_1_oslo_map.png"))
    show()
    