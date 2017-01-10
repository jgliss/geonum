"""geonum example script 4

Elevation profiles and directional topographic analysis
"""
import geonum
from matplotlib.pyplot import show

def load_topo_oslo():
    dat = geonum.topodata.TopoDataAccess()
    topodata = dat.get_data(59.8728, 10.375, 60.026, 10.9144)
    return topodata
    
    
if __name__ == "__main__":
    topo_data = load_topo_oslo()
    basemap = topo_data.plot_2d()
    basemap.ax.set_title("SRTM topo data Oslo")
    basemap.draw_mapscale_auto()
    show()
    