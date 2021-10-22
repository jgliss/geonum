"""geonum example script 5

Low level access to topographic data (here SRTM) using the example of
Oslo (Norway).

Note
----
Very short example, does not provide the deepest insights but gives an idea

"""
import geonum

if __name__ == "__main__":
    topo_access = geonum.TopoDataAccess(mode='srtm')
    topodata = topo_access.get_data(59.8728, 10.375, 60.026, 10.9144)
    my_old_flat = geonum.GeoPoint(59.919386, 10.714970,
                                  name="Somewhere in Oslo",
                                  auto_topo_access=True)

    if geonum.BASEMAP_AVAILABLE:
        basemap = topodata.plot_2d()
        basemap.ax.set_title("SRTM topo data Oslo")
        basemap.draw_mapscale_auto()


        my_old_flat.plot_2d(basemap, add_name=True)