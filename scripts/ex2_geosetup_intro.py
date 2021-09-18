"""geonum example script 2

Introduction into the GeoSetup class
"""

from geonum import GeoPoint, GeoVector3D, GeoSetup, BASEMAP_AVAILABLE
from matplotlib.pyplot import show, close, rcParams
from os.path import join
from os import getcwd
from SETTINGS import OPTPARSE
from numpy import testing as npt
### Set save directory for figures
save_path = join(getcwd(), "scripts_out")

rcParams["font.size"] = 12
def create_geosetup():
    s = GeoSetup()

    #Create two GeoPoints for source and instrument
    source = GeoPoint(37.751005,  14.993435, name="source", auto_topo_access=True)
    camera = GeoPoint(37.73122,  15.1129, name="cam", auto_topo_access=True)

    # Add the two GeoPoints to the GeoSetup
    s.add_geo_points(source, camera)
    s.set_borders_from_points()

    # Create plume vector anchored at source pointing south with horizontal
    # orientation (elevation angle 0)
    plume = GeoVector3D(azimuth=180, dist_hor=s.magnitude,
                        elevation=0, anchor=source, name="plume")
    # Create viewing direction vector anchored at instrument pointing west
    # at elevation angle of 8deg
    view_dir = GeoVector3D(azimuth=270, dist_hor=s.magnitude,
                           elevation=8, anchor=camera, name="cam_view")

    # Add the two GeoVectors to the GeoSetup class
    s.add_geo_vectors(plume, view_dir)
    return s

def plot_geosetup(geosetup):
    # Now plot 2D and 3D overview maps
    if not BASEMAP_AVAILABLE:
        raise ImportError('Basemap module is not available, cannot plot map')
    map2d = geosetup.plot_2d()
    map3d = geosetup.plot_3d()

    return map2d, map3d

if __name__ == "__main__":
    close("all")
    s = create_geosetup()
    try:
        map2d, map3d = plot_geosetup(s)
    except ImportError as e:
        print(repr(e))
    else:
        map2d.ax.figure.savefig(join(save_path, "ex2_out_1_map2D.png"))
        map3d.ax.figure.savefig(join(save_path, "ex2_out_2_map3D.png"))

    # Import script options
    (options, args) = OPTPARSE.parse_args()

    # If applicable, do some tests. This is done only if TESTMODE is active:
    # testmode can be activated globally (see SETTINGS.py) or can also be
    # activated from the command line when executing the script using the
    # option --test 1
    if int(options.test):
        from os.path import basename
        npt.assert_array_equal([4, 2, True],
                               [len(s.points),
                                len(s.vectors),
                                all([x in ['source',
                                           'll',
                                           'tr',
                                           'cam'] for x in list(s.points)])],
                                all([x in ['plume',
                                           'cam_view'] for x in list(s.vectors)]))

        actual = [s.points['source'].altitude,
                  s.points['cam'].altitude]
        npt.assert_allclose(actual=actual,
                            desired=[3264.0,
                                     803.0],
                            rtol=1e-7)
        print(("All tests passed in script: %s" %basename(__file__)))
    try:
        if int(options.show) == 1:
            show()
    except:
        print("Use option --show 1 if you want the plots to be displayed")

