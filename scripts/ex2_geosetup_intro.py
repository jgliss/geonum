"""geonum example script 2

Introduction into the GeoSetup class
"""

from geonum import GeoPoint, GeoVector3D, GeoSetup, BASEMAP_AVAILABLE
from matplotlib.pyplot import close

def create_geosetup():
    s = GeoSetup()

    # Create two GeoPoints called source and cam, for both points
    # retrieve topographic altitude via "auto_topo_access"
    source = GeoPoint(37.751005,  14.993435, name="source",
                      auto_topo_access=True)

    camera = GeoPoint(37.73122,  15.1129, name="cam", auto_topo_access=True)

    # Add the two GeoPoints to the GeoSetup
    s.add_geo_points(source, camera)
    s.set_borders_from_points()

    # Create vector anchored at source pointing south with horizontal
    # orientation (elevation angle 0)
    plume = GeoVector3D(azimuth=180, dist_hor=s.magnitude,
                        elevation=0, anchor=source, name="plume")

    # Create viewing direction vector anchored at cam pointing west
    # at elevation angle of 8deg
    view_dir = GeoVector3D(azimuth=270, dist_hor=s.magnitude,
                           elevation=8, anchor=camera, name="cam_view")

    # Add the two GeoVectors to the GeoSetup class
    s.add_geo_vectors(plume, view_dir)

    # return the GeoSetup
    return s

def plot_geosetup(geosetup):
    # Plot 2D and 3D overview maps of GeoSetup
    map2d = geosetup.plot_2d()
    map3d = geosetup.plot_3d()

    return map2d, map3d

if __name__ == "__main__":
    close("all")
    # see function definition for example code
    s = create_geosetup()

    if BASEMAP_AVAILABLE:
        map2d, map3d = plot_geosetup(s)
