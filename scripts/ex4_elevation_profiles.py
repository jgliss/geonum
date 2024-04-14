"""geonum example script 4

Elevation profiles and directional topographic analysis
"""

from matplotlib.pyplot import subplots

### Imports from other example scripts
from ex2_geosetup_intro import create_geosetup

if __name__ == "__main__":
    # create Etna GeoSetup from example 2
    s = create_geosetup()

    # get camera
    cam = s.points["cam"]
    print(cam)

    # calculate elevation profile between camera and source (Etna summit)
    # the camera altitude is 15 m above topography
    profile = cam.get_elevation_profile(geo_point=s.points["source"])

    (dist,
     dist_err,
     intersect,
     view_elevations,
     ax)= profile.get_first_intersection(elev_angle=10.0,
                                         view_above_topo_m=15.0,
                                         plot=True)

    # Find the elvation angle corresponding to the horizon from the
    # elevation profile
    elev, elev_secs, dist_secs = profile.find_horizon_elev(elev_start=10.0,
                                                           elev_stop=20.0,
                                                           step_deg=0.1,
                                                           view_above_topo_m=15.0)

    fig, ax = subplots(1,1)
    ax.plot(elev_secs, dist_secs, "--x")
    ax.set_xlabel("Elevation angle [deg]")
    ax.set_ylabel("Retrieved distances [km]")
    ax.set_title("Elev horizon: %s" %elev)