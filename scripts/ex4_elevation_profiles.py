"""geonum example script 4

Elevation profiles and directional topographic analysis
"""

from matplotlib.pyplot import show, subplots
from os.path import join
from os import getcwd
from SETTINGS import OPTPARSE
from numpy import testing as npt

### Imports from other example scripts
from ex2_geosetup_intro import create_geosetup

### Set save directory for figures
save_path = join(getcwd(), "scripts_out")

def calc_profile():
    """Not very explanatory function name..."""


    # create Etna GeoSetup from example 2
    s = create_geosetup()

    # get camera
    cam = s.points["cam"]
    print(cam)

    # calculate elevation profile between camera and source (Etna summit)
    # the camera altitude is 15 m above topography
    elev_profile = cam.get_elevation_profile(geo_point=s.points["source"])


    return elev_profile

if __name__ == "__main__":
    profile = calc_profile()

    (dist,
     dist_err,
     intersect,
     view_elevations,
     ax)= profile.get_first_intersection(elev_angle=10.0,
                                         view_above_topo_m=15.0,
                                         plot=True)

    ax.figure.savefig(join(save_path,
                           "ex4_out_1_elev_profile_intersect.png"))

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
    fig.savefig(join(save_path, "ex4_out_2_horizon_search.png"))

    # Import script options
    (options, args) = OPTPARSE.parse_args()

    # If applicable, do some tests. This is done only if TESTMODE is active:
    # testmode can be activated globally (see SETTINGS.py) or can also be
    # activated from the command line when executing the script using the
    # option --test 1
    if int(options.test):
        from os.path import basename
        npt.assert_array_equal([], [])


        actual = [profile.start_point.longitude,
                  profile.start_point.latitude,
                  profile.start_point.altitude,
                  profile.alt_range,
                  dist, dist_err,
                  intersect.longitude,
                  intersect.latitude,
                  intersect.altitude,
                  elev]

        npt.assert_allclose(actual=actual,
                            desired=[15.1129, 37.73122, 803.,
                                     2482.5040909723493,
                                     7.667649e+00,
                                     8.378121e-03,
                                     1.502820e+01,
                                     3.774504e+01,
                                     2.169563e+03,
                                     13.299999999999988],
                            rtol=1e-6)
        print(("All tests passed in script: %s" %basename(__file__)))
    try:
        if int(options.show) == 1:
            show()
    except:
        print("Use option --show 1 if you want the plots to be displayed")
