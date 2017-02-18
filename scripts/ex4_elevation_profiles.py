"""geonum example script 4

Elevation profiles and directional topographic analysis
"""

from matplotlib.pyplot import show, subplots
from os.path import join
from os import getcwd

### Imports from other example scripts
from ex2_geosetup_intro import create_geosetup

### Set save directory for figures
save_path = join(getcwd(), "scripts_out")

def run_ex4():
    """Not very explanatory function name..."""
    
    
    # create Etna GeoSetup from example 2
    s = create_geosetup()
    
    # get camera
    cam = s.points["cam"]
    print cam
    
    # calculate elevation profile between camera and source (Etna summit)
    # the camera altitude is 15 m above topography
    elev_profile = cam.get_elevation_profile(geo_point=s.points["source"])
    dist, dist_err, intersect, view_elevations, ax=\
        elev_profile.get_first_intersection(elev_angle=10.0,
                                            view_above_topo_m=15.0, 
                                            plot=True)
    
    ax.figure.savefig(join(save_path,
                           "ex4_out_1_elev_profile_intersect.png"))
    
    # Find the elvation angle corresponding to the horizon from the 
    # elevation profile
    elev, elev_secs, dist_secs =\
        elev_profile.find_horizon_elev(elev_start=10.0, elev_stop=20.0,
                                       step_deg=0.1,
                                       view_above_topo_m=15.0)
    
    fig, ax = subplots(1,1)
    ax.plot(elev_secs, dist_secs, "--x")
    ax.set_xlabel("Elevation angle [deg]")
    ax.set_ylabel("Retrieved distances [km]")
    ax.set_title("Elev horizon: %s" %elev)
    fig.savefig(join(save_path, "ex4_out_2_horizon_search.png"))
    
    return elev_profile, s
    
if __name__ == "__main__":
    profile, s = run_ex4()
    show()
    