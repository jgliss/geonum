"""geonum example script 4

Elevation profiles and directional topographic analysis
"""

from matplotlib.pyplot import show, subplots
from ex2_geosetup_intro import create_geosetup
from os.path import join
from os import getcwd

### Set save directory for figures
save_path = join(getcwd(), "scripts_out")

def run_ex4():
    s = create_geosetup()
    
    cam = s.points["cam"]
    print cam
    
    # calculate elevation profile between camera and source
    elev_profile = cam.get_elevation_profile(geo_point = s.points["source"])
    dist, dist_err, intersect, view_elevations, ax =\
        elev_profile.get_first_intersection(elev_angle = 10.0,\
            view_above_topo_m = 15.0, plot = True)
    
    ax.figure.savefig(join(save_path, "ex4_out_1_elev_profile_intersect.png"))
    elev, elev_secs, dist_secs = elev_profile.find_horizon_elev(\
            elev_start = 10.0, elev_stop = 20.0, step_deg = 0.1,\
                                            view_above_topo_m = 15.0)
    
    fig, ax = subplots(1,1)
    ax.plot(elev_secs, dist_secs, "--x")
    ax.set_xlabel("Elevation angle [deg]")
    ax.set_ylabel("Retrieved distances [km]")
    ax.set_title("Elev horizon: %s" %elev)
    fig.savefig(join(save_path, "ex4_out_2_horizon_search.png"))
    
if __name__ == "__main__":
    run_ex4()
    show()