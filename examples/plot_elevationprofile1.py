"""
Retrieval of elevation profiles
===============================

"""

import geonum

# Summit region of Mt. Etna
etna = geonum.GeoPoint(latitude=37.751005,
                       longitude=14.993435,
                       altitude=3264.0,
                       name="Etna summit")

# Position of volcanological observatory of Mt. Etna
observatory = geonum.GeoPoint(latitude=37.765755,
                              longitude=15.016696,
                              auto_topo_access=True,
                              name="Observatory")

# Create a lat/lon domain conatining both points
domain = geonum.GeoSetup(points=[etna, observatory])

# Get lower left and upper right corner coordinates of domain
pll = domain.ll
ptr = domain.tr

# Instantiate topo data access class for SRTM access
topo_access = geonum.TopoDataAccess(mode='srtm')
# Retrieve topographic data within domain ranges
topo_data = topo_access.get_data(lat0=pll.latitude, lon0=pll.longitude,
                                 lat1=ptr.latitude, lon1=ptr.longitude)

# Instantiate elevation profile with Etna Observatory as start point and
# Etna summit as end point and provide topographic data
eprof = geonum.ElevationProfile(observer=observatory,
                                endpoint=etna,
                                topo_data=topo_data)

# Calculate elevation profile using
eprof.det_profile(interpolate=False)

# Plot the profile
ax = eprof.plot()
ax.set_title('Etna summit region')

# Annotate names of the 2 points to the plot
ax.annotate(observatory.name,
            xy=(0, observatory.altitude), # location of observatory in plot
            xytext=(0, observatory.altitude+50),
            arrowprops=dict(color='black', lw=1, arrowstyle='->', ),
            ha='left', size=11
            )

# horizontal distance of summit from observatory
summit_dist = (etna-observatory).dist_hor

ax.annotate(etna.name,
            xy=(summit_dist, etna.altitude), # location of observatory in plot
            xytext=(summit_dist-0.5, etna.altitude+50),
            arrowprops=dict(color='black', lw=1, arrowstyle='->', ),
            ha='left', size=11
            )

if __name__=='__main__':
    import matplotlib.pyplot as plt
    plt.show()