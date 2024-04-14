"""
Berlin to Oslo to Paris
========================================

Introduction to GeoPoint and GeoVector3D.
"""

#%%
import matplotlib
matplotlib.rcParams['savefig.dpi'] = 300
import geonum

#%%
#
# Create some GeoPoints
# ---------------------
#
# GeoPoints represent locations on Earth specified via latitude, longtitude
# and altitude relative sea level. Here we create 2, points, one located in
# Paris, one in Oslo, and one in Berlin. The latitudes and longitudes for
# each city were looked up online, the altitudes can be retrieved automatically
# using the SRTM dataset, as shown below.

paris = geonum.GeoPoint(name='Paris',
                        latitude=48.8566,
                        longitude=2.3522,
                        auto_topo_access=True
                        )

#%%
print(paris)

#%%
oslo = geonum.GeoPoint(name='Oslo',
                       latitude=59.9139,
                       longitude=10.7522,
                       auto_topo_access=True
                       )

#%%
print(oslo)

#%%
berlin = geonum.GeoPoint(name='Berlin',
                       latitude=52.5200,
                       longitude=13.4050,
                       auto_topo_access=True
                       )

#%%
print(berlin)

#%%
#
# Compute connection vectors between the 3 locations
# --------------------------------------------------
#
# The connection vectors (dlat, dlon, dh) can be computed by subtracting 2
# GeoPoint instances.

paris_to_berlin = berlin - paris
print(paris_to_berlin)

berlin_to_oslo = oslo - berlin
print(berlin_to_oslo)

oslo_to_paris = paris - oslo
print(oslo_to_paris)

# Assign the corresponding starting points for the vectors
paris_to_berlin.set_anchor(paris)
berlin_to_oslo.set_anchor(berlin)
oslo_to_paris.set_anchor(oslo)

#%%
#
# Draw map of western Europe and add points and vectors
# -----------------------------------------------------
import cartopy.feature as cfeature
from cartopy import crs
import matplotlib.pyplot as plt

# Init map with lat lon range
ax = geonum.plot_helpers.init_figure_with_geoaxes(xlim=(0, 20),
                                                  ylim=(45, 65))

# Add some features to make it look nicer
ax.coastlines(alpha=0.3)
ax.add_feature(cfeature.LAND, color='#fff9e6', alpha=0.4)
ax.add_feature(cfeature.OCEAN, color='#e6f7ff', alpha=0.4)
ax.add_feature(cfeature.BORDERS, color='#808080', ls='--', alpha=0.4)

# Plot paris GeoPoint into map
ax = geonum.plot_helpers.plot_geopoint_into_map(
    ax, paris,
    annotate=True,
    marker='o', facecolor='none',
    color='#5c5c3d',
    annot_kwargs=dict(xytext=(paris.longitude, paris.latitude-2),
                      ha='center'))

# Plot oslo GeoPoint into map
ax = geonum.plot_helpers.plot_geopoint_into_map(
    ax, oslo,
    annotate=True,
    marker='o', facecolor='none',
    color='#5c5c3d',
    annot_kwargs=dict(xytext=(oslo.longitude, oslo.latitude+1.5),
                      ha='center')
)

# Plot berlin GeoPoint into map
ax = geonum.plot_helpers.plot_geopoint_into_map(
    ax, berlin,
    annotate=True,
    marker='o', facecolor='none',
    color='#5c5c3d',
    annot_kwargs=dict(xytext=(berlin.longitude, berlin.latitude-1.5),
                      ha='center')
)

# Draw GeoVector3D Oslo->Paris into map
ax = geonum.plot_helpers.plot_geovector3d_into_map(ax, oslo_to_paris,
                                                   color='#5c5c3d',
                                                   ls='--', lw=1)

# annotate text specifying distance between points:
half_way_vec = oslo_to_paris/2
half_way = oslo + half_way_vec

# distance between Oslo and Paris
distance = oslo_to_paris.dist_hor
ax.annotate(text=f'Oslo/Paris:\n {distance:.0f} km',
            xy=(half_way.longitude, half_way.latitude),
            xytext=(half_way.longitude-1.5, half_way.latitude+1),
            color='#5c5c3d', fontweight='bold', ha='left', fontsize=6)

# Draw GeoVector3D Paris->Berlin into map
ax = geonum.plot_helpers.plot_geovector3d_into_map(ax, paris_to_berlin,
                                                   color='#5c5c3d',
                                                   ls='--', lw=1)

half_way_vec = paris_to_berlin/2
half_way = paris + half_way_vec
distance = paris_to_berlin.dist_hor
ax.annotate(text=f'Paris/Berlin:\n {distance:.0f} km',
            xy=(half_way.longitude, half_way.latitude),
            xytext=(half_way.longitude+2, half_way.latitude-1),
            color='#5c5c3d', fontweight='bold', ha='center', fontsize=6)

# Draw GeoVector3D Berlin->Oslo into map
ax = geonum.plot_helpers.plot_geovector3d_into_map(ax, berlin_to_oslo,
                                                   color='#5c5c3d',
                                                   ls='--', lw=1)

half_way_vec = berlin_to_oslo/2
half_way = berlin + half_way_vec
distance = berlin_to_oslo.dist_hor
ax.annotate(text=f'Berlin/Oslo:\n {distance:.0f} km',
            xy=(half_way.longitude, half_way.latitude),
            xytext=(half_way.longitude+1, half_way.latitude+1.8),
            color='#5c5c3d', fontweight='bold', ha='center', fontsize=6)

# Add x and y ticks
geonum.plot_helpers.set_map_ticks(ax)

# Add title
ax.set_title('Berlin to Oslo to Paris')

plt.show()
