"""
Introduction to GeoSetup class
===============================

"""

#%%
import geonum

#%%
# Define domain
# ---------------
#
# Create instance of Geosetup and specify domain by providing latitude and
# longitude of lower left (ll) and top right (tr) coordinates. Here,
# the summit region of the volcano Mt. Etna, Sicily, is used as an example.
domain = geonum.GeoSetup(id='Etna region (Sicily)',
                         lat_ll=37.5, lat_tr=38.0,
                         lon_ll=14.8, lon_tr=15.30)

#%%
# Add 2 GeoPoints to domain
# -------------------------
#
# Create GeoPoint at summit region of Mt. Etna and add to domain
etna = geonum.GeoPoint(latitude=37.751005,
                       longitude=14.993435,
                       altitude=3264.0,
                       name="Etna summit region")

domain.add_geo_point(etna)

#%%
#
# Create GeoPoint for location of volcanological observatory at Mt. Etna and
# add to domain
observatory = geonum.GeoPoint(latitude=37.765755,
                              longitude=15.016696,
                              auto_topo_access=True,
                              name="Observatory")

domain.add_geo_point(observatory)

#%%
#
# Load topographic data using SRTM database
# -----------------------------------------
topo = domain.load_topo_data(topo_access_mode='srtm')

#%%
#
# Make a beautiful map of the GeoSetup
# ------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import cartopy.mpl.ticker as cticker
import cartopy.crs as crs

# Initiate figure and cartopy GeoAxes for map plot
ax = geonum.plot_helpers.init_figure_with_geoaxes()

# Create meshgrid for domain for contour plot
x = topo.longitude
y = topo.latitude
X,Y = np.meshgrid(x, y)

# Separate land areas from sea areas for plotting below
seamask = np.ma.masked_less_equal(topo.data, 0)
landdata = np.ma.MaskedArray(topo.data, mask=seamask.mask)
seadata = np.ma.MaskedArray(topo.data, mask=~seamask.mask)

# Plot topographic land data
pdata = ax.contourf(X, Y, landdata, 50,
            transform=crs.PlateCarree(),
            cmap='Oranges')

# Add colorbar
cb = ax.figure.colorbar(pdata, label='Altitude [m]')

# Plot sea in a light blue
ax.contourf(X, Y, seadata, 50,
            transform=crs.PlateCarree(),
            colors='#e6f7ff')

# Plot contour lines
ax.contour(X, Y, topo.data, 10,
                linestyles='--',
                linewidths=0.1,
                colors='k',
                transform=crs.PlateCarree())

# Plot the 2 GeoPoints that were added to the GeoSetup
ax.scatter(
    x=[observatory.longitude, etna.longitude],
    y=[observatory.latitude, etna.latitude],
    color='w',
    marker='o',
    facecolor='none',
    s=16,
    transform=crs.PlateCarree()
)

# Annotate names of the 2 GeoPoints
_ = ax.annotate(
    etna.name,
    transform=crs.PlateCarree(),
    xy=(etna.longitude, etna.latitude), # location of observatory in plot
    xytext=(etna.longitude, etna.latitude-0.05),
    arrowprops=dict(color='white', lw=1, arrowstyle='->', shrinkB=4),
    ha='center', size=7, color='w', fontweight='bold'
    )

_ = ax.annotate(
    observatory.name,
    transform=crs.PlateCarree(),
    xy=(observatory.longitude, observatory.latitude), # location of observatory in plot
    xytext=(observatory.longitude, observatory.latitude+0.03),
    arrowprops=dict(color='white', lw=1, arrowstyle='->', shrinkB=4),
    ha='left', size=7, color='w', fontweight='bold'
    )


# Set title
ax.set_title(domain.id)

# Set x and y limits
ax.set_xlim([domain.ll.longitude, domain.tr.longitude])
ax.set_ylim([domain.ll.latitude, domain.tr.latitude])

ax = geonum.plot_helpers.set_map_ticks(ax)
plt.show()
