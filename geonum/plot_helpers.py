"""
Helpers for plotting of maps etc
--------------------------------

.. note::

    BETA stage, code might undergo revisions

"""
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

from geonum.helpers import order_of_magnitude


def init_figure_with_geoaxes(projection=None,
                             xlim=(-180,180),
                             ylim=(-90, 90)):
    """
    Initiate a matplotlib figure containing 1 cartopy GeoAxes instance

    Can be used to plot maps.

    Parameters
    ----------
    projection : cartopy.crs.CRS, optional
        Projection used, if None, then :class:`cartopy.crs.PlateCarree` is
        used. Defaults to None.
    xlim : tuple, optional
        longitude range of map, defaults to whole globe (-180,180).
    ylim : tuple, optional
        latitude range of map, defaults to whole globe (-90,90).

    Returns
    -------
    cartopy.mpl.geoaxes.GeoAxes

    """
    if projection is None:
        projection = crs.PlateCarree()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return ax

def _get_tick_formatter(value):
    """
    Get tick formatting string for map plots

    Parameters
    ----------
    value : float
        value used to estimate formatter (e.g. one of the longitude
        coordinates).

    Returns
    -------
    str
        formatting string

    """
    digits = 2 - order_of_magnitude(value)
    digits = 0 if digits < 0 else digits
    return f'.{digits:d}f'

def _calculate_map_ticks(ax):
    """
    Calculate longitude (x) and latitude (y) ticks for map plots.

    Parameters
    ----------
    ax : cartopy.mpl.geoaxes.GeoAxes
        map axes instances.

    Returns
    -------
    tuple
        2-element tuple containing x and y ticks
    """
    lonleft, lonright = ax.get_xlim()

    num_lonticks = 7 if lonleft == -lonright else 6
    xticks = np.linspace(lonleft, lonright, num_lonticks)

    latleft, latright = ax.get_ylim()
    num_latticks = 7 if latleft == - latright else 6
    yticks = np.linspace(latleft, latright, num_latticks)
    return (xticks, yticks)

def set_map_ticks(ax, xticks=None, yticks=None, tick_format=None):
    """Set or update ticks in instance of GeoAxes object (cartopy)

    Note
    ----
    The input GeoAxes must be setup with :class:`cartopy.crs.PlateCarree`
    projection.

    Parameters
    ----------
    ax : cartopy.GeoAxes
        map axes instance (needs to be set up with PlateCarree projection).
    xticks : iterable, optional
        ticks of x-axis (longitudes). If None, it is retrieved automatically
        based on xlim of input axes. Defaults to None.
    yticks : iterable, optional
        ticks of y-axis (latitudes). If None (or if xticks is None),
        it is retrieved automatically based on ylim of input axes. Defaults
        to None.
    tick_formatter : str, optional
        string formatter for axes labels in the plot.

    Returns
    -------
    cartopy.GeoAxes
        modified axes instance
    """
    if not isinstance(ax.projection, crs.PlateCarree):
        raise AttributeError('Map ticks can only be added to GeoAxes with '
                             'PlateCarree projection.')

    if xticks is None:
        xticks, yticks = _calculate_map_ticks(ax)

    if tick_format is None:
        tick_format = _get_tick_formatter(np.mean(xticks))

    ax.set_xticks(xticks, crs=crs.PlateCarree())
    ax.set_yticks(yticks, crs=crs.PlateCarree())

    lon_formatter = LongitudeFormatter(number_format=tick_format,
                                       degree_symbol='°',
                                       dateline_direction_label=True)
    lat_formatter = LatitudeFormatter(number_format=tick_format,
                                      degree_symbol='°')

    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    return ax

def _infer_text_location(ax,p):
    info = dict(ha='center')
    left,right = ax.get_xlim()
    dx = right - left
    bottom,top = ax.get_ylim()
    dy = top - bottom
    lat, lon = p.latitude, p.longitude
    distleft = abs(left - lon) / dx
    if distleft < 0.3:
        info['ha'] = 'left'
    elif distleft > 0.7:
        info['ha'] = 'right'

    yoffs = dy*0.08 # put text 5% of range above or below the points latitude
    distbottom = abs(bottom - lat) / dx
    if distbottom > 0.6: # point is in upper 40% of map, put text below
        yoffs *= -1
    info['xytext'] = (lon, lat+yoffs)
    return info

def plot_topo_contourf(ax, topo, oceans_separate=True, levels=50, cmap=None):
    if cmap is None:
        cmap = 'Oranges'

    # Create meshgrid for domain for contour plot
    X,Y = topo.init_mesh()

    # Separate land areas from sea areas for plotting below
    if oceans_separate:
        seamask = np.ma.masked_less_equal(topo.data, 0)
        data = np.ma.MaskedArray(topo.data, mask=seamask.mask)
        if np.any(seamask):
            seadata = np.ma.MaskedArray(topo.data, mask=~seamask.mask)
            # Plot sea in a light blue
            ax.contourf(X, Y, seadata, 50,
                        colors='#e6f7ff')
    else:
        data = topo.data

    # Plot topographic land data
    pdata = ax.contourf(X, Y, data, levels=levels,
                cmap=cmap)

    return ax, pdata

def plot_geopoint_into_map(ax, p, annotate=True, annot_kwargs=None, **kwargs):

    if annot_kwargs is None:
        annot_kwargs = {}

    ax.scatter(x=[p.longitude],y=[p.latitude], **kwargs)

    if annotate:
        annot_loc = _infer_text_location(ax,p)
        try:
            c = kwargs['color']
        except KeyError:
            c = 'k'
        annot = dict(
            text = p.name,
            xy = (p.longitude,
                  p.latitude),  # location of observatory in plot
            arrowprops = dict(color=c, lw=1, arrowstyle='->', shrinkB=4),
            size=7, color=c, fontweight='bold', **annot_loc
        )
        annot.update(annot_kwargs)
        ax.annotate(**annot)
    return ax

def plot_geovector3d_into_map(ax, vec, **kwargs):
    if vec.anchor is None:
        raise AttributeError(
            f'input vector {vec} does not have an anchor point assigned. '
            f'Please set anchor location using vec.add_anchor')

    start = vec.anchor
    end = start + vec
    ax.plot([start.longitude, end.longitude],
            [start.latitude, end.latitude], **kwargs)
    return ax