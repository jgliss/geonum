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
        tick_format = _get_tick_formatter(xticks[0])

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