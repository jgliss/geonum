# -*- coding: utf-8 -*-
"""
Helper methods for geonum library
"""

from numpy import radians, cos, arcsin, sin, sqrt, array, linspace, hstack
import matplotlib.cm as colormaps
import matplotlib.colors as colors

def haversine_formula(lon0, lat0, lon1, lat1, radius=6371.0):
    """Haversine formula
    
    Approximate horizontal distance between 2 points assuming a spherical 
    earth
    
    :param float lon0: longitude of first point in decimal degrees
    :param float lat0: latitude of first point in decimal degrees
    :param float lon1: longitude of second point in decimal degrees
    :param float lat1: latitude of second point in decimal degrees
    :param float radius (6371.0): average earth radius in km
    """
    hav = lambda d_theta: sin(d_theta / 2.0) ** 2
    
    d_lon = radians(lon1 - lon0)
    d_lat = radians(lat1 - lat0)
    lat0 = radians(lat0)
    lat1 = radians(lat1)
 
    a = hav(d_lat) + cos(lat0) * cos(lat1) * hav(d_lon)
    c = 2 * arcsin(sqrt(a))
 
    return radius * c

def approximate_connection_vector(lon0, lat0, lon1, lat1, len_lat_km=111.20):
    """Returns approximate connection vector between two points
    
    :param float lon0: longitude of first point in decimal degrees
    :param float lat0: latitude of first point in decimal degrees
    :param float lon1: longitude of second point in decimal degrees
    :param float lat1: latitude of second point in decimal degrees
    
    Careful: Only approximative, suited for small distance ( < 100km)
    """    
    
        
    lat = (lat0 + lat1) / 2 * 0.01745
    dx = len_lat_km * cos(lat) * (lon1 - lon0)
    dy = len_lat_km * (lat1 - lat0)
    return array((dx,dy))
    
def shifted_color_map(vmin, vmax, cmap=None):
    """Shift center of a diverging colormap to value 0
    
    .. note::
    
        This method was found `here <http://stackoverflow.com/questions/
        7404116/defining-the-midpoint-of-a-colormap-in-matplotlib>`_ 
        (last access: 17/01/2017). Thanks to `Paul H <http://stackoverflow.com/
        users/1552748/paul-h>`_ who provided it.
    
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and if you want the
    middle of the colormap's dynamic range to be at zero level
    
    :param vmin: lower end of data value range
    :param vmax: upper end of data value range
    :param cmap: colormap (if None, use default cmap: seismic)
    
    :return: 
        - shifted colormap
        
    """
    #midpoint = 1 - abs(im.max())/(abs(im.max()) + abs(im.min()))
    if cmap is None:
        cmap = colormaps.seismic
        
    midpoint = 1 - abs(vmax)/(abs(vmax) + abs(vmin))
    
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = linspace(0, 1, 257)

    # shifted index to match the data
    shift_index = hstack([
        linspace(0.0, midpoint, 128, endpoint=False), 
        linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    #newcmap = colors.LinearSegmentedColormap('shiftedcmap', cdict)
    #register_cmap(cmap=newcmap)

    return colors.LinearSegmentedColormap('shiftedcmap', cdict)