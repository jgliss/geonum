# Geonum is a Python library for geographical calculations in 3D
# Copyright (C) 2017 Jonas Gliss (jonasgliss@gmail.com)
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License a
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import os

import numpy as np

from geonum import TOPO_INFO_FILE, LOCAL_TOPO_DIR

def order_of_magnitude(num):
    """
    Get exponent precision of input number

    E.g. 1 if input number is 14.2, 0 if input number is 3 or -1 if input
    number is 0.1.

    Note
    ----
    The order of magnitude here is calculated using floor(log10(num)), thus,
    num=9, would result in a value of 0, in contrary to other definitions (
    e.g. https://en.wikipedia.org/wiki/Order_of_magnitude).

    Parameters
    ----------
    num : float
        input number

    Returns
    -------
    int
        order of magnitude of input number

    """
    return int(np.floor(np.log10(np.abs(num))))


def all_topodata_search_dirs():
    """Returns a list of all directories that are searched for topography data

    Returns
    -------
    list
        list containing all search directories (absolute paths)
    """
    if LOCAL_TOPO_DIR is None:
        raise FileNotFoundError('geonum local topo directory not accessible')
    paths = []
    with open(TOPO_INFO_FILE, "r") as f:
        for line in f:
            p = line.strip()
            if os.path.exists(p) and not p in paths:
                paths.append(os.path.normpath(p))
    return paths


def check_and_add_topodir(local_dir):
    """Check if input directory is registered and if not, register it

    Parameters
    ----------
    local_dir : str
        directory that is supposed to be checked

    """
    fp = os.path.normpath(local_dir)
    if not os.path.exists(fp):
        raise ValueError(f'{local_dir} does not exist...')
    elif not fp in all_topodata_search_dirs():
        with open(TOPO_INFO_FILE, "a") as f:
            f.write(f'{fp}\n')
            print(f'Adding {fp} to file ~/.geonum/LOCAL_TOPO_DATA')


def isnum(val):
    """Checks if input is number (int or float) or and not nan

    Parameters
    ----------
    val
        input object to be checked

    Returns
    -------
    bool
        whether or not input object `val` is of numerical type
    """
    if isinstance(val, (int, float)) and not np.isnan(val):
        return True
    return False

def haversine_formula(lon0, lat0, lon1, lat1, radius=None):
    """Haversine formula to compute distances on a sphere

    Approximate horizontal distance between 2 points assuming a spherical
    earth.

    Parameters
    ----------
    lon0 : float
        longitude of first point in decimal degrees
    lat0 : float
        latitude of first point in decimal degrees
    lon1 : float
        longitude of second point in decimal degrees
    lat1 : float
        latitude of second point in decimal degrees
    radius : float
        average earth radius in km, defaults to 6371 km

    Returns
    -------
    float
        distance of both points in km
    """
    if radius is None:
        radius = 6371

    hav = lambda d_theta: np.sin(d_theta / 2.0) ** 2

    d_lon = np.radians(lon1 - lon0)
    d_lat = np.radians(lat1 - lat0)
    lat0 = np.radians(lat0)
    lat1 = np.radians(lat1)

    a = hav(d_lat) + np.cos(lat0) * np.cos(lat1) * hav(d_lon)
    c = 2 * np.arcsin(np.sqrt(a))

    return radius * c

