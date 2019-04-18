# -*- coding: utf-8 -*-
#
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

"""
Processing module of geonum library
"""
import numpy as np
from scipy.ndimage import map_coordinates

class LineOnGrid(object):
    """A class representing a line on a discrete grid"""
    def __init__(self, x0, y0, x1, y1, id=""):
        """Class initialisation
        
        :param int x0: start x coordinate
        :param int y0: start y coordinate
        :param int x1: stop x coordinate
        :param int y1: stop y coordinate
        :param str id: string for identification (optional)
        """
        self.id = id
        self.start  = [x0,y0]
        self.stop   = [x1,y1]
        
        self.profile_coords = None
        
        #self.check_coordinates()
        self._det_normal_vec()
        self.prepare_profile_coordinates()
    
    def check_coordinates(self):
        """Check if coordinates are in the right order and swap, if not"""
        if any([x < 0 for x in self.start]):
            raise ValueError("Invalid value encountered, coords must not be"
                " smaller than 0")
        if any([x < 0 for x in self.stop]):
            raise ValueError("Invalid value encountered, coords must not be"
                " smaller than 0")
        if self.start[0] > self.stop[0]:
            print("X coordinate of Start point is larger than X of stop point")
            print("Start and Stop will be exchanged")
            self.start, self.stop = self.stop, self.start
    
    """Processing functionality"""
    def prepare_profile_coordinates(self):
        """Prepare the evaluation coordinates as stack"""
        length = self.length
        x0, y0 = self.start     
        x1, y1 = self.stop
        x = np.linspace(x0, x1, length+1)
        y = np.linspace(y0, y1, length+1)
        self.profile_coords = np.vstack((y,x))
    
    @property
    def normal_vector(self):
        """Normal vector of line"""
        return self._det_normal_vec()
        
    @property
    def length(self):
        """Determine the length in grid coordinates"""
        del_x, del_y = self._delx_dely()
        return int(round(np.hypot(del_x, del_y)))
    
    def get_line_profile(self, array, **kwargs):
        """Retrieve line profile of data on a 2D array
        
        Uses method :func:`scipy.ndimage.map_coordinates` 
        
        Note
        ----
        The spline interpolation will fail if there are NaNs in the input data. 
        In this case, try using order=1 as additional input argument or use
        prefilter=False. For details see `here <https://docs.scipy.org/doc/
        scipy-0.14.0/reference/generated/scipy.ndimage.interpolation.
        map_coordinates.html>`__
        
        Parameters
        ----------
        array : ndarray
            2D data array
        **kwargs
            additional keyword args that are passed to
            :func:`scipy.ndimage.map_coordinates`
            
        Returns
        -------
        ndarray
            1D array containing the values along the line profile coordinates
        """
        # Extract the values along the line, using cubic interpolation
        zi = map_coordinates(array, self.profile_coords, **kwargs)
        return zi
    
    def plot_line_on_grid(self, array, ax=None, **kwargs):
        """Plot this line on an input array using imshow
        
        :param ndarray array: the data array
        :param **kwargs: additional keyword arguments for imshow
        
        """
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1,1)
        
        im = ax.imshow(array, cmap="gray", interpolation='none', 
                       origin="upper", **kwargs)
        fig.colorbar(im, ax=ax, shrink=0.9)
        ax.plot([self.start[0], self.stop[0]], [self.start[1], self.stop[1]],
                'co-')
            
        ax.set_xlim([0, array.shape[1] - 1])
        ax.set_ylim([array.shape[0] - 1, 0])
        return ax
    
    def plot_line_profile(self, array, ax=None, **kwargs):
        """Plot the line profile 
        
        :param ndarray array: the data array
        :param ax (None): axes object
        :param **kwargs: additional keyword arguments for imshow
        
        """
        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1,1)
        p = self.get_line_profile(array)
        ax.set_xlim([0, self.length])
        ax.grid()
        ax.plot(p)
        ax.set_title("Line profile")
        return ax
    
    def plot(self, array, fig_width=12):
        """High level plot function, calls:
        
            1. :func:`plot_line_on_grid`
            #. :func:`plot_line_profile`
            
        and puts them next to each other in a subplot
        
        :param ndarray array: the data array
        :param int fig_width: width of figure in inch
        
        """
        import matplotlib.pyplot as plt

        dx, dy = self._delx_dely()
        if dx > dy:
            r = float(dy) / dx
            h = int(fig_width * r * 2.5)
            fig, axes = plt.subplots(2,1, figsize=(fig_width, h))
        else:
            fig, axes = plt.subplots(1,2, figsize = (18, 6))
        self.plot_line_on_grid(array, ax=axes[0])
        self.plot_line_profile(array, ax=axes[1])
        plt.tight_layout()
    
    """'Private' functions"""
    def _delx_dely(self):
        """Returns length of x and y range"""
        return (abs(self.stop[0] - self.start[0]),
                abs(self.stop[1] - self.start[1]))
                        
    def _det_normal_vec(self):
        """Determines the normal vector of the line"""
        delx, dely = self._delx_dely()
        try:
            a = np.arctan(float(dely) / float(delx))
        except ZeroDivisionError:
            a = np.pi / 2
        return (np.sin(a), np.cos(a))
        
    """Magic methods"""
    def __str__(self):
        """String representation"""
        s = ("LineOnGrid %s\nStart (x,y): %s\nStop (x,y): %s\n" 
                                        %(self.id, self.start, self.stop))
        return s