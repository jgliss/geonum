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
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import map_coordinates


class LineOnGrid(object):
    """A class representing a line on a discrete grid

    Attributes
    ----------
    start : list
        start location [x0, y0]
    stop : list
        stop location [x1, y1]

    Parameters
    ----------
    x0 : int
        start x coordinate
    y0 : int
        start y coordinate
    x1 : int
        end x coordinate
    y1 : int
        end y coordinate
    name : int
        string for identification (optional)

    Example
    -------

    >>> import numpy as np
    >>> from geonum import LineOnGrid
    >>> data = np.ones((10,15)) * np.arange(15)
    >>> line = LineOnGrid(2,2,5,5)
    >>> profile = line.get_line_profile(data)
    >>> print(profile)
    [2.         2.74392912 3.50304939 4.24920976 5.        ]

    """
    def __init__(self, x0, y0, x1, y1, name=None):
        if name is None:
            name = 'undefined'
        self.name = name
        self.start = [x0, y0]
        self.stop = [x1, y1]

        self._profile_coords = None

        self._det_normal_vec()

    @property
    def profile_coords(self) -> 'numpy.ndarray':
        """(y,x) profile coordinates"""
        if self._profile_coords is None:
            self._profile_coords = self._prepare_profile_coordinates()
        return self._profile_coords


    def _prepare_profile_coordinates(self):
        """Prepare the evaluation coordinates as stack"""
        length = self.length
        x0, y0 = self.start
        x1, y1 = self.stop
        x = np.linspace(x0, x1, length + 1)
        y = np.linspace(y0, y1, length + 1)
        return np.vstack((y, x))

    @property
    def normal_vector(self) -> tuple:
        """Normal vector of line (x, y)"""
        return self._det_normal_vec()

    @property
    def length(self) -> int:
        """Determine the length in units of grid indices

        Note
        ----
        The length is rounded to integer precision.
        """
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
        array : numpy.ndarray
            2D data array
        **kwargs
            additional keyword args that are passed to
            :func:`scipy.ndimage.map_coordinates`

        Returns
        -------
        numpy.ndarray
            1D array containing the values along the line profile coordinates
        """
        return map_coordinates(array, self.profile_coords, **kwargs)

    def plot_line_on_grid(self, array, ax=None, **kwargs):
        """Plot this line on an input array using imshow

        Parameters
        ----------
        array : numpy.ndarray
            the data array to which this line is mapped
        ax : matplotlib.axes.Axes, optional
            axes object
        **kwargs
            additional keyword arguments for imshow

        Returns
        -------
        matplotlib.axes.Axes
            axes instance used for plotting

        """
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        else:
            fig = ax.figure

        imshow_args = dict(cmap="gray", interpolation='none', origin="upper")
        imshow_args.update(kwargs)
        im = ax.imshow(array, **imshow_args)
        fig.colorbar(im, ax=ax, shrink=0.9)
        ax.plot([self.start[0], self.stop[0]],
                [self.start[1], self.stop[1]],
                'co-')

        ax.set_xlim([0, array.shape[1] - 1])
        ax.set_ylim([array.shape[0] - 1, 0])
        return ax

    def plot_line_profile(self, array, ax=None):
        """Plot profile of line in input data array

        Parameters
        ----------
        array : numpy.ndarray
            the data array to which this line is mapped
        ax : matplotlib.axes.Axes, optional
            axes object

        Returns
        -------
        matplotlib.axes.Axes
            axes instance used for plotting
        """
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        p = self.get_line_profile(array)
        ax.set_xlim([0, self.length])
        ax.grid()
        ax.plot(p)
        ax.set_title("Line profile")
        return ax

    def plot(self, array):
        """High level plot function, calls:

            1. :func:`plot_line_on_grid`
            2. :func:`plot_line_profile`

        and puts them next to each other in a subplot

        Parameters
        ----------
        array : numpy.ndarray
            2D data array to which this line is mapped

        Returns
        -------
        matplotlib.figure.Figure
            figure instance containing the plots.
        """
        fig, axes = plt.subplots(2, 1)

        self.plot_line_on_grid(array, ax=axes[0])
        self.plot_line_profile(array, ax=axes[1])
        plt.tight_layout()
        return fig

    def _delx_dely(self):
        """Returns length of x and y range"""
        return (abs(self.stop[0] - self.start[0]),
                abs(self.stop[1] - self.start[1]))

    def _det_normal_vec(self):
        """Determines the normal vector (orientation) of the line

        Returns
        -------
        tuple
            2-element tuple containing (x, y) elements of normal vector
        """
        delx, dely = self._delx_dely()
        try:
            a = np.arctan(float(dely) / float(delx))
        except ZeroDivisionError:
            a = np.pi / 2
        return (np.sin(a), np.cos(a))

    def __str__(self):
        """String representation"""
        s = (
            f"LineOnGrid {self.name}. Start (x,y): {self.start}. Stop (x,"
            f"y): {self.stop}")
        return s
