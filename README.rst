|build-status|

`geonum <https://github.com/jgliss/geonum>`__ provides functionality for geographical calculations in three dimensions and includes interfaces for accessing and processing of topographic data. Most of the features (e.g. distance retrievals) are based on the two fundamental objects *GeoPoint* and *GeoVector3D* which are inherited from the respective 2D base classes of the `LatLon23 module <https://pypi.org/project/LatLon23>`_ and were expanded including the 3rd dimension (altitude).
Geonum features online access to topographic data from the SRTM dataset, using the
`SRTM module <https://pypi.python.org/pypi/SRTM.py/>`_. Furthermore, topographic data from the `ETOPO1 Dataset <https://www.ngdc.noaa.gov/mgg/global/global.html>`_ is supported.

News / Notifications
====================

- Now also available and tested in Python > 3
- Changed requirement **LatLon** -> **LatLon23** to support Python 3

Planned changes
===============

- Refactoring of ``basemap`` dependency to ``cartopy``
- ``TopoData`` should be based on ``xarray.DataArray``
- Support for more topographic datasets, interpolation of gaps in topodata

Copyright
=========

Copyright (C) 2017 Jonas Gliss (jonasgliss@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License a published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see `here <http://www.gnu.org/licenses/>`__.

Requirements
============

It is recommended to use `Anaconda <https://www.continuum.io/downloads>`_ as Python package manager as it includes many of the required dependencies and makes life easier when it comes to the installation or upgrade of non-straight forward installations of additional requirements such as OpenCV or basemap.

  - numpy
  - matplotlib >= 1.4.3
  - matplotlib `basemap <https://pypi.python.org/pypi/basemap/1.0.7>`_ (*installation of this module may not be straight forward - especially on Windows machines -, please follow the instructions provided on the web page*)
  - LatLon23 >= 1.0.7

    - LatLon23 requires installation of `pyproj <https://pypi.python.org/pypi/pyproj/>`_

  - Scipy (including `scipy.ndimage <https://docs.scipy.org/doc/scipy-0.18.1/reference/ndimage.html>`_)

**Optional dependencies (to use extra features)**

  - OpenCV (used for changing resolution of topographic elevation maps, for installation remarks, see `here <http://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_setup/py_setup_in_windows/py_setup_in_windows.html>`_)
  - netCDF4 (needed in case `ETOPO1 <https://www.ngdc.noaa.gov/mgg/global/global.html>`_ data acess is required).

Installation
============

Geonum can be installed from `PyPi <https://pypi.python.org/pypi/geonum>`_ using::

  pip install geonum

or from source by downloading and extracting the latest release. After navigating to the source folder (where the setup.py file is located) call::

  python setup.py install

If the installation fails make sure, that all dependencies (see above) are installed correctly. geonum is tested for Python 2.7.

Instructions and code documentation
===================================

The code documentation of geonum is hosted on `Read The Docs <http://geonum.readthedocs.io/en/latest/index.html>`_

Get started
===========

After installation try running the `example scripts <http://geonum.readthedocs.io/en/latest/examples.html>`_ in order to test the installation. The scripts are also meant to provide an easy start into the main features of geonum.

Supported ETOPO1 files
======================

In order to use topography data from the ETOPO1 dataset, please download and unzip one of the following files to the package folder *geonum/local_topo_data/*.
Tested and supported are the following two files (grid registered):

  1. Ice surface: ETOPO1_Ice_g_gmt4.grd (download `here <https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz>`__)
  2. Bedrock: ETOPO1_Bed_g_gmt4.grd (download `here <https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gmt4.grd.gz>`__)

The local installation folder can be accessed via::

  import geonum
  print geonum.LOCAL_TOPO_PATH

If a valid data file is stored in this folder, it will be detected automatically. It is also possible to store the topodata at another location (e.g. <data_path>). In this case, the local path to the folder containing the topograph files needs to be provided, e.g.::

  import geonum
  access = geonum.topodata.TopoDataAccess(mode = "etopo1", local_path = <data_path>)

If the path is valid, it will be added to the installation file *LOCAL_TOPO_PATHS.txt*

.. |build-status| image:: https://travis-ci.com/jgliss/geonum.svg?branch=master
    :target: https://travis-ci.com/jgliss/geonum
    
.. |docs| image:: https://readthedocs.org/projects/geonum/badge/?version=latest
:target: https://geonum.readthedocs.io/en/latest/?badge=latest
:alt: Documentation Status
