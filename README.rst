|build-status| |docs| |cov| |code-quality|

`geonum <https://github.com/jgliss/geonum>`__ provides functionality for geographical calculations in three dimensions and includes interfaces for accessing and processing of topographic data. Most of the features (e.g. distance retrievals) are based on the two fundamental objects *GeoPoint* and *GeoVector3D* which are inherited from the respective 2D base classes of the `LatLon23 module <https://pypi.org/project/LatLon23>`_ and were expanded to support also the vertical dimension.
Geonum features online access to topographic data from the SRTM dataset, using the
`SRTM module <https://pypi.python.org/pypi/SRTM.py/>`_. Furthermore, topographic data from the `ETOPO1 Dataset <https://www.ngdc.noaa.gov/mgg/global/global.html>`_ is supported.

Copyright
=========

Copyright (C) 2017 Jonas Gliß (jonasgliss@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License a published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see `here <http://www.gnu.org/licenses/>`__.

Requirements
============

It is recommended to use the *conda* Python package manager (get it e.g. via `miniconda <https://docs.conda.io/en/latest/miniconda.html>`__)

Please see `requirements.txt <https://github.com/jgliss/geonum/blob/master/requirements.txt>`__ for a list of requirements. Note that these are only the minimum required dependencies. Further, optional requirements are listed below. Note, that dependent on your OS and python version, some of these may be tricky to install:

**Optional dependencies (to use extra features)**

- basemap (plotting of maps)
- OpenCV
- netCDF4

Installation
============

geonum is tested and can be installed both for Python >= 3.6 and on all common OS (Windows, linux and OSX).

Installation via conda
----------------------

The easiest way to install geonum is to install the `latest release via the conda-forge channel <https://anaconda.org/conda-forge/geonum>`_::

  conda install -c conda-forge geonum

This will install all requirements as well.

Installation via pip or from source
-------------------------------------

Please make sure to install all requirements beforehand (see above). You may do this by downloading the `requirements file <https://github.com/jgliss/geonum/blob/master/requirements.txt>`__ and calling from the command line::

  pip install -r requirements.txt

or - if you use Python 3 - by creating a new conda environment using the provided `conda environment file <https://github.com/jgliss/geonum/blob/master/geonum_env.yml>`_::

  conda env create -n geonum_env -f geonum_env.yml

After installing the requirements, geonum can be installed from `PyPi <https://pypi.python.org/pypi/geonum>`_ using::

  pip install geonum

or from source by downloading and extracting the latest release or cloning the repository. After navigating to the source folder (where the setup.py file is located) call::

  python setup.py install

Instructions and code documentation
===================================

The code documentation of geonum is hosted on `Read The Docs <http://geonum.readthedocs.io/>`_

Getting started
===============

After installation try running the `example scripts <http://geonum.readthedocs.io/en/latest/examples.html>`_ in order to test the installation. The scripts are also meant to provide an easy start into the main features of geonum.

Supported ETOPO1 files
======================

In order to use topography data from the ETOPO1 dataset, please download and unzip one of the following files to the package folder *geonum/local_topo_data/*.
Tested and supported are the following two files (grid registered):

  1. Ice surface: ETOPO1_Ice_g_gmt4.grd (download `here <https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz>`__)
  2. Bedrock: ETOPO1_Bed_g_gmt4.grd (download `here <https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gmt4.grd.gz>`__)

The local installation folder can be accessed via::

  import geonum
  print(geonum.LOCAL_TOPO_PATH)

If a valid data file is stored in this folder, it will be detected automatically. It is also possible to store the topodata at another location (e.g. <data_path>). In this case, the local path to the folder containing the topograph files needs to be provided, e.g.::

  import geonum
  access = geonum.topodata.TopoDataAccess(mode="etopo1", local_path=<data_path>)

If the path is valid, it will be added to the installation file *LOCAL_TOPO_PATHS.txt*

Planned changes (for v2.0.0)
============================

See below for my (@jgliss) personal wish-list of new features, help is more than welcome as I have to work on geonum mostly in my spare time.

- Refactoring of ``basemap`` dependency to ``cartopy``
- Base ``TopoData`` on ``xarray.DataArray``
- Support for more topographic datasets, interpolation of gaps in topodata, merging and interpolation of different topographic datasets

.. |build-status| image:: https://travis-ci.com/jgliss/geonum.svg?branch=master
    :target: https://travis-ci.com/jgliss/geonum

.. |docs| image:: https://readthedocs.org/projects/geonum/badge/?version=latest
    :target: https://geonum.readthedocs.io/en/latest/?badge=latest

.. |cov| image:: https://codecov.io/gh/jgliss/geonum/branch/main-dev/graph/badge.svg?token=802DAZA1W9
    :target: https://codecov.io/gh/jgliss/geonum

.. |code-quality| image:: https://www.codefactor.io/repository/github/jgliss/geonum/badge
   :target: https://www.codefactor.io/repository/github/jgliss/geonum
   :alt: CodeFactor
