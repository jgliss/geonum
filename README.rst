geonum provides functionality for 3D geo numerical calculations
including online access and handling of topography data. The module is based on two fundamental classes "GeoPoint" and "GeoVector3D" which are inherited from the respective base classes of the `LatLon module <https://pypi.python.org/pypi/LatLon>`_ and were expanded for the 3rd (elevation) dimension.
Topography access is supported for SRTM data, which is handled by the 
`SRTM module <https://pypi.python.org/pypi/SRTM.py/>`_, and the `Etopo1 Dataset <https://www.ngdc.noaa.gov/mgg/global/global.html>`_ for which the data needs to be downloaded and stored locally.

Requirements
------------

It is recommended to use the package manager `Anaconda <https://www.continuum.io/downloads>`_ since it includes many of the required dependencies and makes life easier when it comes to installation or upgrade of non-straight forward installations of additional requirements (e.g. opencv, basemap)

  - numpy
  - matplotlib
  - matplotlib `basemap <https://pypi.python.org/pypi/basemap/1.0.7>`_ (*installation of this module may not be straight forward - especially on Windows machines -, please follow the instructions provided on the web page*) 
  - `LatLon <https://pypi.python.org/pypi/LatLon>`_
  
    - LatLon requires installation of `pyproj <https://pypi.python.org/pypi/pyproj/>`_
    
  - Scipy (including `scipy.ndimage <https://docs.scipy.org/doc/scipy-0.18.1/reference/ndimage.html>`_)

**Optional dependencies (to use extra features)**

  - OpenCV (used for changing resolution of topographic elevation maps, for installation remarks, see `here <http://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_setup/py_setup_in_windows/py_setup_in_windows.html>`_)
  - netCDF4 (needed in case `Etopo1 <https://www.ngdc.noaa.gov/mgg/global/global.html>`_ data acess is required).


Installation
------------
geonum can be installed from source by downloading and extracting the latest release. After navigating to the source folder (where the setup.py file is located), call::

  python setup.py install
  
If the installation fails make sure, that all dependencies (see above) are installed correctly. Geonum was only tested for Python 2.7.

Instructions and code documentation
-----------------------------------

The code documentation of geonum is hosted on `Read The Docs <http://geonum.readthedocs.io/en/latest/index.html>`_

Get started
-----------

After installation try running the `example scripts <http://geonum.readthedocs.io/en/latest/examples.html>`_ in order to test the installation. The scripts are also meant to provide an easy start into the main features of geonum.