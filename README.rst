geonum provides functionality for 3D geo numerical calculations
including online access and handling of topography data (:mod:`srtm.py` module). The base classes :class:`GeoPoint` and :class:`GeoVector3D` are inherited from the :class:`LatLon` and :class:`GeoVector` from the `LatLon module <https://pypi.python.org/pypi/LatLon>`_ and were expanded for the 3rd (elevation) dimension.
Topography access is supported for SRTM data, which is handled by the 
`SRTM module <https://pypi.python.org/pypi/SRTM.py/>`_, and the `Etopo1 Dataset <https://www.ngdc.noaa.gov/mgg/global/global.html>`_ which needs to be downloaded and stored locally.

Requirements
------------
 
  - numpy
  - matplotlib
  - matplotlib `basemap <https://pypi.python.org/pypi/basemap/1.0.7>`_ (installation of this module may not be straight forward - especially on Windows machines -, please follow the instructions provided on the web page) 
  - `LatLon <https://pypi.python.org/pypi/LatLon>`_
  - Scipy (including `scipy.ndimage <https://docs.scipy.org/doc/scipy-0.18.1/reference/ndimage.html>`_)

**Optional dependencies (to use extra features)**

  - OpenCV (used for changing resolution of topographic elevation maps, for installation remarks, see `here <http://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_setup/py_setup_in_windows/py_setup_in_windows.html>`_
  - netCDF4 (needed in case `Etopo1 <https://www.ngdc.noaa.gov/mgg/global/global.html>`_ data acess is required).


Generally strongly recommend using `Anaconda <https://www.continuum.io/downloads>`_ for package management since it includes most of the required dependencies and is updated on a regular basis. 

Installation
------------
geonum can be installed from source by downloading and extracting the latest release. After navigating to the source folder (where the setup.py file is located), call::

  python setup.py install
  
Geonum is registered on Alternatively, you may also use::

  pip install geonum
  
If the installation fails, make sure, that all dependencies (see above) are installed correctly. Geonum was only tested for Python 2.7.

Get started
-----------

After installation try running the example scripts, which also provides an easy start into the functionality of geonum.

.. literalinclude:: ../scripts/ex1_basic_objects.py

.. literalinclude:: ../scripts/ex2_geosetup_intro.py

.. literalinclude:: ../scripts/ex3_mapping_basics.py