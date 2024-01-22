Installation
############

You can install geonum via conda, pip or from source.

.. _install_conda:

Via conda
----------

The easiest way to install geonum and all requirements is via the conda
package manager. You can install geonum into a pre-existing conda
environment via::

	conda install -c conda-forge geonum

Or into a new conda environment (recommended) named *geonum* via::

	conda create -c conda-forge --name geonum geonum

This will install the latest release of geonum including all requirements.
Alternatively, you may install via pip or from source as described in the following.


Via pip
-------

This will install the latest released version of geonum.

.. note::

    This is the same pacakge as distributed via *conda-forge* (see prev.
    Sect. :ref:`install_conda`)

::

	pip install geonum


Install from source (into a conda environment)
----------------------------------------------

Please make sure to install all requirements (see Section :ref:`dependencies`)
before installing geonum from source.

To install geonum from source, please download and extract the
`latest release <https://github.com/jgliss/geonum/releases>`__ or clone the
`repo <https://github.com/jgliss/geonum/>`__) and install from the project
root directory (that contains a file *setup.py*) using::

	pip install --no-deps .

The `--no-deps` option will ensure that only the pyearocom package is
installed, preserving the conda environment.

Alternatively, if you plan to apply local changes to the geonum source code,
you may install in editable mode (i.e. setuptools "develop mode")::

	pip install --no-deps -e .[test]

You may also download and extract (or clone) the
`GitHub repo <https://github.com/jgliss/geonum>`__ to install the very
latest (not yet released) version of geonum. Note, if you install in editable
mode, make sure you do not have geonum installed already in the site packages
directory, check e.g. `conda list geonum`.


.. _dependencies:

Dependencies
------------

The list of required libraries is provided via the file `geonum_env.yml`,
located in the project root directory:

.. literalinclude:: ../geonum_env.yml

Install dependencies using conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to install all requirements into a new conda environment, you can
call (from project root)::

	conda env create -n geonum -f geonum_env.yml

This will create a new conda environment called *geonum* which can be
activated using::

	conda activate geonum

Alternatively, you can update an existing environment. First, activate the
existing environment, and then install the dependencies using::

	conda env update -f=geonum_env.yml
