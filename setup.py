import os
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
    name        =   'geonum',
    version     =   os.environ['GITHUB_REF_NAME'],
    author      =   'Jonas Gliss',
    author_email=   'jonasgliss@gmail.com',
    license     =   'GPLv3',
    url         =   'https://github.com/jgliss/geonum',
    package_dir =   {'geonum'     :   'geonum'},
    packages    =   find_packages(exclude=['contrib', 'docs', 'tests*']),
    include_package_data    =   True,
    package_data=   {'geonum'     :   ['local_topo_data/*.rst',
                                       'local_topo_data/*.txt']},
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.,
        'Programming Language :: Python :: 3'
    ],
    description = 'Toolbox for 3D geonumerical calculations and atmospheric composition',
    long_description = readme,
    long_description_content_type='text/x-rst'
)