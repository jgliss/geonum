# -*- coding: utf-8 -*-
import pytest
import numpy as np
from contextlib import contextmanager
from geonum import NETCDF_AVAILABLE


# custom skipif marker that is used below for test functions that
# require SRTM.py to be installed
def srtm_works():
    """SRTM access is broken, see
    https://github.com/tkrajina/srtm.py/issues/51
    """
    import srtm
    try:
        srtm.get_data().get_file(50,8)
        return True
    except Exception as e:
        print('failed to access SRTM data...')
        if not str(e).startswith('Cannot retrieve https'):
            raise
        return False


skip_srtm = pytest.mark.skipif(not srtm_works(),
                               reason=('Skipping SRTM database tests since '
                                       'SRTM data cannot be accessed.'))

skip_netcdf4 = pytest.mark.skipif(NETCDF_AVAILABLE==False,
                                  reason='Skipping tests requiring NetCDF4 library.')


@contextmanager
def does_not_raise_exception():
    yield

@pytest.fixture
def topodata_2x2():
    pass


if __name__ == '__main__':
    print(srtm_works())