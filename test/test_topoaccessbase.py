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
Access and handling of topographic data
"""
import pytest
import numpy as np
import numpy.testing as npt
from geonum.conftest import does_not_raise_exception
from geonum.exceptions import TopoAccessError
from geonum import LOCAL_TOPO_DIR
from geonum import topoaccessbase as mod

@pytest.fixture(scope='module')
def etopo_acc():
    return mod.Etopo1Access(None, None, False, False)

@pytest.mark.parametrize('local_path, file_name, check_access,search_database,'
                         'raises', [
    (None, None, False, True,does_not_raise_exception()),
    (None, None, False, False,does_not_raise_exception()),
    (None, None, True, False,pytest.raises(TopoAccessError)),

])
def test_Etopo1Access___init__(tmpdir, local_path, file_name, check_access,
                               search_database,raises):

    with raises:
        val = mod.Etopo1Access(local_path=local_path,
                               file_name=file_name,
                               check_access=check_access,
                               search_database=search_database)
        if file_name is None:
            assert val._file_name == "ETOPO1_Ice_g_gmt4.grd"
        else:
            assert val._file_name == file_name

        if local_path is None:
            assert val._local_path == LOCAL_TOPO_DIR
        else:
            assert val._local_path == local_path

def test_Etopo1Access_check_access(etopo_acc):
    assert not etopo_acc.check_access()

def test_Etopo1Access_check_access_MOCK_GET_DATA(
        monkeypatch):
    def mock_get_data(self, lat0, lon0, lat1=None, lon1=None):
        return 42

    monkeypatch.setattr(mod.Etopo1Access, 'get_data', mock_get_data)

    acc = mod.Etopo1Access(local_path=None,
              file_name=None,
              check_access=False,
              search_database=False)
    assert acc.check_access()

@pytest.mark.parametrize('lat0,lon0,lat1,lon1,res', [
    (6,1,None,None,
     (np.asarray([6]), np.asarray([1]), [1,1], [1,1])),
    (8,4,-1,-1,
     (np.asarray([5,6,7]), np.asarray([0,1,2,3]), [0, 2], [0, 3]))
])
def test__init_lons_lats(etopo_acc,lat0,lon0,lat1,lon1,res):
    lats_all = np.asarray([5,6,7])
    lons_all=np.asarray([0,1,2,3])
    result = etopo_acc._init_lons_lats(
        lats_all, lons_all,lat0,lon0,lat1,lon1)
    for i, item in enumerate(res):
        npt.assert_array_equal(result[i], item)

def test_init_Etopo1_invalid_file(tmpdir):
    file_name = 'invalid.grd'
    path = tmpdir / file_name
    with open(path, 'w') as f:
        f.write('bla')

    assert path.exists()

    acc = mod.Etopo1Access(local_path=tmpdir,
                           file_name=file_name,
                           check_access=True,
                           search_database=False)





