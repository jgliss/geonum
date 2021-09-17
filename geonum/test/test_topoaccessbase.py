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
import os
from geonum.conftest import does_not_raise_exception
from geonum.exceptions import TopoAccessError
from geonum import LOCAL_TOPO_DIR
from geonum import topoaccessbase as mod

def test_TopoAccessBase___init__():
    with pytest.raises(TypeError):
        mod.TopoAccessBase()

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



