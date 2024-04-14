import pytest
import os, shutil
import numpy as np
import numpy.testing as npt
from geonum import LOCAL_TOPO_DIR, TOPO_INFO_FILE
from geonum import helpers as h

@pytest.mark.parametrize('num,result', [
    (10, 1), (1000, 3), (1.2,0), (0.023, -2)
])
def test_order_of_magnitude(num,result):
    assert h.order_of_magnitude(num) == result

def test_all_topodata_search_dirs():
    all_dirs = h.all_topodata_search_dirs()
    assert isinstance(all_dirs, list)
    assert LOCAL_TOPO_DIR in all_dirs
    for dirloc in all_dirs:
        assert os.path.exists(dirloc)

def test_check_and_add_topodir():
    tmp_file = os.path.join(LOCAL_TOPO_DIR, 'TMP.file')
    topofile_tmp = shutil.copyfile(TOPO_INFO_FILE, tmp_file)
    all_dirs = h.all_topodata_search_dirs()
    num_init = len(all_dirs)
    h.check_and_add_topodir(LOCAL_TOPO_DIR)
    assert len(h.all_topodata_search_dirs()) == num_init # it is already in there and should thus not be added

    homedir = os.path.expanduser('~')
    if not homedir in all_dirs:
        h.check_and_add_topodir(homedir)
        assert len(h.all_topodata_search_dirs()) == num_init + 1

    invaliddir = 'bla/blub.xuspf/nla'
    if os.path.exists(invaliddir):
        print('What a surprise indeed... path {invaliddir} exists :-O')
    else:
        with pytest.raises(ValueError):
            h.check_and_add_topodir(invaliddir)

    # go back to initial state
    topofile_tmp = shutil.copyfile(tmp_file, TOPO_INFO_FILE)
    os.remove(tmp_file)

@pytest.mark.parametrize('val,expectation', [
    (0, True), (0.1, True), ('1', False),
    ({}, False), ([1], False), ((1,2,), False)
    ])
def test_isnum(val,expectation):
    assert h.isnum(val) == expectation


@pytest.mark.parametrize('lon0,lat0,lon1,lat1,radius,value', [
    (-77.037852,38.898556,-77.043934,38.897147,None,0.549156),
    (-90,0,90,0,None,6371 * np.pi), # half around the equator
    (-90,90,90,-90,None,6371 * np.pi), # pole to pole
    (-90,0,90,0,500,500 * np.pi), # half around the equator
    ])
def test_haversine_formula(lon0,lat0,lon1,lat1,radius,value):
    val = h.haversine_formula(lon0, lat0, lon1, lat1, radius)
    npt.assert_allclose(val, value, rtol=1e-5)


if __name__ == "__main__":
    import sys
    pytest.main(sys.argv)


