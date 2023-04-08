from geonum import GeoPoint, delete_local_srtm_files
from geonum.utils import get_local_srtm_files


def test_get_local_srtm_files():
    # make sure there is one file at least
    p = GeoPoint(45, 15, auto_topo_access=True)
    assert p.altitude > 0
    files = get_local_srtm_files()
    assert isinstance(files, list)
    assert len(files) >= 1


def test_delete_local_srtm_files():
    delete_local_srtm_files()
    assert len(get_local_srtm_files()) == 0