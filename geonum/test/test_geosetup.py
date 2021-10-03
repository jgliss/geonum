import numpy as np
import pytest
import numpy.testing as npt
from numpy import ndarray

from geonum.conftest import does_not_raise_exception
from geonum import GeoPoint, GeoVector3D, TopoDataAccess, TopoData
from geonum import geosetup as mod
from geonum.conftest import skip_srtm
from geonum.exceptions import OutOfDomain

@pytest.mark.parametrize(
    'points,vectors,lat_ll,lon_ll,lat_tr,lon_tr,id,topo_access_mode,'
    'local_topo_path,cmap_vecs,raises',[
        (None,None,None,None,None,None,None,None,None,None,
         does_not_raise_exception()),
        (GeoPoint(),None,None,None,None,None,None,None,None,None,
         does_not_raise_exception()),
        ('Blaaaa',None,None,None,None,None,None,None,None,None,
         pytest.raises(ValueError)),
        (None,GeoVector3D(0,0),None,None,None,None,None,None,None,None,
         does_not_raise_exception()),
        (None, 'Blaaaa',None,None,None,None,None,None,None,None,
         pytest.raises(ValueError)),
        (None, None, 42, 14,43, 15, None, None, None, None,
         does_not_raise_exception()),
        (None, None, 42, None, None, None, None, None, None, None,
         does_not_raise_exception()),

    ])
def test_GeoSetup__init__(points,vectors,lat_ll,lon_ll,lat_tr,lon_tr,
                          id,topo_access_mode, local_topo_path,
                          cmap_vecs,raises):
    with raises:
        stp=mod.GeoSetup(points, vectors, lat_ll, lon_ll, lat_tr, lon_tr,
                         id, topo_access_mode, local_topo_path,
                         cmap_vecs)
        if id is None:
            id ='MyGeoSetup'
        if topo_access_mode is None:
            topo_access_mode = "srtm"
        if cmap_vecs is None:
            cmap_vecs = 'Greens'
        assert stp.id == id
        assert stp.topo_access_mode == topo_access_mode
        assert stp.local_topo_path == local_topo_path
        assert stp._cmap_vecs == cmap_vecs

@pytest.fixture(scope='module')
def empty_setup():
    return mod.GeoSetup()

@skip_srtm
def test_GeoSetup_create_test_setup():
    stp = mod.GeoSetup.create_test_setup()
    assert isinstance(stp, mod.GeoSetup)
    assert len(stp.points) == 4
    assert len(stp.vectors) == 2

@pytest.fixture(scope='module')
def dummy_setup():
    stp = mod.GeoSetup.create_test_setup()
    return stp

def test_GeoSetup_ll():
    stp = mod.GeoSetup()
    with pytest.raises(AttributeError):
        stp.ll
    stp.ll = GeoPoint(45, 15, name='ll')
    stp.set_borders_from_points()
    stp.load_topo_data()
    assert isinstance(stp.ll, GeoPoint)
    with pytest.raises(ValueError):
        stp.ll = 42
    stp.ll = GeoPoint(44,15,name='ll')

def test_GeoSetup_tr():
    stp = mod.GeoSetup()
    with pytest.raises(AttributeError):
        stp.tr
    stp.tr = GeoPoint(45, 15, name='tr')
    stp.set_borders_from_points()
    stp.load_topo_data()
    assert isinstance(stp.tr, GeoPoint)
    with pytest.raises(ValueError):
        stp.ll = 42
    stp.tr = GeoPoint(44,15,name='tr')

def test_GeoSetup_lon_ll(dummy_setup):
    npt.assert_allclose(dummy_setup.lon_ll, 14.985, atol=0.001)

def test_GeoSetup_lat_ll(dummy_setup):
    npt.assert_allclose(dummy_setup.lat_ll, 37.743, atol=0.001)

def test_GeoSetup_lon_tr(dummy_setup):
    npt.assert_allclose(dummy_setup.lon_tr, 15.025, atol=0.001)

def test_GeoSetup_lat_tr(dummy_setup):
    npt.assert_allclose(dummy_setup.lat_tr, 37.774, atol=0.001)

def test_GeoSetup_delta_lon(dummy_setup):
    npt.assert_allclose(dummy_setup.delta_lon, 0.03, atol=0.01)

def test_GeoSetup_delta_lat(dummy_setup):
    npt.assert_allclose(dummy_setup.delta_lat, 0.03, atol=0.01)

def test_GeoSetup_borders_set():
    assert mod.GeoSetup().borders_set == False
    assert mod.GeoSetup(points=GeoPoint(0,0)).borders_set == True

def test_GeoSetup_center_coordinates():
    gs = mod.GeoSetup()
    with pytest.raises(AttributeError):
        coords = gs.center_coordinates
    gs.add_geo_point(GeoPoint(-0.5,-0.5,name='ll'))
    with pytest.raises(AttributeError):
        coords = gs.center_coordinates
    gs.add_geo_point(GeoPoint(0.5, 0.5, name='tr'))
    coords = gs.center_coordinates
    assert coords == (0,0)

def test_GeoSetup_magnitude(dummy_setup):
    npt.assert_allclose(dummy_setup.magnitude, 4.898, atol=0.01)

def test_GeoSetup_topo_access(empty_setup):
    assert isinstance(empty_setup.topo_access, TopoDataAccess)

def test_GeoSetup_cmap_vecs(empty_setup):
    assert empty_setup._cmap_vecs == 'Greens'
    from matplotlib.colors import LinearSegmentedColormap
    assert isinstance(empty_setup.cmap_vecs, LinearSegmentedColormap)

def test_GeoSetup_has_points(empty_setup):
    assert not empty_setup.has_points()
    empty_setup.add_geo_point(GeoPoint(name='null'), assert_in_domain=False)
    assert empty_setup.has_points()
    empty_setup.delete_geo_point('null')
    assert not empty_setup.has_points()

@skip_srtm
def test_GeoSetup_has_topo_data(dummy_setup):
    assert dummy_setup.has_topo_data() == False
    dummy_setup.get_topo()
    assert dummy_setup.has_topo_data() == True
    dummy_setup.topo_data = None

def test_GeoSetup_set_local_topo_path(tmpdir):
    gs = mod.GeoSetup()
    fp = str(tmpdir)
    gs.set_local_topo_path(tmpdir)
    assert gs.local_topo_path == fp
    with pytest.raises(FileExistsError):
        gs.set_local_topo_path('/blaaaaaa/blub')

@pytest.mark.parametrize('new_mode,local_path,raises', [
    ('bla', None, does_not_raise_exception()),
    ('bla', 'blub', pytest.raises(FileExistsError)),
    ('bla', '/', does_not_raise_exception()),

])
def test_GeoSetup_change_topo_mode(new_mode,local_path,raises):
    gs = mod.GeoSetup()
    with raises:
        gs.change_topo_mode(new_mode,local_path)

@skip_srtm
def test_GeoSetup_get_topo():
    stp = mod.GeoSetup.create_test_setup()
    stp.get_topo()
    assert isinstance(stp.topo_data, TopoData)

@skip_srtm
def test_GeoSetup_load_topo_data():
    stp = mod.GeoSetup()
    stp.add_geo_point(GeoPoint(45,15))
    stp.load_topo_data()
    assert isinstance(stp.topo_data, TopoData)

def test_GeoSetup_has_point():
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    assert stp.has_point('ll') == True
    assert stp.has_point('bla') == False
    stp.add_geo_point(GeoPoint(10, 20, name='bla'), assert_in_domain=False)
    assert stp.has_point('bla') == True

def test_GeoSetup_has_vector():
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    assert stp.has_vector('bla') == False
    stp.add_geo_vector(GeoVector3D(1,1,1,name='bla'))
    assert stp.has_vector('bla') == True

def test_GeoSetup_contains_coordinate():
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    assert stp.contains_coordinate(46, 16) == True
    assert stp.contains_coordinate(30, 16) == False

@pytest.mark.parametrize('pt,assert_in_domain,overwrite_existing,raises', [
    (42,False,True,pytest.raises(ValueError)),
    (GeoPoint(0,0,name='bla'),False,True,does_not_raise_exception()),
    (GeoPoint(0,0,name='bla'),False,False,pytest.raises(ValueError)),
    (GeoPoint(0,0,name='bla'),True,True,pytest.raises(OutOfDomain)),
    (GeoPoint(0,0,name='ll'),True,False,does_not_raise_exception()),
    (GeoPoint(0,0,name='tr'),True,False,does_not_raise_exception()),
])
def test_GeoSetup_add_geo_point(pt,assert_in_domain,overwrite_existing,raises):
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    stp.points['bla'] = GeoPoint(47,16,name='bla')
    stp.topo_data = TopoData(lats=[45,50], lons=[15,20], data=np.ones((2,2)))
    with raises:
        stp.add_geo_point(pt=pt,assert_in_domain=assert_in_domain,
                          overwrite_existing=overwrite_existing)
        assert pt.name in stp.points


@pytest.mark.parametrize('pts,assert_in_domain,raises', [
    ([GeoPoint(0,0)], True, pytest.raises(OutOfDomain)),
    ([GeoPoint(48,16)], True, does_not_raise_exception()),
    ([GeoPoint(48,16), GeoPoint(49,16)], True, pytest.raises(ValueError)),
    ([GeoPoint(48,16,name='p1'), GeoPoint(49,16,name='p2')], True,
    does_not_raise_exception()),
    ([GeoPoint(48,16,name='p1'), GeoPoint(49,21)], True, pytest.raises(
        OutOfDomain)),
    ([GeoPoint(48,16,name='p1'), GeoPoint(49,21)], False, does_not_raise_exception()),
])
def test_GeoSetup_add_geo_points(pts,assert_in_domain,raises):
    stp = mod.GeoSetup(lat_ll=45,lon_ll=15, lat_tr=50, lon_tr=20)
    with raises:
        stp.add_geo_points(*pts, assert_in_domain=assert_in_domain)
        assert len(stp.points) == len(pts) + 2


@pytest.mark.parametrize('vec,overwrite_existing,raises', [
    (42, False, pytest.raises(ValueError)),
    (GeoVector3D(1,1,1,name='bla'), False, pytest.raises(ValueError)),
    (GeoVector3D(1,1,1,name='bla'), True, does_not_raise_exception()),
    (GeoVector3D(1,1,1,name='blub'), False, does_not_raise_exception()),
])
def test_GeoSetup_add_geo_vector(vec,overwrite_existing,raises):
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    stp.vectors['bla'] = GeoVector3D(1,1,1,name='bla')
    with raises:
        stp.add_geo_vector(vec=vec,
                          overwrite_existing=overwrite_existing)
        assert vec.name in stp.vectors


@pytest.mark.parametrize('vecs,raises', [
    ([GeoVector3D(1,1,1,name='bla')], does_not_raise_exception()),
    ([GeoVector3D(1,1,1,name='bla'),
      GeoVector3D(1,1,1,name='bla')], does_not_raise_exception()),
])
def test_GeoSetup_add_geo_vector(vecs,raises):
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    with raises:
        stp.add_geo_vectors(*vecs)

def test_GeoSetup_delete_geo_point():
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    with pytest.raises(ValueError):
        stp.delete_geo_point('bla')
    stp.add_geo_point(GeoPoint(0,0,name='bla'),assert_in_domain=False)
    assert 'bla' in stp.points
    stp.delete_geo_point('bla')
    assert not 'bla' in stp.points

def test_GeoSetup_delete_geo_vector():
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    with pytest.raises(ValueError):
        stp.delete_geo_vector('bla')
    stp.add_geo_vector(GeoVector3D(1,1,1,name='bla'))
    assert 'bla' in stp.vectors
    stp.delete_geo_vector('bla')
    assert not 'bla' in stp.vectors

def test_GeoSetup_new_geo_point():
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    with pytest.raises(OutOfDomain):
        stp.new_geo_point(latitude=10, longitude=20)

    stp.new_geo_point(latitude=46, longitude=19, name='name')
    stp.new_geo_point(47, 19, name='name2')

def test_GeoSetup__all_lats_lons():
    stp = mod.GeoSetup(lat_ll=45, lon_ll=15, lat_tr=50, lon_tr=20)
    result = stp._all_lats_lons()
    assert isinstance(result[0], ndarray)
    assert isinstance(result[1], ndarray)
    assert list(result[0]) == [45, 50]
    assert list(result[1]) == [15, 20]

@pytest.mark.parametrize('pts,extend_km,to_square,raises,delta_lat,'
                         'delta_lon', [
    ([],1,True,pytest.raises(AttributeError),None,None),
    ([GeoPoint(0,0)],1,True,does_not_raise_exception(),0.013,0.013),
    ([GeoPoint(0,0)],111.15/np.sqrt(2),True,does_not_raise_exception(),
     1.005,0.999),
    ([GeoPoint(0, 0, name='p1'), GeoPoint(1, 10, name='p2')], 0, True,
     does_not_raise_exception(),
     10.061, 10),
    ([GeoPoint(0, 0, name='p1'), GeoPoint(1, 10, name='p2')], 0, False,
     does_not_raise_exception(),
     1, 10),
    ([GeoPoint(0, 0, name='p1'), GeoPoint(1, 10, name='p2')], 100,
     False,
     does_not_raise_exception(),
     2.279, 11.271),
])
def test_set_borders_from_points(pts,extend_km,to_square,raises,delta_lat,
                                 delta_lon):
    stp = mod.GeoSetup(points=pts, init_borders=False)
    with raises:
        stp.set_borders_from_points(extend_km,to_square)
        npt.assert_allclose(delta_lat, stp.delta_lat, atol=0.001)
        npt.assert_allclose(delta_lon, stp.delta_lon, atol=0.001)



near_Etna = GeoPoint(latitude=37.751, longitude=14.993,
                     name="near Etna", auto_topo_access=True)

@skip_srtm
@pytest.mark.parametrize('p,radius,result,raises', [
    (None,None,[],pytest.raises(ValueError)),
    (near_Etna,None,['Etna'],does_not_raise_exception()),
    (near_Etna,0.01,[],does_not_raise_exception()),
    (near_Etna,20,['Etna', 'Observatory', 'll', 'tr'],does_not_raise_exception()),
])
def test_GeoSetup_points_close(p,radius,result,raises):
    stp = mod.GeoSetup.create_test_setup()
    with raises:
        val = stp.points_close(p,radius)
        assert val == result


