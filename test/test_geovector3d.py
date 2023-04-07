import numpy as np
import pytest
import geonum.geovector3d as mod

def test_init_NO_INPUT():
    with pytest.raises(ValueError):
        mod.GeoVector3D()

def test_init_DOUBLE_INPUT():
    v = mod.GeoVector3D(dx=1, dy=1, azimuth=223, dist_hor=1.2456)
    assert v.dx == 1
    assert v.dy == 1
@pytest.mark.parametrize('kwargs,mag,azim,elev', [

    (dict(dx=2, dy=0, dz=0),2,90,0),
    (dict(dx=2, dy=0, elevation=0),2,90,0),
    (dict(dx=2, dy=0, elevation=90),2,90,0),
    (dict(dx=2, dy=0, elevation=45),2**(1.5),90,45),
    (dict(dx=1, dy=-1, dz=0),2**(0.5),135,0),
    (dict(azimuth=140,dist_hor=111,dz=0),111,140,0),
    (dict(azimuth=140,dist_hor=1,dz=1000),2**0.5,140,45),
])
def test_init(kwargs,mag,azim,elev):
    v = mod.GeoVector3D(**kwargs)
    _mag = v.magnitude
    np.testing.assert_allclose(v.magnitude, mag, rtol=1e-5)
    np.testing.assert_allclose(v.azimuth, azim, rtol=1e-5)
    np.testing.assert_allclose(v.elevation, elev, rtol=1e-5)
def test_truediv_CASE1():
    v1 = mod.GeoVector3D(dx=1,dy=1,dz=1)
    v2 = mod.GeoVector3D(dx=1, dy=1, dz=1)

    res = v1.__truediv__(v2)
    assert isinstance(res, mod.GeoVector3D)
    assert res.magnitude == v1.magnitude

def test_truediv_CASE2():
    v1 = mod.GeoVector3D(dx=1,dy=1,dz=1)

    res = v1.__truediv__(1)
    assert isinstance(res, mod.GeoVector3D)
    assert res.magnitude == v1.magnitude