import cartopy
import numpy as np
import pytest

from geonum import plot_helpers as mod
from geonum.conftest import does_not_raise_exception


@pytest.mark.parametrize('projection', [
    None, cartopy.crs.TransverseMercator(), cartopy.crs.PlateCarree(),
    cartopy.crs.RotatedPole()
])
def test_init_figure_with_geoaxes(projection):
    ax = mod.init_figure_with_geoaxes()
    assert isinstance(ax, cartopy.mpl.geoaxes.GeoAxes)


@pytest.mark.parametrize('ax,xticks,yticks,xres,yres,raises', [
    (mod.init_figure_with_geoaxes(), None, None,
     [-180, -120, -60, 0., 60., 120., 180.],
     [-90., -60., -30., 0., 30., 60., 90.],
     does_not_raise_exception()),
    (mod.init_figure_with_geoaxes(xlim=(15.05, 15.08), ylim=(-0.023, 0.012)),
     None, None,
     [15.05, 15.056, 15.062, 15.068, 15.074, 15.08],
     [-0.023, -0.016, -0.009, -0.002, 0.005, 0.012],
     does_not_raise_exception()),
])
def test_set_map_ticks(ax, xticks, yticks, xres, yres, raises):
    with raises:
        ax = mod.set_map_ticks(ax, xticks, yticks)
        assert ax is ax
        xt = ax.get_xticks()
        yt = ax.get_yticks()
        np.testing.assert_allclose(xt, xres)
        np.testing.assert_allclose(yt, yres)
