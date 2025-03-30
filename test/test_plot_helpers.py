import cartopy
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import numpy as np
import pytest

from geonum import plot_helpers as mod
from test.conftest import does_not_raise_exception

# Use a non-rendering backend for matplotlib
plt.switch_backend('Agg')

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


@pytest.mark.parametrize('deg,ha', [
    (45, "center"), (90, "left"), (30, "right")
])
def test_rotate_xtick_labels(deg, ha):
    _, ax = plt.subplots(1,1) 
    ax.set_xticks(range(5))
    ax.set_xticklabels([f"Label {i}" for i in range(5)])

    # Call the function to rotate xtick labels
    updated_ax = mod.rotate_xtick_labels(ax, deg=deg, ha=ha)

    # Verify the rotation and horizontal alignment
    for label in updated_ax.get_xticklabels():
        assert label.get_rotation() == deg
        assert label.get_ha() == ha

@pytest.mark.parametrize('vmin, vmax, cmap, test_values, expected_colors', [
    (-10, 10, None, [-10, 0, 10], [(1.0, 0.0, 0.0, 1.0), (0.5, 0.5, 0.5, 1.0), (0.0, 0.0, 1.0, 1.0)]),  # Default cmap
    (-5, 15, plt.cm.coolwarm, [-5, 5, 15], [plt.cm.coolwarm(0.0), plt.cm.coolwarm(0.5), plt.cm.coolwarm(1.0)]),  # Custom cmap
])
def test_shifted_color_map_colors(vmin, vmax, cmap, test_values, expected_colors):
    shifted_cmap = mod.shifted_color_map(vmin, vmax, cmap)
    assert isinstance(shifted_cmap, LinearSegmentedColormap)