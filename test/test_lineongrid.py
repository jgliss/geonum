import matplotlib.figure
import numpy as np
import pytest

from geonum import LineOnGrid


@pytest.mark.parametrize('line,result', [
    (LineOnGrid(0, 0, 0, 1), (1, 0)),
    (LineOnGrid(0, 0, 1, 0), (0, 1)),
    (LineOnGrid(0, 0, 1, 1), (1 / np.sqrt(2), 1 / np.sqrt(2))),
])
def test_LineOnGrid__det_normal_vec(line, result):
    vec = line._det_normal_vec()
    assert isinstance(vec, tuple)
    np.testing.assert_allclose(vec, result, atol=1e-4)


@pytest.mark.parametrize('line,result', [
    (LineOnGrid(0, 0, 0, 1), (0, 1)),
    (LineOnGrid(0, 0, 1, 0), (1, 0)),
    (LineOnGrid(2, 12, 44, 14), (42, 2)),
])
def test_LineOnGrid__delx_dely(line, result):
    vec = line._delx_dely()
    assert isinstance(vec, tuple)
    np.testing.assert_allclose(vec, result, atol=1e-4)


def test_LineOnGrid___str__():
    line = LineOnGrid(1, 1, 2, 2, name='bla')
    st = str(line)
    assert st == 'LineOnGrid bla. Start (x,y): [1, 1]. Stop (x,y): [2, 2]'


def test_LineOnGrid_plot():
    val = LineOnGrid(1, 1, 4, 6).plot(np.ones((20, 30)))
    assert isinstance(val, matplotlib.figure.Figure)