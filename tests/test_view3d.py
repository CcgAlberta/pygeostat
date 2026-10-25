#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the 3D meshes in pygeostat (needs PyVista)"""

import pytest
import numpy as np
import pandas as pd
import pygeostat as gs

pv = pytest.importorskip('pyvista')


@pytest.fixture
def intervals():
    """Two drill hole intervals with start and end coordinates."""
    return pd.DataFrame({
        'startx': [0.0, 0.0], 'starty': [0.0, 0.0], 'startz': [100.0, 90.0],
        'endx': [0.0, 0.0], 'endy': [0.0, 0.0], 'endz': [90.0, 80.0],
        'Au': [1.5, np.nan],
    })


def test_drillhole_mesh_points_from_datafile():
    """Without start/end, each row is one point at the DataFile's own x, y, z columns."""
    data = gs.ExampleData('oilsands')
    mesh = gs.drillhole_mesh(data, var='Bitumen')
    xyz = data.data[[data.x, data.y, data.z]].to_numpy(dtype=float)
    assert mesh.n_points == len(data)
    assert mesh.n_cells == len(data)
    np.testing.assert_array_equal(mesh.points, xyz)
    np.testing.assert_array_equal(mesh.cell_data['Bitumen'], data.data['Bitumen'].to_numpy(dtype=float))


def test_drillhole_mesh_cylinders(intervals):
    """Each interval becomes a cylinder from start to end, carrying its own value."""
    mesh = gs.drillhole_mesh(intervals, var='Au', start=('startx', 'starty', 'startz'),
                             end=('endx', 'endy', 'endz'), radius=2.0, n_sides=20)
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    assert (zmin, zmax) == (80.0, 100.0)
    assert xmax - xmin == pytest.approx(4.0)
    assert ymax - ymin == pytest.approx(4.0, abs=0.1)
    # Cells above z = 90 belong to the first interval, below it to the second
    centres = mesh.cell_centers().points[:, 2]
    values = mesh.cell_data['Au']
    assert np.all(values[centres > 90] == 1.5)
    assert np.all(np.isnan(values[centres < 90]))


def test_drillhole_mesh_errors(intervals):
    """Missing coordinates, or only one of start/end, raise a clear error."""
    with pytest.raises(ValueError, match='could not be'):
        gs.drillhole_mesh(intervals, var='Au')
    with pytest.raises(ValueError, match='both `start` and `end`'):
        gs.drillhole_mesh(intervals, var='Au', start=('startx', 'starty', 'startz'))


def test_drillhole_mesh_zero_length_warns(intervals, capsys):
    """An interval with the same start and end is not drawn, and the user is told."""
    intervals.loc[1, 'endz'] = intervals.loc[1, 'startz']
    mesh = gs.drillhole_mesh(intervals, var='Au', start=('startx', 'starty', 'startz'),
                             end=('endx', 'endy', 'endz'))
    assert '1 interval(s) have the same start and end' in capsys.readouterr().out
    assert mesh.bounds[4:] == (90.0, 100.0)


def test_drillhole_mesh_text_variable(intervals):
    """Text categories get a clear error instead of a failed float conversion."""
    intervals['Rock'] = ['granite', 'diorite']
    with pytest.raises(ValueError, match="'Rock' is not numeric"):
        gs.drillhole_mesh(intervals, var='Rock', start=('startx', 'starty', 'startz'),
                          end=('endx', 'endy', 'endz'))


@pytest.fixture
def samples():
    """Two drill holes, given out of order; hole B has a 50 m stretch with no samples."""
    return pd.DataFrame({
        'dh': ['B', 'A', 'B', 'A', 'B', 'B'],
        'x': [10.0, 0.0, 10.0, 0.0, 10.0, 10.0],
        'y': [0.0] * 6,
        'z': [95.0, 90.0, 100.0, 100.0, 40.0, 35.0],
        'Au': [2.0, 1.0, 3.0, 4.0, 5.0, 6.0],
    })


def test_drillhole_mesh_trace(samples):
    """Samples are joined down each hole, each cylinder takes its upper sample's value, and the
    gap in hole B (50 m, ten times the 5 m median spacing) is not bridged."""
    mesh = gs.drillhole_mesh(samples, var='Au', x='x', y='y', z='z', dh='dh', trace=True,
                             radius=1.0)
    values = mesh.cell_data['Au']
    centres = mesh.cell_centers().points
    hole_a = centres[:, 0] < 5
    assert set(values[hole_a]) == {4.0}         # A: 100 -> 90, value of the sample at 100
    assert set(values[~hole_a]) == {3.0, 5.0}   # B: 100 -> 95 and 40 -> 35
    assert not np.any((centres[:, 2] > 45) & (centres[:, 2] < 90))


def test_drillhole_mesh_trace_max_gap(samples):
    """A larger max_gap bridges the gap, so hole B gets its 95 -> 40 cylinder too."""
    mesh = gs.drillhole_mesh(samples, var='Au', x='x', y='y', z='z', dh='dh', trace=True,
                             max_gap=60.0)
    assert set(mesh.cell_data['Au']) == {2.0, 3.0, 4.0, 5.0}


def test_drillhole_mesh_trace_from_datafile():
    """A DataFile's own drill hole id column is used without naming it."""
    data = gs.ExampleData('oilsands')
    mesh = gs.drillhole_mesh(data, var='Bitumen', trace=True)
    assert mesh.n_cells > 0
    assert mesh.bounds[4] >= data.data[data.z].min() - 2.0


def test_drillhole_mesh_trace_missing_coordinate(capsys):
    """A sample with no elevation is skipped; the rest of the hole is still joined."""
    data = pd.DataFrame({'dh': ['A'] * 4, 'x': [0.0] * 4, 'y': [0.0] * 4,
                         'z': [100.0, 95.0, 90.0, np.nan], 'Au': [1.0, 2.0, 3.0, 4.0]})
    mesh = gs.drillhole_mesh(data, var='Au', x='x', y='y', z='z', dh='dh', trace=True)
    assert '1 sample(s) with missing coordinates' in capsys.readouterr().out
    assert set(mesh.cell_data['Au']) == {1.0, 2.0}
    assert mesh.bounds[4:] == (90.0, 100.0)
