#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the 3D meshes in pygeostat (needs PyVista)"""

import pytest
import numpy as np
import pandas as pd
import pygeostat as gs

pv = pytest.importorskip('pyvista')


@pytest.fixture(autouse=True)
def close_scenes():
    """Close the figures and PyVista plotters each test leaves open."""
    from matplotlib import pyplot as plt
    yield
    plt.close('all')
    pv.close_all()


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


@pytest.fixture
def grid3d():
    """pygeostat's 40 x 40 x 40 example grid and its grid definition."""
    return gs.ExampleData('3d_grid'), gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1])


def test_block_mesh_uniform(grid3d):
    """Each cell holds the value of the block that GridDef itself finds at the cell's centre."""
    data, griddef = grid3d
    mesh = gs.block_mesh(data, var='value', griddef=griddef)
    assert mesh.n_cells == griddef.count()
    assert mesh.bounds == (0.0, 80.0, 0.0, 80.0, 0.0, 40.0)
    centres = mesh.cell_centers().points
    idx, _ = griddef.get_index(centres[:, 0], centres[:, 1], centres[:, 2])
    np.testing.assert_array_equal(mesh.cell_data['value'], data.data['value'].to_numpy()[idx])


def test_block_mesh_griddef_from_datafile(grid3d):
    """A DataFile's own griddef is used without passing it."""
    data, griddef = grid3d
    data.griddef = griddef
    assert gs.block_mesh(data, var='value').n_cells == griddef.count()


def test_block_mesh_errors(grid3d):
    """No griddef, or a griddef that doesn't match the data, raise a clear error."""
    data, griddef = grid3d
    with pytest.raises(ValueError, match='griddef'):
        gs.block_mesh(data.data, var='value')
    with pytest.raises(ValueError, match='64000 rows but the grid has 32000 blocks'):
        gs.block_mesh(data, var='value', griddef=gs.GridDef([40, 1, 2, 40, 1, 2, 20, 0.5, 1]))


@pytest.fixture
def subblocks():
    """Three sub-blocks of different sizes, in the standard format: centre and sizes."""
    return pd.DataFrame({
        'x': [5.0, 12.5, 15.0], 'y': [5.0, 2.5, 5.0], 'z': [5.0, 2.5, 7.5],
        'xsize': [10.0, 5.0, 10.0], 'ysize': [10.0, 5.0, 10.0], 'zsize': [10.0, 5.0, 5.0],
        'Au': [1.0, 2.0, 3.0],
    })


def test_block_mesh_subblocks(subblocks):
    """Each row becomes one box with its own centre, sizes and value."""
    mesh = gs.block_mesh(subblocks, var='Au', x='x', y='y', z='z',
                         size=('xsize', 'ysize', 'zsize'))
    assert mesh.n_cells == 3
    assert mesh.bounds == (0.0, 20.0, 0.0, 10.0, 0.0, 10.0)
    np.testing.assert_allclose(mesh.cell_centers().points, subblocks[['x', 'y', 'z']])
    volumes = mesh.compute_cell_sizes(length=False, area=False)['Volume']
    np.testing.assert_allclose(volumes, subblocks[['xsize', 'ysize', 'zsize']].prod(axis=1))
    np.testing.assert_array_equal(mesh.cell_data['Au'], subblocks['Au'])


def test_block_mesh_subblocks_errors(subblocks, grid3d):
    """Both griddef and size, or sub-blocks without centres, raise a clear error."""
    _, griddef = grid3d
    with pytest.raises(ValueError, match='not both'):
        gs.block_mesh(subblocks, var='Au', griddef=griddef, size=('xsize', 'ysize', 'zsize'))
    with pytest.raises(ValueError, match='could not be'):
        gs.block_mesh(subblocks, var='Au', size=('xsize', 'ysize', 'zsize'))


def layers(viewer):
    """The PyVista actors of the meshes in a scene, in the order they were added."""
    return [a for a in viewer.plotter.actors.values() if isinstance(a, pv.Actor)]


def test_view3d_layers_share_a_viewer(grid3d):
    """Calls with viewer= add to the same scene, like ax= in the 2D plots."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, opacity=0.3, viewer=gs.Viewer(off_screen=True))
    same = gs.view3d(gs.ExampleData('3d_estimate'), 'Estimate', viewer=viewer)
    assert same is viewer
    blocks, points = layers(viewer)
    assert isinstance(blocks.mapper.dataset, pv.ImageData)
    assert blocks.prop.opacity == 0.3
    assert points.mapper.dataset.n_points == 100


def test_view3d_colours(grid3d):
    """clim sets the colour range, the default colormap is gs.Parameters['plotting.cmap'] and the
    colour bar is titled by label."""
    from matplotlib import colormaps
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, clim=(0, 2), label='Grade',
                       viewer=gs.Viewer(off_screen=True))
    (actor,) = layers(viewer)
    assert tuple(actor.mapper.scalar_range) == (0, 2)
    top = [0.0, 0.0, 0.0]
    actor.mapper.lookup_table.GetColor(2, top)
    cmap = colormaps[gs.Parameters['plotting.cmap']]
    np.testing.assert_allclose(top, cmap(1.0)[:3], atol=0.01)
    assert 'Grade' in viewer.plotter.scalar_bars


def test_view3d_mesh_and_colour(subblocks):
    """A PyVista mesh is shown as given; without var it is drawn in one colour."""
    mesh = gs.block_mesh(subblocks, x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))
    viewer = gs.view3d(mesh, color='red', viewer=gs.Viewer(off_screen=True))
    (actor,) = layers(viewer)
    assert actor.mapper.dataset is mesh
    assert actor.prop.color == pv.Color('red')
    assert not actor.mapper.scalar_visibility


def test_show_static(grid3d):
    """The scene is drawn on a matplotlib axis at the axis' own pixel size, from a viewer that was
    never shown on screen."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef)
    fig, ax = viewer.show_static(figsize=(4, 3), dpi=100, title='Blocks')
    (image,) = ax.get_images()
    box = ax.get_position()  # in fractions of the 4 x 3 inch figure
    assert image.get_array().shape[:2] == (round(box.height * 3 * 100), round(box.width * 4 * 100))
    assert image.get_array().std() > 0
    assert ax.get_title() == 'Blocks'
    assert fig is ax.figure


def test_show_static_on_given_axis(grid3d):
    """With ax=, the picture goes on that axis, next to whatever else is in the figure."""
    from matplotlib import pyplot as plt
    data, griddef = grid3d
    fig, (left, right) = plt.subplots(1, 2, figsize=(8, 4))
    out_fig, out_ax = gs.view3d(data, 'value', griddef=griddef).show_static(ax=right)
    assert (out_fig, out_ax) == (fig, right)
    assert len(right.get_images()) == 1 and len(left.get_images()) == 0


def test_show_static_colour_bar(grid3d):
    """The colour bar is a matplotlib one with the layer's label, limits and colormap, the picture
    still fills the (narrower) axis at full resolution, and cbar=False leaves it out."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, clim=(-1, 2), cmap='magma', label='Grade')
    fig, ax = viewer.show_static(figsize=(4, 3), dpi=100)
    (cbar_ax,) = [a for a in fig.axes if a is not ax]
    assert cbar_ax.get_ylabel() == 'Grade'
    assert cbar_ax.get_ylim() == (-1, 2)
    assert cbar_ax._colorbar.cmap.name == 'magma'
    box = ax.get_position()
    assert ax.get_images()[0].get_array().shape[:2] == (round(box.height * 300),
                                                        round(box.width * 400))
    fig, ax = viewer.show_static(figsize=(4, 3), cbar=False)
    assert fig.axes == [ax]


@pytest.mark.parametrize('change', ['none', 'zoom', 'parallel'])
def test_show_static_keeps_the_camera(subblocks, change):
    """show_static draws exactly what the live scene shows, zoom and projection included."""
    from matplotlib import pyplot as plt
    mesh = gs.block_mesh(subblocks, x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))
    viewer = gs.view3d(mesh, color='tan', viewer=gs.Viewer(off_screen=True))
    viewer.plotter.window_size = [400, 300]
    viewer.plotter.camera_position = 'iso'
    if change == 'zoom':
        viewer.plotter.camera.zoom(2)
    elif change == 'parallel':
        viewer.plotter.camera.parallel_projection = True
    live = viewer.plotter.screenshot(return_img=True)
    fig = plt.figure(figsize=(4, 3), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    viewer.show_static(ax=ax)
    np.testing.assert_array_equal(ax.get_images()[0].get_array(), live)


def test_camera_round_trip(subblocks):
    """A view saved as JSON and set on a new scene gives the same picture, live and static."""
    import json
    from matplotlib import pyplot as plt
    mesh = gs.block_mesh(subblocks, x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))

    def scene():
        viewer = gs.view3d(mesh, color='tan', viewer=gs.Viewer(off_screen=True))
        viewer.plotter.window_size = [400, 300]
        return viewer

    first = scene()
    first.plotter.camera_position = 'yz'
    first.plotter.camera.azimuth = 30
    first.plotter.camera.parallel_projection = True
    first.plotter.camera.zoom(1.5)
    saved = json.loads(json.dumps(first.get_camera()))
    assert saved['parallel_projection'] is True

    second = scene()
    second.set_camera(saved)
    assert second.get_camera() == saved
    expected = first.plotter.screenshot(return_img=True)
    np.testing.assert_array_equal(second.plotter.screenshot(return_img=True), expected)
    # A scene that was never drawn on screen still uses the camera in show_static
    third = scene()
    third.set_camera(saved)
    ax = plt.figure(figsize=(4, 3), dpi=100).add_axes([0, 0, 1, 1])
    third.show_static(ax=ax)
    np.testing.assert_array_equal(ax.get_images()[0].get_array(), expected)


def test_ruler_titles(subblocks):
    """Axis titles follow gs.Parameters like the 2D plots; unit='' leaves the unit out."""
    mesh = gs.block_mesh(subblocks, x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))
    viewer = gs.view3d(mesh, viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    unit = gs.Parameters['plotting.unit']
    axes = viewer.plotter.renderer.cube_axes_actor
    assert axes.GetXTitle() == '{} ({})'.format(gs.Parameters['plotting.xname'], unit)
    assert axes.GetZTitle() == '{} ({})'.format(gs.Parameters['plotting.zname'], unit)
    viewer.add_ruler(unit='')
    assert viewer.plotter.renderer.cube_axes_actor.GetYTitle() == gs.Parameters['plotting.yname']


def test_ruler_box_does_not_grow(subblocks):
    """The ruler box is the data plus 2% (not below 0 where the data starts at 0), and drawing it
    again does not make it bigger."""
    mesh = gs.block_mesh(subblocks, x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))
    viewer = gs.view3d(mesh, viewer=gs.Viewer(off_screen=True))
    for _ in range(3):
        viewer.add_ruler()
        np.testing.assert_allclose(viewer.plotter.renderer.cube_axes_actor.bounds,
                                   (0, 20.4, 0, 10.2, 0, 10.2))


def test_show_static_draws_the_ruler(subblocks):
    """show_static includes the ruler: the picture matches the live view pixel for pixel."""
    from matplotlib import pyplot as plt
    mesh = gs.block_mesh(subblocks, x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))
    viewer = gs.view3d(mesh, color='tan', viewer=gs.Viewer(off_screen=True))
    viewer.plotter.window_size = [400, 300]
    viewer.add_ruler()
    viewer.plotter.camera_position = 'iso'
    live = viewer.plotter.screenshot(return_img=True)
    ax = plt.figure(figsize=(4, 3), dpi=100).add_axes([0, 0, 1, 1])
    viewer.show_static(ax=ax)
    np.testing.assert_array_equal(ax.get_images()[0].get_array(), live)
