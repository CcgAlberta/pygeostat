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
    """Drawing the ruler again does not make its box bigger."""
    mesh = gs.block_mesh(subblocks, x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))
    viewer = gs.view3d(mesh, viewer=gs.Viewer(off_screen=True))
    for _ in range(3):
        viewer.add_ruler()
        np.testing.assert_allclose(viewer.plotter.renderer.cube_axes_actor.bounds,
                                   (0, 20, 0, 10, 0, 10))


@pytest.mark.parametrize('bounds, expected', [
    ((0, 80, 0, 80, 0, 40), [['0', '20', '40', '60', '80'], ['0', '20', '40', '60', '80'],
                             ['0', '20', '40']]),
    ((612950, 613420, 5055100, 5061000, 146, 286),
     [['612500', '613000', '613500'], ['5054000', '5056000', '5058000', '5060000', '5062000'],
      ['100', '200', '300']]),
    ((0.12, 0.47, 0, 0.3, -3, -1), [['0.0', '0.2', '0.4', '0.6'], ['0.0', '0.2', '0.4'],
                                    ['-3.0', '-2.5', '-2.0', '-1.5', '-1.0']]),
])
def test_ruler_round_labels(bounds, expected):
    """Labels land on round numbers, each axis with only the decimals its own step needs."""
    viewer = gs.view3d(pv.Box(bounds), viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    viewer.plotter.render()
    axes = viewer.plotter.renderer.cube_axes_actor
    labels = [[axes.GetAxisLabels(i).GetValue(j) for j in range(axes.GetAxisLabels(i).GetNumberOfValues())]
              for i in range(3)]
    assert labels == expected


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


def ruler_labels(viewer, axis):
    """The label strings the ruler draws on one axis (0, 1, 2 for x, y, z)."""
    viewer.plotter.render()
    labels = viewer.plotter.renderer.cube_axes_actor.GetAxisLabels(axis)
    return [labels.GetValue(i) for i in range(labels.GetNumberOfValues())]


def test_vertical_exaggeration():
    """The scene is stretched vertically, the ruler wraps the stretched scene but keeps the real
    elevations, and a layer added afterwards is stretched too."""
    viewer = gs.view3d(pv.Box((0, 100, 0, 100, 0, 10)), viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    viewer.set_vertical_exaggeration(5)
    assert viewer.plotter.bounds[4:] == (0, 50)
    assert ruler_labels(viewer, 2) == ['0', '5', '10']
    assert viewer.plotter.renderer.cube_axes_actor.bounds[4:] == (0, 50)
    gs.view3d(pv.Box((0, 50, 0, 50, 10, 20)), viewer=viewer)
    assert viewer.plotter.bounds[4:] == (0, 100)
    assert ruler_labels(viewer, 2) == ['0', '5', '10', '15', '20']


def test_show_static_vertical_exaggeration(subblocks):
    """show_static draws the stretched scene and ruler exactly as the live view."""
    from matplotlib import pyplot as plt
    mesh = gs.block_mesh(subblocks, x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))
    viewer = gs.view3d(mesh, color='tan', viewer=gs.Viewer(off_screen=True))
    viewer.plotter.window_size = [400, 300]
    viewer.add_ruler()
    viewer.set_vertical_exaggeration(3)
    viewer.plotter.camera_position = 'iso'
    live = viewer.plotter.screenshot(return_img=True)
    ax = plt.figure(figsize=(4, 3), dpi=100).add_axes([0, 0, 1, 1])
    viewer.show_static(ax=ax)
    np.testing.assert_array_equal(ax.get_images()[0].get_array(), live)


def test_ruler_follows_new_layers():
    """A layer added after the ruler is inside the ruler."""
    viewer = gs.view3d(pv.Box((0, 100, 0, 100, 0, 10)), viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    gs.view3d(pv.Box((0, 100, 0, 100, 10, 40)), viewer=viewer)
    assert ruler_labels(viewer, 2) == ['0', '20', '40']


def test_view3d_threshold(grid3d):
    """Only blocks with a value inside the threshold are shown (count checked against the data
    directly), and the colours keep the full range of the variable."""
    data, griddef = grid3d
    values = data.data['value'].to_numpy()
    viewer = gs.view3d(data, 'value', griddef=griddef, threshold=(1.0, 2.0),
                       viewer=gs.Viewer(off_screen=True))
    (actor,) = layers(viewer)
    shown = actor.mapper.dataset.cell_data['value']
    assert len(shown) == np.sum((values >= 1.0) & (values <= 2.0))
    assert shown.min() >= 1.0 and shown.max() <= 2.0
    np.testing.assert_allclose(actor.mapper.scalar_range, (values.min(), values.max()))


def test_view3d_threshold_drillholes(intervals):
    """Intervals without a value are hidden by a threshold, and it needs a variable."""
    mesh = gs.drillhole_mesh(intervals, var='Au', start=('startx', 'starty', 'startz'),
                             end=('endx', 'endy', 'endz'))
    viewer = gs.view3d(mesh, 'Au', threshold=(0, 10), viewer=gs.Viewer(off_screen=True))
    (actor,) = layers(viewer)
    assert actor.mapper.dataset.bounds[4:] == (90.0, 100.0)
    with pytest.raises(ValueError, match='needs a variable'):
        gs.view3d(mesh, threshold=(0, 10), viewer=gs.Viewer(off_screen=True))


def test_ruler_keeps_the_whole_model_under_a_threshold(grid3d):
    """With a high threshold only a few blocks are left, but the ruler still spans the whole
    model, so they can be placed within it."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, threshold=(3.5, 10),
                       viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    (actor,) = layers(viewer)
    assert actor.mapper.dataset.n_cells == np.sum(data.data['value'].to_numpy() >= 3.5) < 20
    np.testing.assert_allclose(viewer.plotter.renderer.cube_axes_actor.bounds,
                               (0, 80, 0, 80, 0, 40))


@pytest.mark.parametrize('orient, thickness', [('xy', 0), ('xy', 2), ('xz', 1), ('yz', 0)])
def test_view3d_slice_grid(grid3d, orient, thickness):
    """A grid slice shows grid layers slice_number - thickness to slice_number + thickness,
    checked against the data reshaped in grid order (x fastest, then y, then z)."""
    data, griddef = grid3d
    cube = data.data['value'].to_numpy().reshape(griddef.nz, griddef.ny, griddef.nx)
    n = 10
    layers_kept = {'xy': cube[n - thickness:n + thickness + 1],
                   'xz': cube[:, n - thickness:n + thickness + 1],
                   'yz': cube[:, :, n - thickness:n + thickness + 1]}[orient]
    viewer = gs.view3d(data, 'value', griddef=griddef, orient=orient, slice_number=n,
                       slice_thickness=thickness, viewer=gs.Viewer(off_screen=True))
    (actor,) = layers(viewer)
    np.testing.assert_array_equal(np.sort(actor.mapper.dataset.cell_data['value']),
                                  np.sort(layers_kept.ravel()))


def test_view3d_slice_drillholes(intervals):
    """Drill holes are cut at pointtol either side of the slice coordinate."""
    viewer = gs.view3d(intervals, 'Au', start=('startx', 'starty', 'startz'),
                       end=('endx', 'endy', 'endz'), orient='xy', slice_coordinate=92.0,
                       pointtol=3.0, viewer=gs.Viewer(off_screen=True))
    (actor,) = layers(viewer)
    assert actor.mapper.dataset.bounds[4:] == pytest.approx((89.0, 95.0))


def test_view3d_slice_points():
    """Points within pointtol of the slice are kept, counted directly from the data."""
    data = gs.ExampleData('oilsands')
    north = data.data[data.y].to_numpy()
    centre = float(np.median(north))
    viewer = gs.view3d(data, 'Bitumen', orient='xz', slice_coordinate=centre, pointtol=100.0,
                       viewer=gs.Viewer(off_screen=True))
    (actor,) = layers(viewer)
    assert actor.mapper.dataset.n_points == np.sum(np.abs(north - centre) <= 100.0)


def test_view3d_slice_errors(grid3d, intervals):
    """A slice needs a position, data without a grid needs pointtol, and orient is checked."""
    data, griddef = grid3d
    tubes = dict(start=('startx', 'starty', 'startz'), end=('endx', 'endy', 'endz'),
                 viewer=gs.Viewer(off_screen=True))
    with pytest.raises(ValueError, match='slice_coordinate'):
        gs.view3d(intervals, 'Au', orient='xy', slice_number=3, **tubes)
    with pytest.raises(ValueError, match='pointtol'):
        gs.view3d(intervals, 'Au', orient='xy', slice_coordinate=90.0, **tubes)
    with pytest.raises(ValueError, match="'xy', 'xz' or 'yz'"):
        gs.view3d(data, 'value', griddef=griddef, orient='zz', slice_number=3,
                  viewer=gs.Viewer(off_screen=True))


def test_view3d_slice_keeps_full_colour_range(grid3d):
    """A slice is coloured on the full range of the variable, so the colours don't change as the
    slice moves."""
    data, griddef = grid3d
    values = data.data['value'].to_numpy()
    for n in (5, 30):
        viewer = gs.view3d(data, 'value', griddef=griddef, orient='xy', slice_number=n,
                           viewer=gs.Viewer(off_screen=True))
        (actor,) = layers(viewer)
        np.testing.assert_allclose(actor.mapper.scalar_range, (values.min(), values.max()))


def colours_2d(data, var, **kwargs):
    """The colour location_plot draws for each code of var, read off its scatter points."""
    from matplotlib import pyplot as plt
    ax = gs.location_plot(data, x='x', y='y', var=var, **kwargs)
    faces = [c for c in ax.collections if len(c.get_facecolors()) == len(data)][0].get_facecolors()
    plt.close(ax.figure)
    return {int(code): tuple(np.round(faces[i, :3], 3)) for i, code in enumerate(data[var])}


def colours_3d(viewer):
    """The colour of each category in the 3D scene, by its label."""
    (actor,) = layers(viewer)
    table = actor.mapper.lookup_table
    names = viewer._layers[0].kwargs['annotations']
    colours = {}
    for index, name in names.items():
        rgb = [0.0, 0.0, 0.0]
        table.GetColor(index, rgb)
        colours[name] = tuple(np.round(rgb, 3))
    return colours


def test_view3d_categories_match_location_plot():
    """With a catdict that names exactly the codes in the data, each category gets the colour
    location_plot gives it."""
    data = pd.DataFrame({'x': [0.0, 1.0, 2.0, 3.0], 'y': [0.0, 1.0, 2.0, 3.0], 'z': [0.0] * 4,
                         'Facies': [3, 1, 3, 1]})
    catdict = {1: 'Sand', 3: 'Shale'}
    expected = colours_2d(data, 'Facies', catdata=True, catdict=catdict)
    got = colours_3d(gs.view3d(data, 'Facies', x='x', y='y', z='z', catdict=catdict,
                               viewer=gs.Viewer(off_screen=True)))
    assert got == {'Sand': expected[1], 'Shale': expected[3]}


def test_view3d_categories_keep_catdict_colours():
    """When catdict also names codes that are not in the data, each category keeps the colour of
    its place in catdict, so a rock type has the same colour whatever subset is shown.
    (location_plot means to do the same but draws every point in one colour in this case.)"""
    from pygeostat.plotting.utils import _get_cmap
    data = pd.DataFrame({'x': [0.0, 1.0], 'y': [0.0, 1.0], 'z': [0.0, 0.0], 'Facies': [3, 1]})
    palette = _get_cmap(None, True, 3).colors
    got = colours_3d(gs.view3d(data, 'Facies', x='x', y='y', z='z',
                               catdict={1: 'Sand', 2: 'Silt', 3: 'Shale'},
                               viewer=gs.Viewer(off_screen=True)))
    assert got == {'Sand': tuple(np.round(palette[0][:3], 3)),
                   'Shale': tuple(np.round(palette[2][:3], 3))}


def test_view3d_categories_from_datafile():
    """A DataFile's cat column is coloured by category without passing anything, and
    show_static labels the colour bar with the categories."""
    data = gs.ExampleData('point3d_ind_mv')
    viewer = gs.view3d(data, data.cat)
    assert viewer._layers[0].kwargs['annotations'] == {0: '1', 1: '2', 2: '3', 3: '4', 4: '5'}
    fig, ax = viewer.show_static(figsize=(4, 3))
    (bar,) = [a for a in fig.axes if a is not ax]
    assert [t.get_text() for t in bar.get_yticklabels()] == ['1', '2', '3', '4', '5']


def test_view3d_category_missing_from_catdict(intervals):
    """A code the catdict doesn't name gives a clear error."""
    intervals['Rock'] = [1, 2]
    with pytest.raises(ValueError, match=r"Codes \[2\] of 'Rock' are not in catdict"):
        gs.view3d(intervals, 'Rock', start=('startx', 'starty', 'startz'),
                  end=('endx', 'endy', 'endz'), catdict={1: 'Granite'},
                  viewer=gs.Viewer(off_screen=True))


def test_ruler_short_axis_gets_fewer_labels():
    """A short axis on screen has fewer labels, so they don't run together, and gets more as
    vertical exaggeration makes it longer on screen."""
    viewer = gs.view3d(pv.Box((0, 1000, 0, 1000, 0, 40)), viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    assert ruler_labels(viewer, 0) == ['0', '500', '1000']
    assert ruler_labels(viewer, 2) == ['0', '20', '40']
    viewer.set_vertical_exaggeration(25)
    assert ruler_labels(viewer, 2) == ['0', '10', '20', '30', '40']


def test_layers_are_named_and_numbered(grid3d):
    """A layer is named by its label or variable, a repeated name gets a number, and a layer
    without a variable is 'Layer N'."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    gs.view3d(data, 'value', griddef=griddef, viewer=viewer)
    gs.view3d(data, 'value', griddef=griddef, label='Copy', viewer=viewer)
    viewer.add(pv.Sphere(), color='red')
    assert viewer.layers == ['value', 'value (2)', 'Copy', 'Layer 4']


def test_set_layer_hides_and_changes_opacity(grid3d):
    """set_layer changes the actor of the named layer only."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    gs.view3d(gs.ExampleData('3d_estimate'), 'Estimate', viewer=viewer)
    blocks, points = layers(viewer)
    viewer.set_layer('value', visible=False, opacity=0.4)
    assert not blocks.visibility and points.visibility
    assert blocks.prop.opacity == 0.4 and points.prop.opacity == 1
    viewer.set_layer('value', visible=True)
    assert blocks.visibility and blocks.prop.opacity == 0.4


def test_set_layer_errors(grid3d):
    """An unknown layer names the layers there are; opacity outside 0 to 1 is refused."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    with pytest.raises(ValueError, match=r"No layer 'x'.*\['value'\]"):
        viewer.set_layer('x', visible=False)
    with pytest.raises(ValueError, match='opacity'):
        viewer.set_layer('value', opacity=1.5)
    assert layers(viewer)[0].prop.opacity == 1


def test_show_static_follows_set_layer(grid3d):
    """A hidden layer is left out of the picture, with its colour bar, and the opacity is kept."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    gs.view3d(gs.ExampleData('3d_estimate'), 'Estimate', viewer=viewer)
    fig, ax = viewer.show_static()
    both, bars = np.asarray(ax.images[0].get_array()), len(fig.axes)
    viewer.set_layer('Estimate', visible=False)
    fig, ax = viewer.show_static()
    assert len(fig.axes) == bars - 1
    assert not np.array_equal(both, np.asarray(ax.images[0].get_array()))
    viewer.set_layer('Estimate', visible=True)
    fig, ax = viewer.show_static()
    opaque = np.asarray(ax.images[0].get_array())
    viewer.set_layer('value', opacity=0.1)
    fig, ax = viewer.show_static()
    assert not np.array_equal(opaque, np.asarray(ax.images[0].get_array()))


def test_set_layer_clim(grid3d):
    """clim changes the colour range of the actor and of the colour bar in show_static; it is
    refused for a layer without a continuous variable, and a refused call changes nothing."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    viewer.add(pv.Sphere(), color='red')
    gs.view3d(gs.ExampleData('point3d_ind_mv'), 'Lithofacies', catdict={1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e'},
              viewer=viewer)
    viewer.set_layer('value', clim=(0.5, 1.5))
    assert tuple(layers(viewer)[0].mapper.scalar_range) == (0.5, 1.5)
    fig, ax = viewer.show_static()
    assert fig.axes[1].get_ylim() == (0.5, 1.5)
    for name in ('Layer 2', 'Lithofacies'):
        with pytest.raises(ValueError, match='continuous'):
            viewer.set_layer(name, opacity=0.5, clim=(0, 1))
    with pytest.raises(ValueError, match='min below max'):
        viewer.set_layer('value', opacity=0.5, clim=(2, 1))
    assert layers(viewer)[0].prop.opacity == 1 and layers(viewer)[1].prop.opacity == 1


def test_set_layer_threshold(grid3d):
    """A threshold set afterwards gives the same layer as view3d(threshold=), also inside a
    slice, keeps the colours, and False shows the whole layer again."""
    data, griddef = grid3d
    options = dict(griddef=griddef, orient='xz', slice_number=20)
    viewer = gs.view3d(data, 'value', viewer=gs.Viewer(off_screen=True), opacity=0.4, **options)
    whole = layers(viewer)[0].mapper.dataset.n_cells
    viewer.set_layer('value', threshold=(0.5, 1.0))
    expected = gs.view3d(data, 'value', threshold=(0.5, 1.0), viewer=gs.Viewer(off_screen=True), **options)
    (actor,) = layers(viewer)
    assert 0 < actor.mapper.dataset.n_cells < whole
    assert actor.mapper.dataset.n_cells == layers(expected)[0].mapper.dataset.n_cells
    assert actor.prop.opacity == 0.4
    assert tuple(actor.mapper.scalar_range) == tuple(data.data['value'].agg(['min', 'max']))
    assert list(viewer.plotter.scalar_bars.keys()) == ['value']
    viewer.set_layer('value', threshold=False)
    assert layers(viewer)[0].mapper.dataset.n_cells == whole


def test_set_layer_threshold_errors(grid3d):
    """A threshold that is reversed, leaves nothing, or is on a plain mesh is refused and changes
    nothing."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    viewer.add(pv.Sphere(), color='red')
    with pytest.raises(ValueError, match='min not above max'):
        viewer.set_layer('value', threshold=(2, 1))
    with pytest.raises(ValueError, match='Nothing is left'):
        viewer.set_layer('value', opacity=0.5, threshold=(100, 200))
    with pytest.raises(ValueError, match='added as a mesh'):
        viewer.set_layer('Layer 2', threshold=(0, 1))
    assert layers(viewer)[0].mapper.dataset.n_cells == 40 ** 3
    assert layers(viewer)[0].prop.opacity == 1


def test_set_layer_slice(grid3d):
    """A slice set afterwards gives the same layer as view3d(orient=), moving it needs only what
    changes, the threshold stays, and orient=False shows the whole layer."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, threshold=(0.5, 1.0),
                       viewer=gs.Viewer(off_screen=True))
    whole = layers(viewer)[0].mapper.dataset.n_cells
    viewer.set_layer('value', orient='xz', slice_number=20)
    first = layers(viewer)[0].mapper.dataset.n_cells
    expected = gs.view3d(data, 'value', griddef=griddef, threshold=(0.5, 1.0), orient='xz', slice_number=20,
                         viewer=gs.Viewer(off_screen=True))
    assert first == layers(expected)[0].mapper.dataset.n_cells < whole
    viewer.set_layer('value', slice_number=21, slice_thickness=2)
    moved = gs.view3d(data, 'value', griddef=griddef, threshold=(0.5, 1.0), orient='xz', slice_number=21,
                      slice_thickness=2, viewer=gs.Viewer(off_screen=True))
    assert layers(viewer)[0].mapper.dataset.n_cells == layers(moved)[0].mapper.dataset.n_cells
    viewer.set_layer('value', orient=False)
    assert layers(viewer)[0].mapper.dataset.n_cells == whole


def test_set_layer_slice_by_coordinate():
    """Drill holes are sliced by coordinate and distance, and a position given last replaces the
    other kind."""
    data = gs.ExampleData('oilsands')
    viewer = gs.view3d(data, 'Bitumen', viewer=gs.Viewer(off_screen=True))
    viewer.set_layer('Bitumen', orient='xz', slice_coordinate=data.data[data.y].median(), pointtol=500)
    n = layers(viewer)[0].mapper.dataset.n_points
    assert 0 < n < len(data)
    viewer.set_layer('Bitumen', pointtol=1000)
    assert layers(viewer)[0].mapper.dataset.n_points > n


def test_set_layer_slice_errors(grid3d):
    """A slice that is not defined, wrong, or leaves nothing is refused and changes nothing."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    with pytest.raises(ValueError, match='Pass `orient`'):
        viewer.set_layer('value', slice_number=3)
    with pytest.raises(ValueError, match='orient'):
        viewer.set_layer('value', orient='ab', slice_number=3)
    with pytest.raises(ValueError, match='Nothing is left'):
        viewer.set_layer('value', opacity=0.5, orient='xz', slice_coordinate=1e6)
    assert layers(viewer)[0].mapper.dataset.n_cells == 40 ** 3
    assert layers(viewer)[0].prop.opacity == 1
    assert viewer._layers[0].section is None


def far_box(viewer=None):
    """A box far from the origin, 100 by 100 by 40, in a scene with a ruler."""
    viewer = viewer or gs.Viewer(off_screen=True)
    viewer.add(pv.Box(bounds=(1000, 1100, 5000, 5100, 300, 340)), color='red')
    viewer.add_ruler()
    return viewer


def test_axis_origin_counts_labels_from_the_origin():
    """With an origin the labels are those of the same box moved to the origin, the box stays where
    it is, and the titles say where the labels start. None gives the coordinates back."""
    near = gs.Viewer(off_screen=True)
    near.add(pv.Box(bounds=(0, 100, 0, 100, 0, 40)), color='red')
    near.add_ruler()
    viewer = far_box()
    far_labels = [ruler_labels(viewer, axis) for axis in range(3)]
    viewer.set_axis_origin((1000, 5000, 300))
    assert [ruler_labels(viewer, axis) for axis in range(3)] == [ruler_labels(near, axis) for axis in range(3)]
    assert far_labels != [ruler_labels(viewer, axis) for axis in range(3)]
    axes = viewer.plotter.renderer.cube_axes_actor
    np.testing.assert_allclose(axes.bounds, (1000, 1100, 5000, 5100, 300, 340))
    assert axes.GetXTitle().endswith(' from 1000') and axes.GetZTitle().endswith(' from 300')
    viewer.set_axis_origin(None)
    assert [ruler_labels(viewer, axis) for axis in range(3)] == far_labels
    assert ' from ' not in viewer.plotter.renderer.cube_axes_actor.GetXTitle()


def test_axis_origin_from_zero_and_new_layers():
    """from_zero starts each axis at the lowest edge of the scene, also after a layer is added,
    and show_static draws the same."""
    near = gs.Viewer(off_screen=True)
    near.add(pv.Box(bounds=(0, 100, 0, 100, 0, 40)), color='red')
    near.add_ruler()
    viewer = far_box()
    viewer.set_axis_origin(from_zero=True)
    assert [ruler_labels(viewer, axis) for axis in range(3)] == [ruler_labels(near, axis) for axis in range(3)]
    viewer.add(pv.Box(bounds=(1000, 1100, 4900, 5100, 300, 340)), color='blue')
    assert ruler_labels(viewer, 1)[0] == '0'
    assert viewer.plotter.renderer.cube_axes_actor.bounds[2] == 4900
    fig, ax = viewer.show_static()
    assert ax.images


def test_axis_origin_needs_three_numbers():
    with pytest.raises(ValueError, match='easting, northing, elevation'):
        gs.Viewer(off_screen=True).set_axis_origin((1, 2))


def test_axis_origin_from_zero_without_layers():
    """Before any layer there is nothing to count from, so the origin is 0 and nothing fails."""
    viewer = gs.Viewer(off_screen=True)
    viewer.set_axis_origin(from_zero=True)
    assert viewer._origin_values() == (0.0, 0.0, 0.0)


def two_layer_viewer(grid3d):
    """A block model and drill holes with a ruler, to set up and compare views."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    gs.view3d(gs.ExampleData('oilsands'), 'Bitumen', trace=True, radius=15, viewer=viewer)
    viewer.add_ruler()
    return viewer


def test_view_round_trip(grid3d):
    """get_view survives json and gives the same scene back on a fresh viewer: the layers'
    visibility, opacity, colour limits, threshold and slice, the ruler origin, the exaggeration,
    the background, the ruler's ticks, what is drawn around the layers and the camera."""
    import json
    viewer = two_layer_viewer(grid3d)
    viewer.set_layer('value', opacity=0.4, clim=(-1, 1), threshold=(-2, 2), orient='xz', slice_number=7,
                     slice_thickness=2)
    viewer.set_layer('Bitumen', visible=False)
    viewer.set_axis_origin((1, 2, 3))
    viewer.set_vertical_exaggeration(4)
    viewer.set_background('black')
    viewer.add_ruler()
    viewer.set_axis_ticks(xy=20)
    viewer.set_scene(outline=True, axis_arrows=True, corner_axes=True)
    viewer.plotter.camera.position = (300, -400, 200)
    view = json.loads(json.dumps(viewer.get_view()))
    fresh = two_layer_viewer(grid3d)
    fresh.set_view(view)
    assert fresh.get_view() == view
    assert fresh._layers[0].mesh.n_cells == viewer._layers[0].mesh.n_cells
    assert not fresh._layers[1].actor.visibility
    assert fresh._ruler is not None and fresh._tick_spacing == (20, None) and 'outline' in fresh.plotter.renderer.actors


def test_set_view_with_other_layers(grid3d):
    """A view that names a layer the scene lacks is refused and changes nothing; a scene with more
    layers than the view keeps them as they are."""
    viewer = two_layer_viewer(grid3d)
    view = viewer.get_view()
    other = gs.view3d(grid3d[0], 'value', griddef=grid3d[1], opacity=0.3, viewer=gs.Viewer(off_screen=True))
    with pytest.raises(ValueError, match="'Bitumen'"):
        other.set_view(view)
    assert other._layers[0].kwargs['opacity'] == 0.3
    view['layers'] = view['layers'][:1]
    view['layers'][0]['opacity'] = 0.2
    viewer.set_layer('Bitumen', opacity=0.6)
    viewer.set_view(view)
    assert viewer._layers[0].kwargs['opacity'] == 0.2 and viewer._layers[1].kwargs['opacity'] == 0.6


def test_set_view_that_fails_midway_changes_nothing(grid3d):
    """A view whose second layer cannot be shown (the threshold leaves nothing) is refused after the
    first layer was changed, and the first layer goes back."""
    viewer = two_layer_viewer(grid3d)
    view = viewer.get_view()
    view['layers'][0]['opacity'] = 0.2
    view['layers'][1]['threshold'] = [1e9, 2e9]
    with pytest.raises(ValueError, match='Nothing is left'):
        viewer.set_view(view)
    assert viewer.get_view()['layers'][0]['opacity'] == 1.0


def test_drillhole_mesh_keeps_what_a_slice_needs():
    """Tubes carry their interval, drill hole, ends and radius, so a plane can find the intervals
    it crosses: with the ends of an interval, its hole, and for each cell the interval it is part of."""
    from pygeostat.plotting3d.meshes import SEGMENT, HOLE, UPPER, LOWER, RADIUS
    data = pd.DataFrame({'x': [0.0, 0.0, 0.0, 50.0, 50.0], 'y': 0.0, 'z': [30.0, 20.0, 10.0, 25.0, 15.0],
                         'hole': ['a', 'a', 'a', 'b', 'b'], 'Au': [1.0, 2.0, 3.0, 4.0, 5.0]})
    mesh = gs.drillhole_mesh(data, var='Au', x='x', y='y', z='z', dh='hole', trace=True, radius=3)
    np.testing.assert_array_equal(mesh.field_data[UPPER][:, 2], [30, 20, 25])
    np.testing.assert_array_equal(mesh.field_data[LOWER][:, 2], [20, 10, 15])
    np.testing.assert_array_equal(mesh.field_data[RADIUS], [3])
    segment, hole = mesh.cell_data[SEGMENT], mesh.cell_data[HOLE]
    assert set(segment) == {0, 1, 2}
    # Each cell's value is the upper sample of its interval; intervals 0 and 1 are hole a, 2 is hole b
    np.testing.assert_array_equal(mesh.cell_data['Au'], np.array([1.0, 2.0, 4.0])[segment])
    hole_of = {s: set(hole[segment == s]) for s in (0, 1, 2)}
    assert hole_of[0] == hole_of[1] and len(hole_of[0]) == 1 and hole_of[2] != hole_of[0]


def test_drillhole_mesh_intervals_without_a_hole_id_are_their_own_hole(intervals):
    from pygeostat.plotting3d.meshes import SEGMENT, HOLE
    mesh = gs.drillhole_mesh(intervals, var='Au', start=('startx', 'starty', 'startz'),
                             end=('endx', 'endy', 'endz'))
    np.testing.assert_array_equal(np.unique(mesh.cell_data[HOLE]), [0, 1])
    np.testing.assert_array_equal(mesh.cell_data[HOLE], mesh.cell_data[SEGMENT])


def test_slicer_cuts_a_block_layer(grid3d):
    """The slicer shows the cut through its plane, with the values of the blocks the plane runs
    through (the same as the slice of those blocks), and off shows the whole layer again."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    viewer.set_slicer(enabled=True, origin=(40, 40, 10.25), normal=(0, 0, 1))
    cut = viewer._layers[0].mesh
    np.testing.assert_allclose(cut.points[:, 2], 10.25)
    expected = gs.view3d(data, 'value', griddef=griddef, orient='xy', slice_number=10,
                         viewer=gs.Viewer(off_screen=True))._layers[0].mesh
    assert cut.n_cells == expected.n_cells == 40 * 40
    np.testing.assert_array_equal(np.sort(cut.cell_data['value']), np.sort(expected.cell_data['value']))
    viewer.set_slicer(enabled=False)
    assert viewer._layers[0].mesh.n_cells == 40 ** 3


def test_slicer_plane_can_be_tilted(grid3d):
    """A tilted plane gives a cut whose points all lie on the plane."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    viewer.set_slicer(enabled=True, origin=(40, 40, 20), normal=(1, 1, 0))
    points = viewer._layers[0].mesh.points
    assert len(points)
    np.testing.assert_allclose((points - (40, 40, 20)) @ (1, 1, 0), 0, atol=1e-9)


def test_slicer_defaults_and_errors(grid3d):
    """The plane starts at the centre of the scene facing up, a zero normal means up, and wrong
    arguments are refused."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    viewer.set_slicer(enabled=True)
    xmin, xmax, ymin, ymax, zmin, zmax = viewer._layers[0].source.mesh.bounds
    assert viewer._slicer['origin'] == ((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
    viewer.set_slicer(normal=(0, 0, 0))
    assert viewer._slicer['normal'] == (0.0, 0.0, 1.0)
    with pytest.raises(ValueError, match='x, y, z'):
        viewer.set_slicer(origin=(1, 2))
    with pytest.raises(ValueError, match='nx, ny, nz'):
        viewer.set_slicer(normal=(1, 2))
    with pytest.raises(ValueError, match='0 or more'):
        viewer.set_slicer(tolerance=-1)


def test_slicer_accepts_a_tolerance_of_zero(grid3d):
    """A tolerance of 0 is allowed (only a negative one is refused)."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    viewer.set_slicer(tolerance=0)
    assert viewer._slicer['tolerance'] == 0.0


def test_slicer_with_threshold_and_a_plane_that_misses(grid3d):
    """A threshold applies to the cut too; a plane outside a layer leaves nothing to draw, which is
    not an error, and the layer comes back when the plane does."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    viewer.set_slicer(enabled=True, origin=(40, 40, 10.25), normal=(0, 0, 1))
    viewer.set_layer('value', threshold=(0.5, 1.0))
    values = viewer._layers[0].mesh.cell_data['value']
    assert 0 < len(values) < 40 * 40 and values.min() >= 0.5 and values.max() <= 1.0
    viewer.set_slicer(origin=(40, 40, 1000))
    assert viewer._layers[0].actor is None and viewer._layers[0].mesh.n_cells == 0
    viewer.set_layer('value', opacity=0.5, clim=(0, 1), visible=True)
    fig, ax = viewer.show_static()
    assert ax.images
    viewer.set_slicer(origin=(40, 40, 10.25))
    assert viewer._layers[0].actor is not None
    assert viewer._layers[0].actor.prop.opacity == 0.5


def test_slicer_leaves_layers_made_from_a_mesh(grid3d):
    """A layer added as a ready mesh has nothing to cut from and stays as it is."""
    viewer = gs.Viewer(off_screen=True)
    viewer.add(pv.Sphere(radius=5), color='red')
    viewer.set_slicer(enabled=True, origin=(0, 0, 0))
    assert viewer._layers[0].mesh.n_cells == pv.Sphere(radius=5).n_cells


def test_slicer_leaves_samples_drawn_as_points():
    """Samples without intervals (points) are neither cut nor hidden by the slicer."""
    data = gs.ExampleData('oilsands')
    viewer = gs.view3d(data, 'Bitumen', viewer=gs.Viewer(off_screen=True))
    assert viewer._layers[0].kind == 'points'
    viewer.set_slicer(enabled=True)
    assert viewer._layers[0].mesh.n_points == len(data) and viewer._layers[0].actor is not None
    holes = gs.view3d(data, 'Bitumen', trace=True, radius=15, viewer=gs.Viewer(off_screen=True))
    assert holes._layers[0].kind == 'holes'


@pytest.fixture
def two_holes():
    """Two vertical drill holes with samples every 10 down to 0, at x = 0 and x = 50."""
    return pd.DataFrame({'x': [0.0] * 4 + [50.0] * 4, 'y': 0.0, 'z': [30.0, 20.0, 10.0, 0.0] * 2,
                         'hole': ['a'] * 4 + ['b'] * 4, 'Au': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]})


def hole_viewer(two_holes):
    return gs.view3d(two_holes, 'Au', x='x', y='y', z='z', dh='hole', trace=True, radius=3,
                     viewer=gs.Viewer(off_screen=True))


def test_slicer_hides_the_drill_holes_and_shows_them_again(two_holes, grid3d):
    """Switching the slicer on hides drill hole layers (blocks are cut instead), and switching it
    off shows them again."""
    data, griddef = grid3d
    viewer = hole_viewer(two_holes)
    gs.view3d(data, 'value', griddef=griddef, viewer=viewer)
    holes, blocks = viewer._layers
    viewer.set_slicer(enabled=True, origin=(10, 0, 10), normal=(0, 0, 1))
    assert not holes.visible and not holes.actor.visibility and blocks.visible
    viewer.set_slicer(tolerance=5)
    assert not holes.visible
    viewer.set_slicer(enabled=False)
    assert holes.visible and holes.actor.visibility


def test_slicer_shows_the_holes_that_cross_the_plane(two_holes):
    """Only holes whose tube reaches the plane are shown, with the intervals inside the tolerance
    band, cut to the band."""
    viewer = hole_viewer(two_holes)
    layer = viewer._layers[0]
    viewer.set_slicer(enabled=True, origin=(0, 0, 0), normal=(1, 0, 0), tolerance=4, show_holes=True)
    # The plane x = 0 runs along hole a (all three intervals) and misses hole b
    assert set(layer.crossing.cell_data['segment index']) == {0, 1, 2}
    assert layer.crossing.bounds[0] >= -2 - 1e-5 and layer.crossing.bounds[1] <= 2 + 1e-5
    # A plane across both holes within a thin band: the interval from z = 20 to 10 of each
    viewer.set_slicer(origin=(0, 0, 15), normal=(0, 0, 1), tolerance=2)
    assert set(layer.crossing.cell_data['segment index']) == {1, 4}
    assert layer.crossing.bounds[4] >= 14 - 1e-5 and layer.crossing.bounds[5] <= 16 + 1e-5
    assert layer.crossing_actor is not None
    # A hole that only passes near the plane, further than its radius, is left out however wide the band
    viewer.set_slicer(origin=(10, 0, 0), normal=(1, 0, 0), tolerance=100)
    assert layer.crossing is None and layer.crossing_actor is None


def test_slicer_holes_are_shown_while_the_plane_is_inside_the_tube(two_holes):
    """Hole a (radius 3) is at x = 0: a plane 2 away cuts its tube, a plane 4 away does not,
    whatever the tolerance."""
    viewer = hole_viewer(two_holes)
    layer = viewer._layers[0]
    viewer.set_slicer(enabled=True, origin=(2, 0, 15), normal=(1, 0, 0), tolerance=100, show_holes=True)
    assert layer.crossing is not None and set(layer.crossing.cell_data['segment index']) == {0, 1, 2}
    viewer.set_slicer(origin=(4, 0, 15))
    assert layer.crossing is None


def test_slicer_holes_follow_the_threshold_and_switch_off(two_holes):
    """The holes on the slice follow the layer's threshold, and leave with show_holes or the slicer."""
    viewer = hole_viewer(two_holes)
    layer = viewer._layers[0]
    viewer.set_layer('Au', threshold=(3, 8))
    viewer.set_slicer(enabled=True, origin=(0, 0, 15), normal=(0, 0, 1), tolerance=100, show_holes=True)
    assert layer.crossing.cell_data['Au'].min() >= 3
    viewer.set_slicer(show_holes=False)
    assert layer.crossing is None
    viewer.set_slicer(show_holes=True)
    assert layer.crossing is not None
    viewer.set_slicer(enabled=False)
    assert layer.crossing is None and layer.crossing_actor is None


def test_slicer_colour_bar_and_static_picture_follow_the_holes(two_holes):
    """A hidden layer loses its colour bar unless it has holes on the slice; show_static draws the
    holes on the slice and their colour bar."""
    viewer = hole_viewer(two_holes)
    bar = viewer.plotter.scalar_bars['Au']
    viewer.set_slicer(enabled=True, origin=(0, 0, 15), normal=(0, 0, 1), tolerance=5)
    assert not bar.GetVisibility()
    fig, ax = viewer.show_static()
    without, bars = np.asarray(ax.images[0].get_array()), len(fig.axes)
    viewer.set_slicer(show_holes=True)
    assert bar.GetVisibility()
    fig, ax = viewer.show_static()
    assert len(fig.axes) == bars + 1 and not np.array_equal(without, np.asarray(ax.images[0].get_array()))


def test_view_round_trip_with_the_slicer(two_holes, grid3d):
    """get_view holds the slicer and set_view gives it back, with the holes hidden as saved."""
    import json
    data, griddef = grid3d
    viewer = hole_viewer(two_holes)
    gs.view3d(data, 'value', griddef=griddef, viewer=viewer)
    viewer.set_slicer(enabled=True, origin=(10, 0, 10.25), normal=(0, 1, 1), tolerance=7, show_holes=True)
    view = json.loads(json.dumps(viewer.get_view()))
    fresh = hole_viewer(two_holes)
    gs.view3d(data, 'value', griddef=griddef, viewer=fresh)
    fresh.set_view(view)
    assert fresh.get_view() == view
    assert not fresh._layers[0].visible and fresh._layers[0].crossing is not None
    assert fresh._layers[1].mesh.n_cells == viewer._layers[1].mesh.n_cells < 40 ** 3


def test_slicer_moves_do_not_rebuild_the_legend_or_remove_the_slice(grid3d, monkeypatch):
    """While the plane moves, the layer is replaced in place (same name, nothing removed first) and its
    colour bar is the same object with the same limits, so the slice and the legend do not blink."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, clim=(-1, 1), viewer=gs.Viewer(off_screen=True))
    viewer.set_slicer(enabled=True, origin=(40, 40, 10.25))
    plotter = viewer.plotter
    bar = plotter.scalar_bars['value']
    n_actors = len(plotter.renderer.actors)
    removed = []
    remove = plotter.remove_actor
    monkeypatch.setattr(plotter, 'remove_actor', lambda *a, **k: (removed.append(a), remove(*a, **k))[1])
    for z in (11.25, 12.25, 13.25, 12.25, 11.25):
        viewer.set_slicer(origin=(40, 40, z))
        assert plotter.scalar_bars['value'] is bar
        assert len(plotter.renderer.actors) == n_actors
        assert list(plotter.renderer.actors).count('layer value') == 1
        assert tuple(viewer._layers[0].actor.mapper.scalar_range) == (-1, 1)
        assert tuple(bar.GetLookupTable().GetRange()) == (-1, 1)
    assert removed == []


def test_slicer_moves_keep_the_holes_on_the_slice_in_place(two_holes):
    """The holes on the slice are also replaced in place, without leaving actors behind."""
    viewer = hole_viewer(two_holes)
    viewer.set_slicer(enabled=True, origin=(0, 0, 25), normal=(0, 0, 1), tolerance=8, show_holes=True)
    n_actors = len(viewer.plotter.renderer.actors)
    for z in (24, 22, 18, 12):
        viewer.set_slicer(origin=(0, 0, z))
        assert len(viewer.plotter.renderer.actors) == n_actors
        assert viewer._layers[0].crossing is not None
    viewer.set_slicer(show_holes=False)
    assert len(viewer.plotter.renderer.actors) == n_actors - 1


def test_colour_bars_of_different_titles_have_their_own_places(grid3d):
    """Each colour bar gets its own place, side by side, so bars of several layers do not overlap."""
    data, griddef = grid3d
    viewer = gs.view3d(data, 'value', griddef=griddef, viewer=gs.Viewer(off_screen=True))
    gs.view3d(gs.ExampleData('3d_estimate'), 'Estimate', viewer=viewer)
    gs.view3d(data, 'value', griddef=griddef, viewer=viewer)          # a second layer with the same title
    first, second = (viewer.plotter.scalar_bars[t].GetPosition() for t in ('value', 'Estimate'))
    assert first[0] != second[0] and list(viewer.plotter.scalar_bars.keys()) == ['value', 'Estimate']


def test_colour_bars_stop_at_the_room_there_is():
    """As PyVista does, more colour bars than there is room for is refused, not drawn off screen."""
    viewer = gs.Viewer(off_screen=True)
    with pytest.raises(RuntimeError, match='room for 10 colour bars'):
        for i in range(12):
            mesh = pv.Sphere(center=(i * 3, 0, 0), radius=1)
            mesh.cell_data['v{}'.format(i)] = np.arange(mesh.n_cells, dtype=float)
            gs.view3d(mesh, 'v{}'.format(i), viewer=viewer)
    assert len(viewer.plotter.scalar_bars.keys()) == 10


def test_set_background(grid3d):
    """The background and the text on it change together, also for a colour bar made afterwards,
    and a view saved before the background existed shows with the white one."""
    viewer = gs.view3d(grid3d[0], 'value', griddef=grid3d[1], viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    viewer.set_background('Black')
    assert viewer._background == 'black'
    assert pv.Color(viewer.plotter.background_color) == pv.Color('#3a3a3a')
    gs.view3d(gs.ExampleData('oilsands'), 'Bitumen', viewer=viewer)
    for bar in viewer.plotter.scalar_bars.values():
        assert bar.GetTitleTextProperty().GetColor() == (1.0, 1.0, 1.0)
    view = viewer.get_view()
    del view['background']
    viewer.set_view(view)
    assert viewer._background == 'white'
    assert pv.Color(viewer.plotter.background_color) == pv.Color('white')
    with pytest.raises(ValueError, match="'white', 'black'"):
        viewer.set_background('red')
    assert viewer._background == 'white'


def test_show_static_uses_the_background(grid3d):
    """The picture has the scene's background colour in its corner."""
    viewer = gs.view3d(grid3d[0], 'value', griddef=grid3d[1], viewer=gs.Viewer(off_screen=True))
    viewer.set_background('black')
    fig, ax = viewer.show_static(figsize=(3, 3), cbar=False)
    np.testing.assert_array_equal(ax.images[0].get_array()[0, 0][:3], [0x3a, 0x3a, 0x3a])
    viewer.set_background('white')
    fig, ax = viewer.show_static(figsize=(3, 3), cbar=False)
    np.testing.assert_array_equal(ax.images[0].get_array()[0, 0][:3], [255, 255, 255])


def test_set_scene(grid3d):
    """The box, the arrows and the corner marker are drawn when asked for and taken away when not,
    and the ruler is added by asking for it and stays hidden when switched off."""
    viewer = gs.view3d(grid3d[0], 'value', griddef=grid3d[1], viewer=gs.Viewer(off_screen=True))
    def actors():
        return viewer.plotter.renderer.actors

    assert 'outline' not in actors() and 'axis arrows' not in actors() and viewer._ruler is None
    viewer.set_scene(outline=True, axis_arrows=True, corner_axes=True, ruler=True)
    assert viewer._ruler is not None and viewer.plotter.renderer.cube_axes_actor is not None
    np.testing.assert_allclose(actors()['outline'].mapper.dataset.bounds, viewer._layers[0].mesh.bounds)
    assert 'axis arrows' in actors() and 'axis arrow labels-labels' in actors()
    assert viewer.plotter.renderer.axes_widget.GetEnabled()
    viewer.set_scene(outline=False, axis_arrows=False, corner_axes=False, ruler=False)
    assert 'outline' not in actors() and 'axis arrows' not in actors() and 'axis arrow labels-labels' not in actors()
    assert not viewer.plotter.renderer.axes_widget.GetEnabled()
    assert viewer.plotter.renderer.cube_axes_actor is None
    # Drawing again (a new layer, the exaggeration) does not bring the ruler back
    viewer.set_vertical_exaggeration(2)
    assert viewer.plotter.renderer.cube_axes_actor is None


def test_set_scene_follows_the_layers_and_the_scale(grid3d):
    """The box wraps every layer added after it was switched on; the z arrow has the scene's
    exaggeration in its length but not in its shape (its actor is not scaled again)."""
    viewer = gs.view3d(grid3d[0], 'value', griddef=grid3d[1], viewer=gs.Viewer(off_screen=True))
    viewer.set_scene(outline=True, axis_arrows=True)
    viewer.set_vertical_exaggeration(3)
    gs.view3d(gs.ExampleData('oilsands'), 'Bitumen', viewer=viewer)
    lower, upper = gs.plotting3d.viewer._corners(viewer._extents)
    np.testing.assert_allclose(viewer.plotter.renderer.actors['outline'].mapper.dataset.bounds,
                               np.stack([lower, upper], axis=1).ravel())
    assert viewer.plotter.renderer.actors['axis arrows'].GetScale() == (1.0, 1.0, 1.0)
    bounds = viewer.plotter.renderer.actors['axis arrows'].mapper.dataset.bounds
    assert bounds[5] - bounds[4] == pytest.approx(bounds[1] - bounds[0], rel=0.2)


def test_set_axis_ticks(grid3d):
    """A chosen spacing sets the labels of the ruler; one that VTK does not label is replaced by a
    multiple of it, with a note; a spacing of 0 or less is refused and changes nothing."""
    viewer = gs.view3d(grid3d[0], 'value', griddef=grid3d[1], viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    cube = viewer.plotter.renderer.cube_axes_actor
    viewer.set_axis_ticks(xy=20)
    assert cube.GetBounds()[:2] == (0.0, 80.0) and cube.GetXAxisRange() == (0.0, 80.0)
    assert viewer.tick_note == ''
    # 9 labels are asked for along the 40 high elevation, but VTK puts 5 and not 0, 5, 10 ... 40
    viewer.set_axis_ticks(z=5)
    assert viewer.tick_note == 'Elevation: 5 does not fit, using 10'
    assert viewer.plotter.renderer.cube_axes_actor.GetZAxisRange() == (0.0, 40.0)
    # Only a multiple of the spacing replaces it: 10 would fit, 20 is the multiple of 4 that does
    viewer.set_axis_ticks(z=4)
    assert viewer.tick_note == 'Elevation: 4 does not fit, using 20'
    viewer.set_axis_ticks(xy=1, z=1)
    assert viewer.tick_note.startswith('Easting: 1 does not fit, using 10; Northing:')
    viewer.set_axis_ticks()
    assert viewer.tick_note == '' and viewer._tick_spacing == (None, None)
    viewer.set_axis_ticks(xy=20)
    with pytest.raises(ValueError, match='above 0'):
        viewer.set_axis_ticks(xy=0)
    assert viewer._tick_spacing == (20, None)


def test_view_without_a_ruler_does_not_add_one(grid3d):
    """A view saved with no ruler says so, and showing it on a scene with none adds none."""
    viewer = gs.view3d(grid3d[0], 'value', griddef=grid3d[1], viewer=gs.Viewer(off_screen=True))
    view = viewer.get_view()
    assert view['scene']['ruler'] is False
    viewer.set_view(view)
    assert viewer._ruler is None


def test_set_axis_ticks_with_an_origin(grid3d):
    """The labels still land on multiples of the chosen spacing when they count from another origin
    (the box starts at -1240 in those counts, which is -5.5 in the scene)."""
    viewer = gs.view3d(grid3d[0], 'value', griddef=grid3d[1], viewer=gs.Viewer(off_screen=True))
    viewer.add_ruler()
    viewer.set_axis_origin((1234.5, 0, 0))
    viewer.set_axis_ticks(xy=20)
    cube = viewer.plotter.renderer.cube_axes_actor
    assert cube.GetXAxisRange() == (-1240.0, -1140.0) and cube.GetBounds()[0] == pytest.approx(-5.5)
