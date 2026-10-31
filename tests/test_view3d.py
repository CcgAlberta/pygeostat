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
