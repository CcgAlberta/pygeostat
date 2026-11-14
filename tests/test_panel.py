#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the interactive layer panel of the 3D viewer (needs PyVista and trame)"""

import re
import numpy as np
import pytest
import pygeostat as gs

pv = pytest.importorskip('pyvista')
pytest.importorskip('trame_pyvista')


@pytest.fixture(autouse=True)
def close_scenes():
    """Close the PyVista plotters each test leaves open."""
    yield
    pv.close_all()


@pytest.fixture
def viewer():
    """Two layers: a block model and a set of points."""
    data = gs.ExampleData('3d_grid')
    viewer = gs.view3d(data, 'value', griddef=gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1]),
                       viewer=gs.Viewer(off_screen=True))
    return gs.view3d(gs.ExampleData('3d_estimate'), 'Estimate', viewer=viewer)


def rows(panel):
    """(name, visible, opacity) of each layer in the page."""
    return [(l['name'], l['visible'], l['opacity']) for l in panel.state.layers]


def triggered(panel, handler):
    """The function the page calls for `handler`, found from the page's own HTML."""
    controller = panel.server.controller
    names = re.findall(r"trigger\('(trigger__\d+)', \[[^\]]*\]\)", panel.layout.html)
    return next(controller.trigger_fn(n) for n in names if controller.trigger_fn(n) == handler)


def test_panel_lists_the_layers(viewer):
    """The page starts with one entry per layer, with its visibility and opacity."""
    viewer.set_layer('Estimate', opacity=0.5)
    panel = viewer.panel()
    assert rows(panel) == [('value', True, 1.0), ('Estimate', True, 0.5)]


def test_panel_controls_change_the_viewer(viewer):
    """The checkbox and slider of the page call the handlers, which change the layer in the
    viewer and the page."""
    panel = viewer.panel()
    triggered(panel, panel.set_visible)('Estimate', False)
    triggered(panel, panel.set_opacity)('value', 0.25)
    blocks, points = [layer.actor for layer in viewer._layers]
    assert not points.visibility and blocks.prop.opacity == 0.25
    assert rows(panel) == [('value', True, 0.25), ('Estimate', False, 1.0)]


def test_panel_without_trame_says_what_to_install(viewer, monkeypatch):
    """A missing trame gives the install message, not a bare ImportError."""
    import sys
    monkeypatch.setitem(sys.modules, 'trame_pyvista.ui', None)
    with pytest.raises(ImportError, match=r'pygeostat\[view3d-panel\]'):
        viewer.panel()


def test_panel_follows_the_script(viewer):
    """Layers added or changed in a script after panel() show in the page."""
    panel = viewer.panel()
    gs.view3d(gs.ExampleData('oilsands'), 'Bitumen', viewer=viewer)
    viewer.set_layer('value', visible=False)
    assert [(l['name'], l['visible']) for l in panel.state.layers] == [
        ('value', False), ('Estimate', True), ('Bitumen', True)]


def test_panel_colour_limits(viewer):
    """The typed boxes set one end of the limits; the page shows the data range beside each, and a
    value that does not work leaves the limits and says why."""
    panel = viewer.panel()
    low, high = viewer._layers[0].mesh.get_data_range('value')
    assert panel.state.layers[0]['data_range'] == [low, high]
    set_clim = triggered(panel, panel.set_clim)
    set_clim('value', 1, '1.5')
    set_clim('value', 0, '0.5')
    assert viewer._layers[0].clim == (0.5, 1.5)
    assert panel.state.layers[0]['clim'] == [0.5, 1.5] and panel.state.message == ''
    set_clim('value', 0, '2')
    assert viewer._layers[0].clim == (0.5, 1.5)
    assert 'min below max' in panel.state.message and panel.state.layers[0]['clim'] == [0.5, 1.5]
    set_clim('value', 1, 'abc')
    assert viewer._layers[0].clim == (0.5, 1.5) and panel.state.message


def test_panel_no_colour_limits_for_categories():
    """A layer coloured by category, or not by a variable, has no limits boxes."""
    viewer = gs.view3d(gs.ExampleData('point3d_ind_mv'), 'Lithofacies', catdict={1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e'},
                       viewer=gs.Viewer(off_screen=True))
    viewer.add(pv.Sphere(), color='red')
    assert [l['clim'] for l in viewer.panel().state.layers] == [None, None]


def test_panel_data_range_is_the_whole_layer():
    """After a threshold or slice the label still shows the range of the whole variable, which
    the colours keep."""
    data = gs.ExampleData('3d_grid')
    viewer = gs.view3d(data, 'value', griddef=gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1]),
                       threshold=(0.5, 1.0), viewer=gs.Viewer(off_screen=True))
    row = viewer.panel().state.layers[0]
    assert row['data_range'] == [data.data['value'].min(), data.data['value'].max()]
    assert row['clim'] == row['data_range']


def test_panel_threshold(viewer):
    """A checkbox and a typed minimum. On shows the values from the minimum
    up, which starts at the lowest of the data; off shows everything. The minimum typed while it
    is off is kept, and one that does not work is explained."""
    panel = viewer.panel()
    layer = viewer._layers[0]
    low, high = layer.mesh.get_data_range('value')
    row = panel.state.layers[0]
    assert row['threshold_range'] == [low, high] and row['threshold_on'] is False and row['variable'] == 'value'
    triggered(panel, panel.set_threshold_min)('value', '0.5')
    assert layer.threshold is None and layer.mesh.n_cells == 40 ** 3
    triggered(panel, panel.set_threshold_on)('value', True)
    assert layer.threshold == (0.5, high)
    values = layer.mesh.cell_data['value']
    assert 0 < len(values) < 40 ** 3 and values.min() >= 0.5
    assert panel.state.layers[0]['threshold_on'] and panel.state.layers[0]['threshold_min'] == 0.5
    triggered(panel, panel.set_threshold_min)('value', '1')
    assert layer.threshold == (1.0, high) and layer.mesh.cell_data['value'].min() >= 1
    triggered(panel, panel.set_threshold_min)('value', 'abc')
    assert "'abc' is not a number" in panel.state.message and layer.threshold == (1.0, high)
    triggered(panel, panel.set_threshold_min)('value', str(high + 1))
    assert 'min not above max' in panel.state.message and layer.threshold == (1.0, high)
    triggered(panel, panel.set_threshold_on)('value', False)
    assert layer.threshold is None and layer.mesh.n_cells == 40 ** 3
    triggered(panel, panel.set_threshold_on)('value', True)
    assert layer.threshold == (1.0, high)


def test_panel_threshold_starts_at_the_lowest_value(viewer):
    """With nothing typed, switching the threshold on keeps every value (the minimum is the lowest)."""
    panel = viewer.panel()
    triggered(panel, panel.set_threshold_on)('value', True)
    assert viewer._layers[0].threshold == tuple(viewer._layers[0].mesh.get_data_range('value'))
    assert viewer._layers[0].mesh.n_cells == 40 ** 3


def test_panel_axis_origin():
    """The boxes set where the ruler counts from, one axis at a time; All axes from 0 uses the
    lowest edge of the scene; the card is only there when the viewer has a ruler."""
    viewer = gs.Viewer(off_screen=True)
    viewer.add(pv.Box(bounds=(1000, 1100, 5000, 5100, 300, 340)), color='red')
    panel = viewer.panel()
    assert panel.state.axes['has_ruler'] is False
    viewer.add_ruler()
    assert panel.state.axes == {'has_ruler': True, 'origin': [0.0, 0.0, 0.0], 'from_zero': False}

    def handler(function):
        controller = panel.server.controller
        names = re.findall(r"trigger\('(trigger__\d+)', \[[^\]]*\]\)", panel.layout.html)
        return next(controller.trigger_fn(n) for n in names if controller.trigger_fn(n) == function)

    handler(panel.set_origin)(0, '1000')
    handler(panel.set_origin)(2, '300')
    assert viewer._axis_origin == (1000.0, 0.0, 300.0)
    assert panel.state.axes['origin'] == [1000.0, 0.0, 300.0]
    handler(panel.set_origin)(1, 'abc')
    assert "'abc' is not a number" in panel.state.message and viewer._axis_origin == (1000.0, 0.0, 300.0)
    handler(panel.set_from_zero)(True)
    assert panel.state.axes == {'has_ruler': True, 'origin': [1000.0, 5000.0, 300.0], 'from_zero': True}
    handler(panel.set_from_zero)(False)
    assert panel.state.axes['origin'] == [1000.0, 0.0, 300.0]


def test_panel_saves_and_loads_a_view(viewer, tmp_path):
    """Save writes the view to a json file in the folder and lists it; Load gives it back, also on
    a fresh viewer; a name with a path in it, a missing file and a file that is not a view are
    explained and change nothing."""
    import json
    folder = str(tmp_path / 'views')
    panel = viewer.panel(views_folder=folder)
    assert panel.state.view_files == []
    viewer.set_layer('value', opacity=0.3, threshold=(0, 1))
    panel.save_view('first look')
    assert panel.state.view_files == ['first look'] and panel.state.view_file == 'first look'
    with open(folder + '/first look.json') as f:
        assert json.load(f) == viewer.get_view()
    viewer.set_layer('value', opacity=1, threshold=False)
    panel.load_view('first look')
    assert viewer._layers[0].kwargs['opacity'] == 0.3 and viewer._layers[0].threshold == (0.0, 1.0)
    assert panel.state.note == 'Loaded {}/first look.json'.format(folder)
    for bad in ('../escape', 'a/b'):
        panel.save_view(bad)
        assert 'letters, digits' in panel.state.message
    assert sorted(p.name for p in tmp_path.iterdir()) == ['views']
    panel.load_view('../views/first look')
    assert 'letters, digits' in panel.state.message
    panel.load_view('nothing')
    assert 'Could not load' in panel.state.message
    (tmp_path / 'views' / 'broken.json').write_text('{"layers": 3}')
    panel.load_view('broken')
    assert 'Could not load' in panel.state.message
    fresh = gs.view3d(gs.ExampleData('3d_grid'), 'value', griddef=gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1]),
                      viewer=gs.Viewer(off_screen=True))
    gs.view3d(gs.ExampleData('3d_estimate'), 'Estimate', viewer=fresh)
    other = fresh.panel(views_folder=folder)
    assert other.state.view_files == ['broken', 'first look']
    other.load_view('first look')
    assert fresh._layers[0].kwargs['opacity'] == 0.3


def test_panel_shows_itself_in_a_notebook_without_await(viewer, monkeypatch):
    """Displaying the panel gives a placeholder, which becomes the page (an iframe) once the server
    is up, so the user needs no `await layout.ready`."""
    import asyncio
    import IPython.display
    shown = []
    monkeypatch.setattr(IPython.display, 'display', lambda widget, *a, **k: shown.append(widget))
    panel = viewer.panel()

    async def cell():
        panel._ipython_display_()
        assert 'Starting' in shown[0].value
        await panel._starting
        assert '<iframe' in shown[0].value and 'not running' not in shown[0].value
        await panel.server.stop()
    asyncio.run(cell())


def test_panel_view_uses_the_panels_server(viewer):
    """The 3D view is drawn through the panel's own server (on another one the area stays empty),
    also for a second panel of the same viewer."""
    first = viewer.panel()
    second = viewer.panel()
    assert first._view.server is first.server
    assert second._view.server is second.server
    assert first.server is not second.server


def drag(panel, **values):
    """Move sliders of the page, as the mouse does: the page sets its values, and the server hears it."""
    import asyncio
    for key, value in values.items():
        panel.state[key] = value
    panel.state.flush()
    return asyncio.sleep(0.02)


def test_panel_slicer(viewer):
    """The slicer card of the page drives Viewer.set_slicer: it starts at the centre of the scene
    facing up, and each control changes the plane, the cut layers and the drill holes."""
    import asyncio
    import pandas as pd
    two_holes = pd.DataFrame({'x': [0.0] * 4 + [50.0] * 4, 'y': 0.0, 'z': [30.0, 20.0, 10.0, 0.0] * 2,
                              'hole': ['a'] * 4 + ['b'] * 4, 'Au': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]})
    gs.view3d(two_holes, 'Au', x='x', y='y', z='z', dh='hole', trace=True, radius=3, viewer=viewer)
    panel = viewer.panel()
    state = panel.state.slicer
    assert state['enabled'] is False and state['normal'] == [0.0, 0.0, 1.0] and state['tolerance'] == 10.0
    assert len(state['bounds']) == 6 and len(state['origin']) == 3
    blocks, holes = viewer._layers[0], viewer._layers[2]

    async def scenario():
        await panel.layout.ready
        triggered(panel, panel.set_slicer_enabled)(True)
        assert panel.state.slicer['enabled'] and not holes.visible
        assert [row['visible'] for row in panel.state.layers] == [True, True, False]
        await drag(panel, slicer_origin_2=10.25)
        np.testing.assert_allclose(blocks.mesh.points[:, 2], 10.25)
        assert panel.state.slicer['origin'][2] == 10.25
        triggered(panel, panel.snap_slicer)('xz')
        await asyncio.sleep(0.02)
        assert panel.state.slicer['normal'] == [0.0, 1.0, 0.0] and panel.state.slicer_normal_1 == 1.0
        await drag(panel, slicer_normal_0=1.0)
        assert viewer._slicer['normal'] == (1.0, 1.0, 0.0)
        await drag(panel, slicer_origin_0=0.0, slicer_origin_1=0.0, slicer_tolerance=6.0)
        triggered(panel, panel.set_show_holes)(True)
        assert panel.state.slicer['tolerance'] == 6.0 and panel.state.slicer['show_holes']
        assert holes.crossing is not None
        await drag(panel, slicer_tolerance=-5.0)
        assert '0 or more' in panel.state.message and viewer._slicer['tolerance'] == 6.0
        assert panel.state.slicer_tolerance == 6.0
        triggered(panel, panel.set_slicer_enabled)(False)
        assert holes.visible and blocks.mesh.n_cells == 40 ** 3 and holes.crossing is None
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_slicer_follows_the_slider_while_it_moves(viewer, monkeypatch):
    """As the slider moves up and down, the cut moves with every step, not only when it stops; and
    while the slicer is off, moving it draws nothing again. The sliders are plain values of the
    page (v-model on their own key)."""
    import asyncio
    panel = viewer.panel()
    html = panel.layout.html
    for key in ['slicer_origin_{}'.format(i) for i in range(3)] + ['slicer_normal_{}'.format(i) for i in range(3)] \
            + ['slicer_tolerance']:
        slider = re.search(r'<VSlider v-model="{}"[^>]*>'.format(key), html).group(0)
        assert '@' not in slider
    blocks = viewer._layers[0]
    placed = []
    place = viewer._place
    monkeypatch.setattr(viewer, '_place', lambda layer: (placed.append(layer.name), place(layer)))

    async def scenario():
        await panel.layout.ready
        await drag(panel, slicer_origin_2=7.25)
        assert placed == [] and blocks.mesh.n_cells == 40 ** 3       # slicer off: nothing is redrawn
        triggered(panel, panel.set_slicer_enabled)(True)
        for z in (2.25, 3.25, 4.25, 5.25, 4.25, 3.25, 9.25):          # up, up, up, then down, down
            await drag(panel, slicer_origin_2=z)
            np.testing.assert_allclose(blocks.mesh.points[:, 2], z)
            assert blocks.mesh.n_cells == 40 * 40
        # The keys follow the viewer when a script moves the plane
        viewer.set_slicer(origin=(1, 2, 3), normal=(0, 1, 0), tolerance=7)
        assert (panel.state.slicer_origin_0, panel.state.slicer_origin_1, panel.state.slicer_origin_2) == (1, 2, 3)
        assert panel.state.slicer_normal_1 == 1.0 and panel.state.slicer_tolerance == 7.0
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_vertical_exaggeration_follows_the_slider(viewer):
    """The slider is a plain value of the page that stretches the scene as it moves, and it
    follows the viewer when a script sets the exaggeration."""
    import asyncio
    panel = viewer.panel()
    slider = re.search(r'<VSlider v-model="vertical_exaggeration"[^>]*>', panel.layout.html).group(0)
    assert '@' not in slider and panel.state.vertical_exaggeration == 1.0

    async def scenario():
        await panel.layout.ready
        await drag(panel, vertical_exaggeration=4.0)
        assert viewer._vertical_exaggeration == 4.0
        assert viewer.plotter.scale[2] == 4.0
        viewer.set_vertical_exaggeration(2.5)
        assert panel.state.vertical_exaggeration == 2.5
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_background(viewer):
    """Choosing the background changes the scene's background and the text colour of the ruler and
    the colour bars, and the page follows a script."""
    import asyncio
    viewer.add_ruler()
    panel = viewer.panel()
    assert panel.state.background == 'White'

    async def scenario():
        await panel.layout.ready
        await drag(panel, background='Black')
        assert viewer._background == 'black'
        assert pv.Color(viewer.plotter.background_color) == pv.Color('#3a3a3a')
        assert viewer.plotter.renderer.cube_axes_actor.GetXAxesLinesProperty().GetColor() == (1.0, 1.0, 1.0)
        assert all(bar.GetTitleTextProperty().GetColor() == (1.0, 1.0, 1.0)
                   for bar in viewer.plotter.scalar_bars.values())
        viewer.set_background('white')
        assert panel.state.background == 'White'
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_scene_and_ticks(viewer):
    """The boxes of the scene and the tick selects call the viewer, the page follows a script, and a
    spacing that does not fit is explained in the page."""
    import asyncio
    viewer.add_ruler()
    panel = viewer.panel()
    state = panel.state
    assert (state.ruler_visible, state.outline_visible, state.axis_arrows, state.axes_visible) == (
        True, False, False, False)
    assert state.tick_spacing_xy == 'Auto' and state.tick_options_xy[0] == {'title': 'Auto', 'value': 'Auto'}
    assert [o['value'] for o in state.tick_options_xy[1:]] == [10.0, 20.0]

    async def scenario():
        await panel.layout.ready
        await drag(panel, outline_visible=True, axis_arrows=True, axes_visible=True)
        assert (viewer._outline, viewer._axis_arrows, viewer._corner_axes) == (True, True, True)
        await drag(panel, ruler_visible=False)
        assert viewer.plotter.renderer.cube_axes_actor is None and not viewer._ruler_visible
        await drag(panel, ruler_visible=True)
        await drag(panel, tick_spacing_xy=20)
        assert viewer._tick_spacing == (20.0, None) and viewer.plotter.renderer.cube_axes_actor is not None
        await drag(panel, tick_spacing_z=1.0)
        assert state.tick_note.startswith('Elevation: 1 does not fit')
        viewer.set_scene(outline=False)
        viewer.set_axis_ticks()
        assert (state.outline_visible, state.tick_spacing_xy, state.tick_note) == (False, 'Auto', '')
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_camera_buttons(viewer):
    """The view buttons, the projection box and Fit view call the viewer, and the page follows a
    script; each button sends the page one frame."""
    import asyncio
    panel = viewer.panel()
    assert panel.state.parallel_projection is False
    sent = []
    update = panel._view.update
    panel._view.update = lambda *args, **kwargs: (sent.append(1), update(*args, **kwargs))
    viewer.plotter._on_render_callbacks.add(lambda *args: sent.append(1))

    async def scenario():
        await panel.layout.ready
        triggered(panel, panel.show_view)('top')
        assert sent == [1]
        assert viewer.plotter.camera.parallel_projection and panel.state.parallel_projection
        assert viewer._hidden_axis == 'z'
        triggered(panel, panel.show_view)('iso')
        assert not panel.state.parallel_projection and viewer._hidden_axis is None
        await drag(panel, parallel_projection=True)
        assert viewer.plotter.camera.parallel_projection
        viewer.plotter.camera.zoom(5)
        triggered(panel, panel.fit_view)()
        assert viewer.plotter.camera.parallel_projection
        viewer.set_camera_view('x')
        assert panel.state.parallel_projection is False
        with pytest.raises(ValueError):
            panel.viewer.set_camera_view('up')
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_colormap_and_cell_edges(viewer):
    """The colormap list and the edges box call the viewer; the page lists the colormaps in groups
    with a preview bar each, adds one made by hand, and follows a script."""
    import asyncio
    from matplotlib.colors import ListedColormap
    gs.view3d(gs.ExampleData('3d_grid'), 'value', griddef=gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1]),
              cmap=ListedColormap(['red', 'blue'], name='mine'), label='mine', viewer=viewer)
    panel = viewer.panel()
    state = panel.state
    assert [(row['name'], row['cmap'], row['edges']) for row in state.layers] == [
        ('value', 'viridis', False), ('Estimate', 'viridis', None), ('mine', 'mine', False)]
    items = state.cmap_items
    groups = [item['title'] for item in items if item.get('type') == 'subheader']
    assert groups[:2] == ['In the scene', 'pygeostat'] and 'Perceptually uniform' in groups
    names = {item['value']: item['gradient'] for item in items if 'value' in item}
    assert {'mine', 'topo1', 'viridis', 'plasma'} <= set(names)
    assert names['mine'].startswith('linear-gradient(to right, #ff0000') and names['mine'].endswith('#0000ff)')
    assert panel._cmap_items() == items

    async def scenario():
        await panel.layout.ready
        triggered(panel, panel.set_cmap)('value', 'plasma')
        assert viewer._layers[0].cmap == 'plasma' and state.layers[0]['cmap'] == 'plasma'
        triggered(panel, panel.set_edges)('value', True)
        assert viewer._layers[0].actor.prop.show_edges and state.layers[0]['edges'] is True
        triggered(panel, panel.set_cmap)('Estimate', 'nonsense')
        assert "Unknown colormap 'nonsense'" in state.message and state.layers[1]['cmap'] == 'viridis'
        viewer.set_layer('value', cmap='magma', edges=False)
        assert (state.layers[0]['cmap'], state.layers[0]['edges']) == ('magma', False)
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_tube_radius(viewer):
    """The radius slider of a drill hole layer draws the tubes again when it is let go, from a quarter
    to four times the radius the layer had; the page follows a script, even past that range."""
    import asyncio
    viewer.set_layer('Estimate', visible=False)
    gs.view3d(gs.ExampleData('oilsands'), 'Bitumen', trace=True, radius=8, viewer=viewer)
    panel = viewer.panel()
    state = panel.state
    rows = {row['name']: row for row in state.layers}
    assert rows['value']['radius'] is None and rows['Bitumen']['radius'] == 8.0
    assert rows['Bitumen']['radius_range'] == [2.0, 32.0, 2.0] and state.radius_ui == {'Bitumen': 8.0}
    slider = re.search(r'<VSlider[^>]*Tube radius[^>]*>', panel.layout.html).group(0)
    assert 'v-model="radius_ui[layer.name]"' in slider and '@end=' in slider and '@keyup=' in slider

    async def scenario():
        await panel.layout.ready
        triggered(panel, panel.set_radius)('Bitumen', 12)
        assert viewer._layers[2].radius == 12.0 and state.radius_ui == {'Bitumen': 12.0}
        triggered(panel, panel.set_radius)('Bitumen', 'wide')
        assert "'wide' is not a number" in state.message and viewer._layers[2].radius == 12.0
        viewer.set_layer('Bitumen', radius=100)
        assert state.layers[2]['radius_range'] == [2.0, 100.0, 2.0] and state.radius_ui == {'Bitumen': 100.0}
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_colour_by_and_palette():
    """A layer with other variables has a Colour by list; a layer coloured by categories has a palette
    list with only the palettes that have enough colours; both call the viewer."""
    import asyncio
    import pandas as pd
    two_holes = pd.DataFrame({'x': [0.0] * 4 + [50.0] * 4, 'y': 0.0, 'z': [30.0, 20.0, 10.0, 0.0] * 2,
                              'hole': ['a'] * 4 + ['b'] * 4, 'Au': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]})
    two_holes['Cu'] = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    two_holes['code'] = [1, 1, 3, 3] * 2
    viewer = gs.view3d(two_holes, 'Au', x='x', y='y', z='z', dh='hole', trace=True, radius=3,
                       other_vars={'Cu': None, 'code': {1: 'a', 3: 'b'}}, viewer=gs.Viewer(off_screen=True))
    gs.view3d(two_holes, 'Au', x='x', y='y', z='z', dh='hole', trace=True, label='alone', other_vars=[], viewer=viewer)
    panel = viewer.panel()
    state = panel.state
    first, alone = state.layers
    assert first['color_by'] == 'Au' and first['cmap_kind'] == 'continuous' and first['palette_items'] == []
    assert [o['title'] for o in first['color_options']] == ['Au', 'Cu', 'code (categories)']
    assert len(alone['color_options']) == 1

    async def scenario():
        await panel.layout.ready
        triggered(panel, panel.set_color_by)('Au', 'code')
        row = state.layers[0]
        assert row['color_by'] == 'code' and row['cmap_kind'] == 'categorical' and row['cmap'] is None
        palettes = [i['value'] for i in row['palette_items'] if 'value' in i]
        assert 'cat_vibrant' in palettes and 'tab10' in palettes and 'turbo' in palettes
        triggered(panel, panel.set_cmap)('Au', 'cat_pastel')
        assert state.layers[0]['cmap'] == 'cat_pastel' and viewer._layers[0].cmap == 'cat_pastel'
        triggered(panel, panel.set_color_by)('Au', 'Zn')
        assert "can be coloured by ['Au', 'Cu', 'code']" in state.message and state.layers[0]['color_by'] == 'code'
        viewer.set_layer('Au', color_by='Cu')
        assert state.layers[0]['color_by'] == 'Cu' and state.layers[0]['cmap_kind'] == 'continuous'
        await panel.server.stop()
    asyncio.run(scenario())


def test_panel_palettes_leave_out_short_ones():
    """With 13 categories, the pygeostat palettes (12, 12 and 6 colours) are not offered."""
    codes = list(range(1, 14))
    import pandas as pd
    points = pd.DataFrame({'x': np.arange(13.0), 'y': 0.0, 'z': 0.0, 'code': codes})
    viewer = gs.view3d(points, 'code', x='x', y='y', z='z', catdict={c: str(c) for c in codes},
                       viewer=gs.Viewer(off_screen=True))
    panel = viewer.panel()
    palettes = [i['value'] for i in panel.state.layers[0]['palette_items'] if 'value' in i]
    assert not {'cat_vibrant', 'cat_pastel', 'cat_dark'} & set(palettes) and 'tab20' in palettes


def test_panel_with_a_palette_made_in_view3d():
    """A layer made with a palette of categories does not put it in the list of colormaps."""
    import pandas as pd
    points = pd.DataFrame({'x': [0.0, 1.0], 'y': 0.0, 'z': 0.0, 'code': [1, 3]})
    viewer = gs.view3d(points, 'code', x='x', y='y', z='z', catdict={1: 'a', 3: 'b'}, cmap='cat_pastel',
                       viewer=gs.Viewer(off_screen=True))
    panel = viewer.panel()
    assert panel.state.layers[0]['cmap'] == 'cat_pastel'
    assert 'In the scene' not in [i['title'] for i in panel.state.cmap_items if i.get('type') == 'subheader']


def test_panel_isolates_categories():
    """The buttons of the categories show one, add others and show all again, with a colour each; the
    page follows a script, and a layer without categories has no buttons."""
    import asyncio
    import pandas as pd
    points = pd.DataFrame({'x': np.arange(3.0), 'y': 0.0, 'z': 0.0, 'code': [1, 2, 3], 'Au': [1.0, 2.0, 3.0]})
    viewer = gs.view3d(points, 'code', x='x', y='y', z='z', catdict={1: 'a', 2: 'b', 3: 'c'},
                       other_vars=['Au'], viewer=gs.Viewer(off_screen=True))
    panel = viewer.panel()
    state = panel.state
    row = state.layers[0]
    assert row['category_names'] == ['a', 'b', 'c'] and len(row['category_colours']) == 3
    assert row['category_colours'][0].startswith('#') and row['categories'] == []

    async def scenario():
        await panel.layout.ready
        shown = lambda: sorted(set(viewer._layers[0].mesh.cell_data['code']))
        triggered(panel, panel.toggle_category)('code', 'b')
        assert shown() == [2.0] and state.layers[0]['categories'] == ['b']
        triggered(panel, panel.toggle_category)('code', 'c')
        assert shown() == [2.0, 3.0] and state.layers[0]['categories'] == ['b', 'c']
        triggered(panel, panel.toggle_category)('code', 'b')
        assert shown() == [3.0]
        triggered(panel, panel.toggle_category)('code', 'a')
        triggered(panel, panel.toggle_category)('code', 'b')
        assert shown() == [1.0, 2.0, 3.0] and state.layers[0]['categories'] == []
        triggered(panel, panel.toggle_category)('code', 'a')
        triggered(panel, panel.show_all_categories)('code')
        assert shown() == [1.0, 2.0, 3.0] and state.layers[0]['categories'] == []
        viewer.set_layer('code', categories=['c'])
        assert state.layers[0]['categories'] == ['c']
        viewer.set_layer('code', color_by='Au')
        row = state.layers[0]
        assert row['category_names'] == [] and row['categories'] == [] and row['category_colours'] == []
        await panel.server.stop()
    asyncio.run(scenario())


def hole_scene():
    """Two vertical holes inside a block model that stands in front of them, with a panel."""
    import pandas as pd
    holes = pd.DataFrame({'x': [10.0] * 4 + [30.0] * 4, 'y': 10.0, 'z': [30.0, 20.0, 10.0, 0.0] * 2,
                          'hole': ['a'] * 4 + ['b'] * 4, 'Au': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]})
    viewer = gs.view3d(gs.ExampleData('3d_grid'), 'value', griddef=gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1]),
                       opacity=0.2, viewer=gs.Viewer(off_screen=True))
    gs.view3d(holes, 'Au', x='x', y='y', z='z', dh='hole', trace=True, radius=3, viewer=viewer)
    viewer.set_camera_view('y')
    return viewer, viewer.panel()


def click(viewer, point):
    """What a click on a world point picks, from pyvista's own picker."""
    import vtk
    viewer.plotter.render()
    renderer = viewer.plotter.renderer
    renderer.SetWorldPoint(*point, 1.0)
    renderer.WorldToDisplay()
    x, y = renderer.GetDisplayPoint()[:2]
    picker = vtk.vtkCellPicker()
    picker.SetTolerance(0.05)
    picker.Pick(x, y, 0, renderer)
    return picker


def test_panel_picks_a_drill_hole_through_a_block():
    """A click on a drill hole says which interval it is and marks it, though a block model is in
    front of it; a click beside it, or on a hidden layer, says nothing was picked; anything drawn
    again takes the text and the mark away."""
    viewer, panel = hole_scene()
    state = panel.state
    assert state.pick_info == 'Click a drill hole to see its values (drag to rotate).'
    panel._on_pick(None, click(viewer, (10.0, 10.0, 25.0)))
    assert state.pick_info == 'Au: hole a, elevation 30.0 to 20.0, Au = 1.00'
    mark = viewer.plotter.renderer.actors['pick highlight']
    assert mark.prop.color == (1.0, 0.0, 1.0) and mark.prop.opacity == 0.6 and not mark.pickable
    assert mark.mapper.dataset.bounds[0] == pytest.approx(10 - 5, abs=0.3)
    panel._on_pick(None, click(viewer, (30.0, 10.0, 5.0)))
    assert state.pick_info == 'Au: hole b, elevation 10.0 to 0.0, Au = 7.00'
    panel._on_pick(None, click(viewer, (20.0, 10.0, 25.0)))
    assert 'Nothing picked there' in state.pick_info and 'pick highlight' not in viewer.plotter.renderer.actors
    panel._on_pick(None, click(viewer, (10.0, 10.0, 25.0)))
    viewer.set_layer('Au', opacity=0.5)
    assert state.pick_info.startswith('Click a drill hole') and 'pick highlight' not in viewer.plotter.renderer.actors
    viewer.set_layer('Au', visible=False)
    panel._on_pick(None, click(viewer, (10.0, 10.0, 25.0)))
    assert 'Nothing picked there' in state.pick_info


def test_panel_picks_the_holes_on_the_slice():
    """A hole on the slicer's plane is picked with the note that it is on the slice."""
    viewer, panel = hole_scene()
    viewer.set_slicer(enabled=True, origin=(10, 10, 0), normal=(0, 1, 0), tolerance=20, show_holes=True)
    panel._on_pick(None, click(viewer, (10.0, 10.0, 25.0)))
    assert panel.state.pick_info == 'Au (on the slice): hole a, elevation 30.0 to 20.0, Au = 1.00'


def test_panel_pick_needs_the_click_to_be_on_the_hole():
    """A click that the picker snaps to a hole 60 pixels away picks nothing."""
    viewer, panel = hole_scene()
    picker = click(viewer, (10.0, 10.0, 25.0))
    renderer = viewer.plotter.renderer
    renderer.SetWorldPoint(10.0, 10.0, 25.0, 1.0)
    renderer.WorldToDisplay()
    x, y = renderer.GetDisplayPoint()[:2]
    assert panel._click_distance(picker) == 0.0
    import vtk
    far = vtk.vtkCellPicker()
    far.SetTolerance(0.05)
    far.Pick(x + 60, y, 0, renderer)
    # The picker's own tolerance reaches a hole 60 pixels from the click, as it did in the notebook
    assert far.GetActor() is not None and far.GetCellId() >= 0
    assert panel._click_distance(far) > 50
    panel._on_pick(None, far)
    assert 'Nothing picked there' in panel.state.pick_info


def test_panel_a_click_picks_and_a_drag_does_not():
    """The mouse button released within 4 pixels of where it was pressed picks, a longer move is a
    rotation and does not; a second panel of the viewer replaces the first one's picking."""
    viewer, panel = hole_scene()
    picks = []

    class Interactor:
        def __init__(self):
            self.position = (0, 0)

        def GetEventPosition(self):
            return self.position

        def GetPicker(self):
            return type('Picker', (), {'Pick': lambda self, x, y, z, renderer: picks.append((x, y))})()

    iren, observers = Interactor(), {}
    real = viewer.plotter.iren

    class Recorder:
        def add_observer(self, event, call):
            observers[event] = call
            return len(observers)

        def remove_observer(self, tag):
            pass

        def reset_picker(self):
            pass

        def __getattr__(self, name):
            return getattr(real, name)
    viewer.plotter.iren = Recorder()
    panel._enable_picking()
    viewer.plotter.iren = real
    press, release = observers['LeftButtonPressEvent'], observers['EndInteractionEvent']
    iren.position = (100, 100)
    press(iren, None)
    iren.position = (103, 102)
    release(iren, None)
    assert picks == [(103, 102)]
    iren.position = (100, 100)
    press(iren, None)
    iren.position = (120, 100)
    release(iren, None)
    release(iren, None)             # a release with no press before it
    assert picks == [(103, 102)]
    second = viewer.panel()
    assert second is not panel and len(viewer._pick_observers) == 2


def test_panel_picks_with_a_vertical_exaggeration():
    """With the scene stretched, a click on a drill hole still picks it (the distance on the screen
    is measured on the stretched tube)."""
    viewer, panel = hole_scene()
    viewer.set_vertical_exaggeration(3)
    panel._on_pick(None, click(viewer, (10.0, 10.0, 25.0 * 3)))
    assert panel.state.pick_info == 'Au: hole a, elevation 30.0 to 20.0, Au = 1.00'
