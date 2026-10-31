#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the interactive layer panel of the 3D viewer (needs PyVista and trame)"""

import re
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
    names = re.findall(r"trigger\('(trigger__\d+)', \[layer.name, [^\]]*\]\)", panel.layout.html)
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
