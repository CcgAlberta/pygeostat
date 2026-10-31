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


def triggered(panel, handler):
    """The function the page calls for `handler`, found from the page's own HTML."""
    controller = panel.server.controller
    names = re.findall(r"trigger\('(trigger__\d+)', \[layer.name, \$event\]\)", panel.layout.html)
    return next(controller.trigger_fn(n) for n in names if controller.trigger_fn(n) == handler)


def test_panel_lists_the_layers(viewer):
    """The page starts with one entry per layer, with its visibility and opacity."""
    viewer.set_layer('Estimate', opacity=0.5)
    panel = viewer.panel()
    assert panel.state.layers == [{'name': 'value', 'visible': True, 'opacity': 1.0},
                                  {'name': 'Estimate', 'visible': True, 'opacity': 0.5}]


def test_panel_controls_change_the_viewer(viewer):
    """The checkbox and slider of the page call the handlers, which change the layer in the
    viewer and the page."""
    panel = viewer.panel()
    triggered(panel, panel.set_visible)('Estimate', False)
    triggered(panel, panel.set_opacity)('value', 0.25)
    blocks, points = [layer.actor for layer in viewer._layers]
    assert not points.visibility and blocks.prop.opacity == 0.25
    assert panel.state.layers == [{'name': 'value', 'visible': True, 'opacity': 0.25},
                                  {'name': 'Estimate', 'visible': False, 'opacity': 1.0}]


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
