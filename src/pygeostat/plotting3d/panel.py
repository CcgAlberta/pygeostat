#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""panel is the interactive layer panel of a Viewer, a web page (trame) with the controls beside the 3D view"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from contextlib import contextmanager

import matplotlib as mpl

from .viewer import _corners, _gradient, _MAX_LABELS


# The colormaps offered for a layer coloured by a continuous variable: pygeostat's own two (its
# hillshading maps) and a curated set of matplotlib's. The diverging ones are reversed so that high
# values are red
_CMAP_GROUPS = {
	'pygeostat': ['topo1', 'topo2'],
	'Perceptually uniform': ['viridis', 'plasma', 'inferno', 'magma', 'cividis'],
	'Rainbow': ['turbo', 'jet', 'rainbow'],
	'Diverging': ['coolwarm', 'RdYlBu_r', 'Spectral_r', 'RdBu_r', 'seismic'],
	'Sequential': ['YlOrRd', 'hot', 'Reds', 'Blues', 'Greys'],
	'Earth': ['terrain', 'gist_earth'],
}
# The palettes offered for a layer coloured by categories: pygeostat's own, matplotlib's qualitative
# maps, and continuous maps cut into evenly spaced colours (the only ones with no limit on the
# number of categories)
_PALETTE_GROUPS = {
	'pygeostat': ['cat_vibrant', 'cat_pastel', 'cat_dark'],
	'matplotlib': ['tab10', 'tab20', 'tab20b', 'tab20c', 'Set1', 'Set2', 'Set3', 'Dark2', 'Paired', 'Accent',
				   'Pastel1'],
	'From a continuous map': ['turbo', 'viridis', 'Spectral', 'rainbow'],
}
# How far (screen pixels) from a drill hole a click may land and still pick it, and how far the
# mouse may move between pressing and releasing the button for it to be a click and not a drag
_PICK_SLACK = 8
_CLICK_SLOP = 4
_PICK_DEFAULT = 'Click a drill hole to see its values (drag to rotate).'
_PICK_MISSED = ('Nothing picked there: click right on a drill hole. Zooming in or a bigger tube '
				'radius makes them easier to hit.')
# How a colormap's preview bar looks beside its name in the list and beside the chosen one
_PREVIEW_STYLE = ("{ background: item.raw.gradient, width: '64px', height: '10px', "
				  "borderRadius: '2px', flexShrink: 0 }",)


class Panel:
	"""
	The controls of a :class:`~pygeostat.plotting3d.viewer.Viewer`, shown beside its 3D view.
	Made by :meth:`Viewer.panel() <pygeostat.plotting3d.viewer.Viewer.panel>`.

	Every control calls a method of the viewer (e.g. ``set_layer``), so what the panel does can
	also be done in a script.

	Parameters:
		viewer (gs.Viewer): the scene to control
		views_folder (str): folder of the saved views, see :meth:`Viewer.panel() <pygeostat.plotting3d.viewer.Viewer.panel>`

	Attributes:
		layout: the page. A panel shows itself when it is the last line of a notebook cell
	"""

	def __init__(self, viewer, views_folder='views'):
		import uuid
		try:
			from trame.app import get_server
			from trame.ui.vuetify3 import SinglePageWithDrawerLayout
			from trame.widgets import html, vuetify3
			from trame_pyvista.ui import get_viewer
		except ImportError:
			raise ImportError('The layer panel requires trame: pip install "pygeostat[view3d-panel]"') from None
		self.viewer = viewer
		self.views_folder = views_folder
		# The range of the tube radius slider of each drill hole layer, from the radius it first had
		self._radius_ranges = {}
		# The preview bar of each colormap offered, made once
		self._gradients = {}
		# Set while a button changes the scene in several steps, so the page is sent only the result
		self._paused = False
		# The lowest value typed for the threshold of each layer, kept while the threshold is off
		self._typed_min = {}
		# A new name each time: get_server caches by name, so a fixed one would share state
		self.server = get_server('pygeostat_{}'.format(uuid.uuid4().hex[:8]), client_type='vue3')
		self.state = self.server.state
		self.state.message = ''
		self.state.note = ''
		self.state.pick_info = _PICK_DEFAULT
		self.state.view_name = ''
		self.state.view_file = None
		self._list_views()
		# trame_pyvista keeps one viewer per plotter, tied to the server it was made with, so one
		# left from an earlier panel of this viewer would draw on that panel's server
		import trame_pyvista.ui
		trame_pyvista.ui._VIEWERS.pop(viewer.plotter._id_name, None)
		self._view = get_viewer(viewer.plotter, server=self.server)
		self.refresh()
		self._enable_picking()
		self.state.change(*['slicer_origin_{}'.format(i) for i in range(3)],
						  *['slicer_normal_{}'.format(i) for i in range(3)], 'slicer_tolerance')(
			self._on_slicer_sliders)
		self.state.change('vertical_exaggeration')(self._on_vertical_exaggeration)
		self.state.change('background')(self._on_background)
		self.state.change('outline_visible', 'ruler_visible', 'axis_arrows', 'axes_visible')(self._on_scene)
		self.state.change('tick_spacing_xy', 'tick_spacing_z')(self._on_tick_spacing)
		self.state.change('parallel_projection')(self._on_parallel_projection)

		with SinglePageWithDrawerLayout(self.server) as self.layout:
			self.layout.title.set_text('Layers')
			with self.layout.drawer:
				html.Div('{{ message }}', v_if='message', classes='text-caption text-error ma-2')
				with vuetify3.VCard(v_if='axes.has_ruler', variant='outlined', classes='ma-2'):
					with vuetify3.VCardText():
						html.Div('Axis labels count from', classes='text-caption')
						with html.Div(classes='d-flex ga-2'):
							for axis, label in enumerate(('East', 'North', 'Elev')):
								vuetify3.VTextField(
									label=label, type='number', model_value=('axes.origin[{}]'.format(axis),),
									disabled=('axes.from_zero',), density='compact', hide_details=True,
									variant='outlined',
									change=(self.set_origin, '[{}, $event.target.value]'.format(axis)))
						vuetify3.VCheckbox(label='All axes from 0', model_value=('axes.from_zero',),
										   update_modelValue=(self.set_from_zero, '[$event]'),
										   density='compact', hide_details=True)
				with vuetify3.VCard(v_for='layer in layers', key='layer.name', variant='outlined', classes='ma-2'):
					with vuetify3.VCardText():
						vuetify3.VCheckbox(label=('layer.name',), model_value=('layer.visible',),
										   update_modelValue=(self.set_visible, '[layer.name, $event]'),
										   density='compact', hide_details=True)
						vuetify3.VSlider(label='Opacity', model_value=('layer.opacity',), min=0, max=1, step=0.05,
										 density='compact', hide_details=True, thumb_label=True,
										 update_modelValue=(self.set_opacity, '[layer.name, $event]'))
						# Values tied to the data are typed, with the data's own range in the label
						vuetify3.VSelect(v_if='layer.color_options.length > 1', label='Colour by',
										 model_value=('layer.color_by',), items=('layer.color_options',),
										 density='compact', hide_details=True, classes='mt-2',
										 update_modelValue=(self.set_color_by, '[layer.name, $event]'))
						for kind, label, items in (('continuous', 'Colormap', 'cmap_items'),
												   ('categorical', 'Palette', 'layer.palette_items')):
							with html.Div(v_if="layer.cmap_kind === '{}'".format(kind), classes='mt-2'):
								with vuetify3.VSelect(
										label=label, model_value=('layer.cmap',), items=(items,),
										density='compact', hide_details=True,
										update_modelValue=(self.set_cmap, '[layer.name, $event]')):
									with vuetify3.Template(v_slot_item='{ props, item }'):
										with vuetify3.VListItem(v_bind='props', density='compact'):
											with vuetify3.Template(v_slot_append=True):
												html.Div(style=_PREVIEW_STYLE, classes='ml-4')
									with vuetify3.Template(v_slot_selection='{ item }'):
										html.Span('{{ item.title }}', classes='mr-3')
										html.Div(style=_PREVIEW_STYLE)
						vuetify3.VCheckbox(v_if='layer.edges !== null', label='Show cell edges',
										   model_value=('layer.edges',), density='compact', hide_details=True,
										   update_modelValue=(self.set_edges, '[layer.name, $event]'))
						# A button for each category in its colour: filled when shown, outlined when hidden
						with html.Div(v_if='layer.category_names.length', classes='mt-2'):
							html.Div('Isolate {{ layer.color_by }} (click to toggle)', classes='text-caption mb-1')
							with html.Div(classes='d-flex flex-wrap ga-1'):
								vuetify3.VBtn(
									'{{ category }}', v_for='(category, i) in layer.category_names', key='category',
									size='small', click=(self.toggle_category, '[layer.name, category]'),
									color=('layer.category_colours[i]',),
									variant=("layer.categories.length === 0 || layer.categories.includes(category)"
											 " ? 'flat' : 'outlined'",),
									prepend_icon=("layer.categories.includes(category) ? 'mdi-check' : undefined",))
							vuetify3.VBtn('Show all', size='small', variant='text', prepend_icon='mdi-eye',
										  classes='mt-1', disabled=('layer.categories.length === 0',),
										  click=(self.show_all_categories, '[layer.name]'))
						# The tubes are drawn again when the slider is let go, not at every step of its move
						vuetify3.VSlider(v_if='layer.radius !== null', label='Tube radius', v_model=('radius_ui[layer.name]',),
										 min=('layer.radius_range[0]',), max=('layer.radius_range[1]',),
										 step=('layer.radius_range[2]',), density='compact', hide_details=True,
										 thumb_label=True, classes='mt-2',
										 end=(self.set_radius, '[layer.name, $event]'),
										 keyup=(self.set_radius, '[layer.name, radius_ui[layer.name]]'))
						with html.Div(v_if='layer.clim', classes='mt-2'):
							html.Div('Colour limits', classes='text-caption')
							with html.Div(classes='d-flex ga-2'):
								self._range_boxes('layer.clim', 'layer.data_range', self.set_clim)
						# A checkbox, and the lowest value to show typed in a box
						with html.Div(v_if='layer.threshold_range', classes='mt-2'):
							vuetify3.VCheckbox(label='Threshold', model_value=('layer.threshold_on',),
											   update_modelValue=(self.set_threshold_on, '[layer.name, $event]'),
											   density='compact', hide_details=True)
							vuetify3.VTextField(
								label=("'Min ' + layer.variable + ' (data ' + Number(layer.threshold_range[0]).toPrecision(3)"
									   " + ' to ' + Number(layer.threshold_range[1]).toPrecision(3) + ')'",),
								type='number', model_value=('layer.threshold_min',), hint='Press Enter to apply',
								persistent_hint=True, density='compact', variant='outlined', classes='mt-1',
								change=(self.set_threshold_min, '[layer.name, $event.target.value]'))
				# The slicer is one plane for the whole scene
				with vuetify3.VCard(variant='outlined', classes='ma-2'):
					with vuetify3.VCardText():
						html.Div('Slicer', classes='text-caption')
						vuetify3.VCheckbox(label='Enable slicer', model_value=('slicer.enabled',),
										   update_modelValue=(self.set_slicer_enabled, '[$event]'),
										   density='compact', hide_details=True)
						for axis, name in enumerate('XYZ'):
							self._slider('Origin ' + name, 'slicer_origin_{}'.format(axis),
										 'slicer.bounds[{}]'.format(2 * axis), 'slicer.bounds[{}]'.format(2 * axis + 1))
						for axis, name in enumerate('XYZ'):
							self._slider('Normal ' + name, 'slicer_normal_{}'.format(axis), -1, 1, step=0.05)
						with html.Div(classes='d-flex ga-1 mb-1'):
							for orient in ('xy', 'xz', 'yz'):
								vuetify3.VBtn(orient.upper(), size='x-small', variant='tonal',
											  click=(self.snap_slicer, "['{}']".format(orient)))
						self._slider('Tolerance', 'slicer_tolerance', 0, 50)
						vuetify3.VCheckbox(label='Show drill holes on slice', model_value=('slicer.show_holes',),
										   update_modelValue=(self.set_show_holes, '[$event]'),
										   density='compact', hide_details=True)
				with vuetify3.VCard(variant='outlined', classes='ma-2'):
					with vuetify3.VCardText():
						html.Div('Scene', classes='text-caption')
						self._slider('Vertical exaggeration', 'vertical_exaggeration', 1, 10, step=0.5)
						vuetify3.VSelect(v_model=('background',), items=(['White', 'Black'],), label='Background',
										 density='compact', hide_details=True, classes='mt-2')
						for label, key in (('Outline', 'outline_visible'), ('Axis', 'ruler_visible'),
										   ('X/Y/Z arrows in the scene', 'axis_arrows'),
										   ('Coordinate legend (corner)', 'axes_visible')):
							vuetify3.VCheckbox(label=label, v_model=(key,), density='compact', hide_details=True)
						html.Div('Axis ticks (same settings = comparable pictures)', classes='text-caption mt-2')
						with html.Div(classes='d-flex ga-2'):
							for label, key in (('East/North', 'tick_spacing_xy'), ('Elevation', 'tick_spacing_z')):
								vuetify3.VSelect(label=label, v_model=(key,), items=(key.replace('spacing', 'options'),),
												 density='compact', hide_details=True)
						html.Div('{{ tick_note }}', v_if='tick_note', classes='text-caption text-warning')
				with vuetify3.VCard(variant='outlined', classes='ma-2'):
					with vuetify3.VCardText():
						html.Div('Picking', classes='text-caption')
						html.Div('{{ pick_info }}', classes='text-body-2')
				with vuetify3.VCard(variant='outlined', classes='ma-2'):
					with vuetify3.VCardText():
						html.Div('Camera', classes='text-caption')
						html.Div('2D views (flat, no perspective)', classes='text-caption')
						with html.Div(classes='d-flex flex-wrap ga-1 mb-2'):
							for label, view in (('Top 2D', 'top'), ('East-Elev 2D', 'east'), ('North-Elev 2D', 'north')):
								vuetify3.VBtn(label, size='small', variant='tonal',
											  click=(self.show_view, "['{}']".format(view)))
						html.Div('3D views (perspective)', classes='text-caption')
						with html.Div(classes='d-flex flex-wrap ga-1 mb-2'):
							for label, view in (('Iso 3D', 'iso'), ('X 3D', 'x'), ('Y 3D', 'y'), ('Z 3D', 'z')):
								vuetify3.VBtn(label, size='small', variant='tonal',
											  click=(self.show_view, "['{}']".format(view)))
						vuetify3.VCheckbox(label='Parallel projection', v_model=('parallel_projection',),
										   density='compact', hide_details=True)
						vuetify3.VBtn('Fit view', size='small', variant='tonal', classes='mt-1',
									  click=(self.fit_view, '[]'))
				with vuetify3.VCard(variant='outlined', classes='ma-2'):
					with vuetify3.VCardText():
						html.Div('Export', classes='text-caption')
						with html.Div(classes='d-flex flex-wrap ga-1'):
							for label, icon, handler, filename, mime in (
									('Screenshot (PNG)', 'mdi-file-png-box', self.attach_screenshot, 'scene.png', 'image/png'),
									('Interactive HTML', 'mdi-download', self.attach_export, 'scene.html', 'text/html')):
								vuetify3.VBtn(label, size='small', variant='tonal', prepend_icon=icon,
											  click=("utils.download('{}', trigger('{}'), '{}')".format(
												  filename, self.server.trigger_name(handler), mime)))
				with vuetify3.VCard(variant='outlined', classes='ma-2'):
					with vuetify3.VCardText():
						html.Div('Saved views', classes='text-caption')
						with html.Div(classes='d-flex align-center ga-1'):
							vuetify3.VTextField(v_model=('view_name',), label='Name (optional)', density='compact',
												hide_details=True)
							vuetify3.VBtn('Save', size='small', variant='tonal',
										  click=(self.save_view, '[view_name]'))
						with html.Div(classes='d-flex align-center ga-1 mt-2'):
							vuetify3.VSelect(v_model=('view_file',), items=('view_files',), label='Saved views',
											 density='compact', hide_details=True,
											 no_data_text='None saved yet')
							vuetify3.VBtn('Load', size='small', variant='tonal', disabled=('!view_file',),
										  click=(self.load_view, '[view_file]'))
						html.Div('{{ note }}', v_if='note', classes='text-caption mt-1')
			with self.layout.content:
				with vuetify3.VContainer(fluid=True, classes='pa-0 fill-height'):
					self._view.ui(add_menu=False)

	@staticmethod
	def _slider(label, key, low, high, step=None):
		"""A slider for a number of the slicer: its value is a plain value of
		the page (``key``) that the server follows on every move, so the cut changes while the
		slider is dragged (see ``_on_slicer_sliders``)"""
		from trame.widgets import vuetify3
		options = {} if step is None else {'step': step}
		vuetify3.VSlider(label=label, v_model=(key,), min=(str(low),), max=(str(high),),
						 density='compact', hide_details=True, thumb_label=True, **options)

	@staticmethod
	def _range_boxes(values, data_range, handler):
		"""Two typed boxes, Min and Max, for ``values`` (a pair in the page), each with the matching
		end of ``data_range`` in the label. ``handler(name, end, value)`` is called with end 0 or 1"""
		from trame.widgets import vuetify3
		for end, label in enumerate(('Min', 'Max')):
			vuetify3.VTextField(
				label=("'{} (data ' + Number({}[{}]).toPrecision(3) + ')'".format(label, data_range, end),),
				type='number', model_value=('{}[{}]'.format(values, end),), density='compact',
				hide_details=True, variant='outlined',
				change=(handler, "[layer.name, {}, $event.target.value]".format(end)))

	def _ipython_display_(self):
		"""Show the panel in a notebook. trame's own display needs ``await layout.ready`` first,
		because the server only starts when the cell is over, so a placeholder is shown that is
		replaced by the page as soon as the server is up"""
		import asyncio
		from IPython.display import display
		from ipywidgets import HTML
		widget = HTML('Starting the panel ...')
		display(widget)

		async def show_page():
			await self.layout.ready
			widget.value = self.layout._jupyter_content()
		self._starting = asyncio.ensure_future(show_page())

	def refresh(self):
		"""Show the viewer's layers as they are now. The viewer calls this when a layer is added or
		changed, so the page follows a script."""
		self._clear_pick()
		self.state.layers = [self._row(layer) for layer in self.viewer._layers]
		self.state.radius_ui = {layer.name: layer.radius for layer in self.viewer._layers if layer.radius is not None}
		slicer = self.viewer._slicer
		lower, upper = _corners(self.viewer._extents) if self.viewer._extents else ([0.0] * 3, [1.0] * 3)
		self.state.slicer = {
			'enabled': slicer['enabled'], 'origin': _floats(self._slicer_origin()) if self.viewer._extents else [0.0] * 3,
			'normal': _floats(slicer['normal']), 'tolerance': slicer['tolerance'],
			'show_holes': slicer['show_holes'],
			'bounds': [float(v) for pair in zip(lower, upper) for v in pair]}
		# The sliders' own values, which the page moves while they are dragged
		self.state.update({'slicer_origin_{}'.format(i): v for i, v in enumerate(self.state.slicer['origin'])})
		self.state.update({'slicer_normal_{}'.format(i): v for i, v in enumerate(self.state.slicer['normal'])})
		self.state.slicer_tolerance = self.state.slicer['tolerance']
		self.state.vertical_exaggeration = self.viewer._vertical_exaggeration
		self.state.background = self.viewer._background.capitalize()
		viewer = self.viewer
		self.state.cmap_items = self._cmap_items()
		self.state.update({'outline_visible': viewer._outline, 'axis_arrows': viewer._axis_arrows,
						   'axes_visible': viewer._corner_axes, 'ruler_visible': self._ruler_shown(),
						   'tick_note': viewer._tick_note,
						   'parallel_projection': bool(viewer.plotter.camera.parallel_projection)})
		spans = (_corners(viewer._extents)[1] - _corners(viewer._extents)[0]) if viewer._extents else None
		for key, chosen, span in zip(('xy', 'z'), viewer._tick_spacing,
									 (None if spans is None else spans[:2].max(), None if spans is None else spans[2])):
			self.state['tick_spacing_' + key] = 'Auto' if chosen is None else chosen
			self.state['tick_options_' + key] = _tick_options(span, chosen)
		self.state.axes = {'has_ruler': self.viewer._ruler is not None,
						   'origin': _floats(self.viewer._origin_values()),
						   'from_zero': self.viewer._axes_from_zero}
		if not self._paused:
			self._view.update()

	@contextmanager
	def _one_frame(self):
		"""Change the scene in several steps and send the page only the result. trame_pyvista draws
		the page again on every render of the plotter, and pyvista renders inside several of its
		camera calls, so a button used to send half-done frames first (no ruler, not yet fitted).
		pyvista has no public way to pause that, so its set of callbacks is swapped out."""
		plotter = self.viewer.plotter
		callbacks = plotter._on_render_callbacks
		plotter._on_render_callbacks = set()
		self._paused = True
		try:
			yield
		finally:
			plotter._on_render_callbacks = callbacks
			self._paused = False
		self.refresh()

	def _enable_picking(self):
		"""Let a click on a drill hole show its values. The picker is pyvista's, which calls
		``_on_pick`` with what was hit. A left click picks too, and a drag does not: a left press
		also starts a rotation, so the pick is made on release, if the mouse hardly moved. The
		release is caught as the end of an interaction, because trame's mouse handler sends no
		button release event"""
		plotter = self.viewer.plotter
		# A panel made before this one on the same viewer picks no more
		plotter.disable_picking()
		for tag in self.viewer._pick_observers:
			plotter.iren.remove_observer(tag)
		# The picker's own tolerance would snap a click that misses to a hole 40 to 80 pixels
		# away, so ``_on_pick`` measures the distance on the screen itself
		plotter.enable_point_picking(callback=self._on_pick, picker='cell', use_picker=True, show_message=False,
									 show_point=False, tolerance=0.05, pickable_window=True)
		pressed = {}

		def on_press(iren, event):
			pressed['at'] = iren.GetEventPosition()

		def on_release(iren, event):
			start = pressed.pop('at', None)
			x, y = iren.GetEventPosition()
			if start and max(abs(x - start[0]), abs(y - start[1])) <= _CLICK_SLOP:
				iren.GetPicker().Pick(x, y, 0, plotter.renderer)
		self.viewer._pick_observers = [plotter.iren.add_observer('LeftButtonPressEvent', on_press),
									   plotter.iren.add_observer('EndInteractionEvent', on_release)]

	def _click_distance(self, picker):
		"""The distance on the screen from the click to the picked cell (0 on the cell). It is measured
		to the corners of the cell: the position the picker reports is on the click's own line of
		sight, which made every distance 0"""
		import numpy as np
		renderer = self.viewer.plotter.renderer
		actor = picker.GetActor()
		corners = actor.mapper.dataset.get_cell(picker.GetCellId()).points * np.array(actor.GetScale())
		screen = []
		for corner in corners:
			renderer.SetWorldPoint(*corner, 1.0)
			renderer.WorldToDisplay()
			screen.append(renderer.GetDisplayPoint()[:2])
		screen = np.array(screen)
		click = np.array(picker.GetSelectionPoint()[:2])
		outside = np.maximum(screen.min(axis=0) - click, 0) + np.maximum(click - screen.max(axis=0), 0)
		return float(np.hypot(*outside))

	def _on_pick(self, point, picker):
		"""Something was clicked: if it is a drill hole, say what it is and mark it, else say nothing
		was picked. pyvista calls this for a click that hits nothing too, so that this can say so"""
		actor, cell = picker.GetActor(), picker.GetCellId()
		found = None
		if actor is not None and cell >= 0 and self._click_distance(picker) <= _PICK_SLACK:
			for layer in self.viewer._layers:
				if layer.kind == 'holes' and actor in (layer.actor, layer.crossing_actor):
					found = (layer, actor is layer.crossing_actor)
		self._clear_pick()
		if found is None:
			self.state.pick_info = _PICK_MISSED
		else:
			layer, on_slice = found
			self._show_pick(layer, self.viewer.interval_info(layer.name, cell, on_slice), on_slice)
		self.state.flush()
		self._view.update()

	def _show_pick(self, layer, info, on_slice):
		"""Write what is known of the interval in the page, and mark it with a wider tube. The mark is
		not pickable, or the next click on or near the hole would hit it"""
		import pyvista as pv
		values = ', '.join('{} = {}'.format(name, 'no value' if value is None else
											 '{:.2f}'.format(value) if isinstance(value, float) else value)
						   for name, value in info['values'].items())
		where = '' if info['hole'] is None else 'hole {}, '.format(info['hole'])
		self.state.pick_info = '{}{}: {}elevation {:.1f} to {:.1f}, {}'.format(
			layer.name, ' (on the slice)' if on_slice else '', where, info['upper'][2], info['lower'][2], values)
		mark = pv.Line(info['upper'], info['lower']).tube(radius=layer.radius * 5 / 3, n_sides=20)
		self.viewer.plotter.add_mesh(mark, name='pick highlight', color='magenta', opacity=0.6, pickable=False)

	def _clear_pick(self):
		"""Take the mark and the text of a picked interval away. They are of what was drawn when it was
		picked, so they go when anything is drawn again"""
		self.viewer.plotter.remove_actor('pick highlight', render=False)
		self.state.pick_info = _PICK_DEFAULT

	def _layer(self, name):
		return self.viewer._layers[self.viewer.layers.index(name)]

	def _list_views(self):
		"""The saved views in the folder, for the list in the page"""
		import os
		folder = self.views_folder
		self.state.view_files = sorted(f[:-5] for f in os.listdir(folder) if f.endswith('.json')) \
			if os.path.isdir(folder) else []

	def _guard(self, action):
		"""Call ``action()``. A value that does not work is explained in the panel, and the page
		goes back to the viewer as it is"""
		try:
			self.state.note = ''
			action()
			self.state.message = ''
		except ValueError as error:
			self.state.message = str(error)
			self.refresh()

	def _try(self, name, changes):
		"""``Viewer.set_layer(name, **changes())``, guarded"""
		self._guard(lambda: self.viewer.set_layer(name, **changes()))

	def set_visible(self, name, visible):
		"""Show or hide a layer, like ``Viewer.set_layer(name, visible=...)``"""
		self.viewer.set_layer(name, visible=visible)

	def set_opacity(self, name, opacity):
		"""Set the opacity of a layer, like ``Viewer.set_layer(name, opacity=...)``"""
		self.viewer.set_layer(name, opacity=opacity)

	def _cmap_items(self):
		"""The list of colormaps in the page: groups under their names, each colormap with a preview
		bar drawn from the colours pyvista gives it. A colormap in the scene that is not in a group
		is added, so that it can be shown as the chosen one"""
		from ..plotting.utils import _get_cmap

		def entry(name):
			if name not in self._gradients:
				cmap = scene_cmaps.get(name) or _get_cmap(name, False, None)
				self._gradients[name] = _gradient(cmap)
			return {'title': name, 'value': name, 'gradient': self._gradients[name]}
		groups = dict(_CMAP_GROUPS)
		known = {name for names in groups.values() for name in names}
		scene_cmaps = {layer.cmap: layer.kwargs['cmap'] for layer in self.viewer._layers
					   if layer.clim is not None and layer.cmap not in known}
		if scene_cmaps:
			groups = {'In the scene': sorted(scene_cmaps), **groups}
		items = []
		for group, names in groups.items():
			items += [{'type': 'subheader', 'title': group}] + [entry(name) for name in names]
		return items

	def _row(self, layer):
		"""What the page shows of a layer: :func:`_row` and the range of its tube radius slider, from a
		quarter to four times the radius the layer first had"""
		row = _row(layer)
		source = layer.source
		row['color_by'] = layer.color_var
		row['color_options'] = [] if source is None else [
			{'title': name if catdict is None else name + ' (categories)', 'value': name}
			for name, catdict in source.colourings.items()]
		row['cmap_kind'] = 'continuous' if layer.clim is not None else 'categorical' if layer.categorical else None
		names = list(layer.kwargs['annotations'].values()) if layer.categorical else []
		row['category_names'] = names
		row['category_colours'] = [mpl.colors.to_hex(c) for c in layer.kwargs['cmap'].colors] if names else []
		row['categories'] = list(layer.categories)
		row['palette_items'] = self._palette_items(len(layer.kwargs['annotations'])) if layer.categorical else []
		row['radius'] = layer.radius
		if layer.radius is not None:
			first = self._radius_ranges.setdefault(layer.name, layer.radius)
			row['radius_range'] = [first / 4, max(4 * first, layer.radius), first / 4]
		return row

	def set_radius(self, name, radius):
		"""Set the radius of the tubes of a drill hole layer, like ``Viewer.set_layer(name, radius=...)``"""
		self._try(name, lambda: {'radius': _number(radius)})

	def _palette_items(self, n):
		"""The palettes that have colours for ``n`` categories, in groups, each with a preview of its
		colours as blocks. A palette that is the layer's own but not in a group is not listed"""
		from ..plotting.cmaps import cat_palettes
		from ..plotting.utils import _get_cmap
		items = []
		for group, names in _PALETTE_GROUPS.items():
			entries = []
			for palette in names:
				if palette in cat_palettes and len(cat_palettes[palette]) < n:
					continue
				key = (palette, n)
				if key not in self._gradients:
					colours = [mpl.colors.to_hex(c) for c in _get_cmap(palette, True, n).colors]
					width = 100 / n
					stops = ['{} {:.2f}% {:.2f}%'.format(c, i * width, (i + 1) * width) for i, c in enumerate(colours)]
					self._gradients[key] = 'linear-gradient(to right, {})'.format(', '.join(stops))
				entries.append({'title': palette, 'value': palette, 'gradient': self._gradients[key]})
			if entries:
				items += [{'type': 'subheader', 'title': group}] + entries
		return items

	def toggle_category(self, name, category):
		"""Isolate buttons: the first click shows only that category, the next ones add or remove
		others. Every category shown is the same as none chosen, which shows all"""
		names = list(self._layer(name).kwargs['annotations'].values())
		chosen = set(self._layer(name).categories) ^ {category}
		self._try(name, lambda: {'categories': [] if len(chosen) == len(names) else [n for n in names if n in chosen]})

	def show_all_categories(self, name):
		"""Show every category of a layer again"""
		self._try(name, lambda: {'categories': []})

	def set_color_by(self, name, variable):
		"""Colour a layer by another variable, like ``Viewer.set_layer(name, color_by=...)``"""
		self._try(name, lambda: {'color_by': variable})

	def set_cmap(self, name, cmap):
		"""Set the colormap of a layer, like ``Viewer.set_layer(name, cmap=...)"""
		self._try(name, lambda: {'cmap': cmap})

	def set_edges(self, name, edges):
		"""Show or hide the cell edges of a block layer, like ``Viewer.set_layer(name, edges=...)"""
		self.viewer.set_layer(name, edges=edges)

	def set_clim(self, name, end, value):
		"""Set one end of the colour limits of a layer (``end`` 0 for min, 1 for max), like
		``Viewer.set_layer(name, clim=...)``"""
		def changes():
			clim = list(self._layer(name).clim)
			clim[end] = _number(value)
			return {'clim': clim}
		self._try(name, changes)

	def _threshold(self, name):
		"""The threshold a layer shows with its typed minimum: from the minimum up to the maximum
		the layer already has (from a script), else the highest value of the data. The minimum is
		the lowest of the data until one is typed"""
		layer = self._layer(name)
		low, high = layer.source.mesh.get_data_range(layer.source.var)
		if layer.threshold is not None:
			high = layer.threshold[1]
		return (self._typed_min.get(name, low), high)

	def set_threshold_on(self, name, on):
		"""Show only the values from the typed minimum up (on), or all values (off), like
		``Viewer.set_layer(name, threshold=...)``"""
		self._try(name, lambda: {'threshold': self._threshold(name) if on else False})

	def set_threshold_min(self, name, value):
		"""Type the lowest value to show. It is applied when the threshold is on, and kept
		when it is not; a value that does not work is explained and not kept"""
		def action():
			number = _number(value)
			if self._layer(name).threshold is not None:
				self.viewer.set_layer(name, threshold=(number, self._threshold(name)[1]))
			self._typed_min[name] = number
		self._guard(action)
		self.refresh()

	def set_slicer_enabled(self, enabled):
		"""Switch the slicer on or off, like ``Viewer.set_slicer(enabled=...)``"""
		self._guard(lambda: self.viewer.set_slicer(enabled=enabled))

	def _on_slicer_sliders(self, **changed):
		"""The slicer's sliders moved: cut again with their values, as the slider is dragged"""
		state, slicer = self.state, self.viewer._slicer
		origin = tuple(float(state['slicer_origin_{}'.format(i)]) for i in range(3))
		normal = tuple(float(state['slicer_normal_{}'.format(i)]) for i in range(3))
		tolerance = float(state.slicer_tolerance)
		if (origin, normal, tolerance) != (tuple(self._slicer_origin()), slicer['normal'], slicer['tolerance']):
			self._guard(lambda: self.viewer.set_slicer(origin=origin, normal=normal, tolerance=tolerance))

	def _on_vertical_exaggeration(self, vertical_exaggeration, **changed):
		"""The slider moved: stretch the scene with its value, as the slider is dragged"""
		factor = float(vertical_exaggeration)
		if factor != self.viewer._vertical_exaggeration:
			self._guard(lambda: self.viewer.set_vertical_exaggeration(factor))

	def _on_background(self, background, **changed):
		"""The background was chosen: change it, like ``Viewer.set_background``"""
		if background.lower() != self.viewer._background:
			self._guard(lambda: self.viewer.set_background(background))

	def _ruler_shown(self):
		"""Whether the scene has a ruler that is drawn"""
		return self.viewer._ruler is not None and self.viewer._ruler_visible

	def _on_scene(self, **changed):
		"""A box of the scene was ticked or unticked: change what the viewer draws around the
		layers, like ``Viewer.set_scene``"""
		viewer, state = self.viewer, self.state
		wanted = {'outline': state.outline_visible, 'ruler': state.ruler_visible,
				  'axis_arrows': state.axis_arrows, 'corner_axes': state.axes_visible}
		current = {'outline': viewer._outline, 'ruler': self._ruler_shown(),
				   'axis_arrows': viewer._axis_arrows, 'corner_axes': viewer._corner_axes}
		differ = {key: bool(value) for key, value in wanted.items() if bool(value) != current[key]}
		if differ:
			self._guard(lambda: viewer.set_scene(**differ))

	def _on_tick_spacing(self, **changed):
		"""A tick spacing was chosen: like ``Viewer.set_axis_ticks``"""
		chosen = tuple(None if self.state['tick_spacing_' + key] == 'Auto' else float(self.state['tick_spacing_' + key])
					   for key in ('xy', 'z'))
		if chosen != tuple(self.viewer._tick_spacing):
			self._guard(lambda: self.viewer.set_axis_ticks(*chosen))

	def show_view(self, view):
		"""Look from a direction, like ``Viewer.set_camera_view``"""
		with self._one_frame():
			self._guard(lambda: self.viewer.set_camera_view(view))

	def fit_view(self):
		"""Frame the whole scene, like ``Viewer.fit_view``"""
		with self._one_frame():
			self.viewer.fit_view()

	def attach_screenshot(self):
		"""The scene as the page shows it, as a picture, for the browser to download"""
		return self.server.protocol.addAttachment(self._view.screenshot())

	def attach_export(self):
		"""The scene as a web page, like ``Viewer.export_html``, for the browser to download"""
		try:
			page = self.viewer._export_page()
		finally:
			self._view.update()
		return self.server.protocol.addAttachment(page)

	def _on_parallel_projection(self, parallel_projection, **changed):
		"""The box was ticked: like ``Viewer.set_parallel_projection``"""
		if bool(parallel_projection) != bool(self.viewer.plotter.camera.parallel_projection):
			with self._one_frame():
				self.viewer.set_parallel_projection(bool(parallel_projection))

	def snap_slicer(self, orient):
		"""Face the plane along an axis: 'xy' (plan, normal z), 'xz' (normal y) or 'yz' (normal x)"""
		normal = {'xy': (0.0, 0.0, 1.0), 'xz': (0.0, 1.0, 0.0), 'yz': (1.0, 0.0, 0.0)}[orient]
		self._guard(lambda: self.viewer.set_slicer(normal=normal))

	def set_show_holes(self, show):
		"""Show the drill hole intervals that cross the plane"""
		self._guard(lambda: self.viewer.set_slicer(show_holes=show))

	def _slicer_origin(self):
		"""The origin of the plane: where it is, or the centre of the scene before it was moved"""
		origin = self.viewer._slicer['origin']
		if origin is None:
			lower, upper = _corners(self.viewer._extents)
			origin = (lower + upper) / 2
		return origin

	def set_origin(self, axis, value):
		"""Set where the ruler counts from along one axis (0 east, 1 north, 2 elevation), like
		``Viewer.set_axis_origin(origin)``"""
		def action():
			origin = list(self.viewer._axis_origin or (0.0, 0.0, 0.0))
			origin[axis] = _number(value)
			self.viewer.set_axis_origin(origin, self.viewer._axes_from_zero)
		self._guard(action)

	def set_from_zero(self, from_zero):
		"""Count every axis from 0 at the lowest edge of the scene, like
		``Viewer.set_axis_origin(from_zero=...)``"""
		self.viewer.set_axis_origin(self.viewer._axis_origin, from_zero)

	def _view_path(self, name):
		"""The file of a saved view. The name is only a name, never a path"""
		import os
		import re
		if not re.fullmatch(r'[\w\- ]+', name):
			raise ValueError("A view name can only have letters, digits, spaces, - and _")
		return os.path.join(self.views_folder, name + '.json')

	def save_view(self, name):
		"""Save the view (see ``Viewer.get_view``) as ``<views_folder>/<name>.json``. Without a
		name, one is made from the time. Letters, digits, spaces, - and _ only"""
		def action():
			import json
			import os
			from datetime import datetime
			view_name = (name or '').strip() or datetime.now().strftime('view %Y-%m-%d %H%M%S')
			path = self._view_path(view_name)
			os.makedirs(self.views_folder, exist_ok=True)
			with open(path, 'w') as f:
				json.dump(self.viewer.get_view(), f, indent=1)
			self._list_views()
			self.state.view_file = view_name
			self.state.note = 'Saved ' + path
		self._guard(action)

	def load_view(self, name):
		"""Show a view saved by :meth:`save_view`, like ``Viewer.set_view``"""
		def action():
			import json
			path = self._view_path(name)
			try:
				with open(path) as f:
					view = json.load(f)
				self.viewer.set_view(view)
			except (OSError, KeyError, TypeError) as error:
				raise ValueError("Could not load {}: {}".format(path, error)) from None
			self.state.note = 'Loaded ' + path
		self._guard(action)


def _tick_options(span, chosen):
	"""The tick spacings to choose from for an axis ``span`` long, as the page's list: 'Auto' and the
	round spacings (1, 2 or 5 times a power of 10) that give 3 to 10 labels, and the ``chosen``
	one when a script chose another"""
	steps = [] if span is None else [f * 10.0 ** p for p in range(-3, 10) for f in (1, 2, 5)
									 if 2 <= span / (f * 10.0 ** p) <= _MAX_LABELS - 1]
	if chosen is not None and chosen not in steps:
		steps = sorted(steps + [chosen])
	return [{'title': 'Auto', 'value': 'Auto'}] + [{'title': '{:g}'.format(step), 'value': step} for step in steps]


def _number(value):
	"""A number typed into a box"""
	try:
		return float(value)
	except ValueError:
		raise ValueError("'{}' is not a number".format(value)) from None


def _floats(values):
	"""A list of plain floats, which the page can be given, or None"""
	return None if values is None else [float(v) for v in values]


def _row(layer):
	"""What the page shows of a layer. The page is given plain numbers and strings"""
	source = layer.source
	threshold_range = None
	if source is not None and source.var is not None:
		threshold_range = _floats(source.mesh.get_data_range(source.var))
	return {'name': layer.name, 'visible': layer.visible, 'opacity': layer.kwargs.get('opacity', 1.0),
			'cmap': layer.cmap, 'edges': bool(layer.kwargs.get('show_edges', False)) if layer.kind == 'blocks' else None,
			'clim': _floats(layer.clim), 'data_range': _floats(layer.clim and layer.data_range),
			'variable': source and source.var, 'threshold_range': threshold_range,
			'threshold_on': layer.threshold is not None,
			'threshold_min': None if layer.threshold is None else float(layer.threshold[0])}
