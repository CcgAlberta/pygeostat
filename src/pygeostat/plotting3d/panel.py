#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""panel is the interactive layer panel of a Viewer, a web page (trame) with the controls beside the 3D view"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from .viewer import _corners


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
		# The lowest value typed for the threshold of each layer, kept while the threshold is off
		self._typed_min = {}
		# A new name each time: get_server caches by name, so a fixed one would share state
		self.server = get_server('pygeostat_{}'.format(uuid.uuid4().hex[:8]), client_type='vue3')
		self.state = self.server.state
		self.state.message = ''
		self.state.note = ''
		self.state.view_name = ''
		self.state.view_file = None
		self._list_views()
		# trame_pyvista keeps one viewer per plotter, tied to the server it was made with, so one
		# left from an earlier panel of this viewer would draw on that panel's server
		import trame_pyvista.ui
		trame_pyvista.ui._VIEWERS.pop(viewer.plotter._id_name, None)
		self._view = get_viewer(viewer.plotter, server=self.server)
		self.refresh()
		self.state.change(*['slicer_origin_{}'.format(i) for i in range(3)],
						  *['slicer_normal_{}'.format(i) for i in range(3)], 'slicer_tolerance')(
			self._on_slicer_sliders)
		self.state.change('vertical_exaggeration')(self._on_vertical_exaggeration)

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
						self._slider('Vertical exaggeration', 'vertical_exaggeration', 1, 10, step=0.5)
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
		self.state.layers = [_row(layer) for layer in self.viewer._layers]
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
		self.state.axes = {'has_ruler': self.viewer._ruler is not None,
						   'origin': _floats(self.viewer._origin_values()),
						   'from_zero': self.viewer._axes_from_zero}
		self._view.update()

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
			'clim': _floats(layer.clim), 'data_range': _floats(layer.clim and layer.data_range),
			'variable': source and source.var, 'threshold_range': threshold_range,
			'threshold_on': layer.threshold is not None,
			'threshold_min': None if layer.threshold is None else float(layer.threshold[0])}
