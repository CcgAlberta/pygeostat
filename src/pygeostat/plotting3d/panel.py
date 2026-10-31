# !/usr/bin/env python
#  -*- coding: utf-8 -*-
"""panel is the interactive layer panel of a Viewer, a web page (trame) with the controls beside the 3D view"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


class Panel:
	"""
	The controls of a :class:`~pygeostat.plotting3d.viewer.Viewer`, shown beside its 3D view.
	Made by :meth:`Viewer.panel() <pygeostat.plotting3d.viewer.Viewer.panel>`.

	Every control calls a method of the viewer (e.g. ``set_layer``), so what the panel does can
	also be done in a script.

	Parameters:
		viewer (gs.Viewer): the scene to control

	Attributes:
		layout: the page. A panel shows itself when it is the last line of a notebook cell
	"""

	def __init__(self, viewer):
		import uuid
		try:
			from trame.app import get_server
			from trame.ui.vuetify3 import SinglePageWithDrawerLayout
			from trame.widgets import html, vuetify3
			from trame_pyvista.ui import get_viewer
		except ImportError:
			raise ImportError("The layer panel requires trame: pip install pygeostat[view3d-panel]") from None
		self.viewer = viewer
		# A new name each time: get_server caches by name, so a fixed one would share state
		self.server = get_server('pygeostat_{}'.format(uuid.uuid4().hex[:8]), client_type='vue3')
		self.state = self.server.state
		self.state.message = ''
		self._view = get_viewer(viewer.plotter)
		self.refresh()

		with SinglePageWithDrawerLayout(self.server) as self.layout:
			self.layout.title.set_text('Layers')
			with self.layout.drawer:
				html.Div('{{ message }}', v_if='message', classes='text-caption text-error ma-2')
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
						with html.Div(v_if='layer.threshold_range', classes='mt-2'):
							with html.Div(classes='d-flex align-center'):
								html.Div('Show values between', classes='text-caption')
								vuetify3.VBtn('All', size='x-small', variant='text', classes='ml-auto',
											  disabled=('!layer.threshold',),
											  click=(self.clear_threshold, '[layer.name]'))
							with html.Div(classes='d-flex ga-2'):
								self._range_boxes('(layer.threshold || layer.threshold_range)', 'layer.threshold_range',
												  self.set_threshold)
						with html.Div(v_if='layer.can_slice', classes='mt-2'):
							vuetify3.VCheckbox(label='Slice', model_value=('!!layer.slice',),
											   update_modelValue=(self.set_slice_on, '[layer.name, $event]'),
											   density='compact', hide_details=True)
							with html.Div(v_if='layer.slice'):
								vuetify3.VSelect(
									items=([{'title': 'XY (plan)', 'value': 'xy'}, {'title': 'XZ (section)', 'value': 'xz'},
											 {'title': 'YZ (section)', 'value': 'yz'}],), model_value=('layer.slice.orient',),
									update_modelValue=(self.set_orient, '[layer.name, $event]'), density='compact',
									hide_details=True)
								with html.Div(classes='d-flex ga-2 mt-1'):
									for key, label, handler in (('position', 'layer.slice.position_label', self.set_position),
																('width', 'layer.slice.width_label', self.set_width)):
										vuetify3.VTextField(
											label=(label,), type='number', model_value=('layer.slice.' + key,),
											density='compact', hide_details=True, variant='outlined',
											change=(handler, '[layer.name, $event.target.value]'))
			with self.layout.content:
				with vuetify3.VContainer(fluid=True, classes='pa-0 fill-height'):
					self._view.ui(add_menu=False)

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
		from IPython.display import display
		display(self.layout)

	def refresh(self):
		"""Show the viewer's layers as they are now. The viewer calls this when a layer is added or
		changed, so the page follows a script."""
		self.state.layers = [_row(layer) for layer in self.viewer._layers]
		self._view.update()

	def _layer(self, name):
		return self.viewer._layers[self.viewer.layers.index(name)]

	def _try(self, name, changes):
		"""``Viewer.set_layer(name, **changes())``. A value that does not work is explained in the
		panel, and the page goes back to the layer as it is"""
		try:
			self.viewer.set_layer(name, **changes())
			self.state.message = ''
		except ValueError as error:
			self.state.message = str(error)
			self.refresh()

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

	def set_threshold(self, name, end, value):
		"""Set one end of the threshold of a layer (``end`` 0 for min, 1 for max), like
		``Viewer.set_layer(name, threshold=...)``. A layer without one starts from the range of its
		variable"""
		def changes():
			layer = self._layer(name)
			threshold = list(layer.threshold or layer.source.mesh.get_data_range(layer.source.var))
			threshold[end] = _number(value)
			return {'threshold': threshold}
		self._try(name, changes)

	def clear_threshold(self, name):
		"""Show all values of a layer again"""
		self._try(name, lambda: {'threshold': False})

	def set_slice_on(self, name, on):
		"""Show only a slice of a layer, in the middle of the layer, or all of it again"""
		if on:
			self._try(name, lambda: _middle_slice(self._layer(name), 'xz'))
		else:
			self._try(name, lambda: {'orient': False})

	def set_orient(self, name, orient):
		"""Change the orientation of the slice of a layer, and put the slice in the middle"""
		self._try(name, lambda: _middle_slice(self._layer(name), orient))

	def set_position(self, name, value):
		"""Move the slice of a layer: to a slice number if it is set by one, else to a coordinate"""
		def changes():
			by_number = self._layer(name).section['slice_number'] is not None
			return {'slice_number': int(_number(value))} if by_number else {'slice_coordinate': _number(value)}
		self._try(name, changes)

	def set_width(self, name, value):
		"""Set how thick the slice of a layer is: extra slices each side for data on a grid, else
		the distance each side"""
		def changes():
			gridded = self._layer(name).source.griddef is not None
			return {'slice_thickness': int(_number(value))} if gridded else {'pointtol': _number(value)}
		self._try(name, changes)


def _number(value):
	"""A number typed into a box"""
	try:
		return float(value)
	except ValueError:
		raise ValueError("'{}' is not a number".format(value)) from None


def _middle_slice(layer, orient):
	"""The arguments of set_layer for a slice of a layer in the middle of it along the axis not in
	``orient``, as thick as the slice shown now, or 0 or a twentieth of the layer when there is none"""
	axis = {'xy': 2, 'xz': 1, 'yz': 0}[orient]
	source, section = layer.source, layer.section or {}
	if source.griddef is not None:
		size = (source.griddef.nx, source.griddef.ny, source.griddef.nz)[axis]
		return {'orient': orient, 'slice_number': size // 2,
				'slice_thickness': section.get('slice_thickness') or 0}
	lower, upper = source.mesh.bounds[2 * axis:2 * axis + 2]
	return {'orient': orient, 'slice_coordinate': (lower + upper) / 2,
			'pointtol': section.get('pointtol') or (upper - lower) / 20}


def _floats(values):
	"""A list of plain floats, which the page can be given, or None"""
	return None if values is None else [float(v) for v in values]


def _row(layer):
	"""What the page shows of a layer. The page is given plain numbers and strings"""
	source, section = layer.source, layer.section
	threshold_range = None
	if source is not None and source.var is not None:
		threshold_range = _floats(source.mesh.get_data_range(source.var))
	row = {'name': layer.name, 'visible': layer.visible, 'opacity': layer.kwargs.get('opacity', 1.0),
		   'clim': _floats(layer.clim), 'data_range': _floats(layer.clim and layer.data_range),
		   'threshold_range': threshold_range, 'threshold': _floats(layer.threshold),
		   'can_slice': source is not None, 'slice': None}
	if section is not None:
		by_number = section['slice_number'] is not None
		gridded = source.griddef is not None
		row['slice'] = {
			'orient': section['orient'],
			'position': section['slice_number'] if by_number else section['slice_coordinate'],
			'position_label': 'Slice number' if by_number else 'Coordinate',
			'width': (section['slice_thickness'] or 0) if gridded else section['pointtol'],
			'width_label': 'Extra slices each side' if gridded else 'Distance each side'}
	return row
