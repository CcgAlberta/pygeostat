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
						with html.Div(v_if='layer.clim', classes='d-flex ga-2 mt-2'):
							for end, label in enumerate(('Min', 'Max')):
								# The data's own range in the label, the limits typed: a value tied to the data is typed
								vuetify3.VTextField(
									label=("'{} (data ' + Number(layer.data_range[{}]).toPrecision(3) + ')'".format(label, end),),
									type='number', model_value=('layer.clim[{}]'.format(end),), density='compact',
									hide_details=True, variant='outlined',
									change=(self.set_clim, "[layer.name, {}, $event.target.value]".format(end)))
			with self.layout.content:
				with vuetify3.VContainer(fluid=True, classes='pa-0 fill-height'):
					self._view.ui(add_menu=False)

	def _ipython_display_(self):
		from IPython.display import display
		display(self.layout)

	def refresh(self):
		"""Show the viewer's layers as they are now. The viewer calls this when a layer is added or
		changed, so the page follows a script."""
		self.state.layers = [{'name': layer.name, 'visible': layer.visible,
							  'opacity': layer.kwargs.get('opacity', 1.0),
							  'clim': layer.clim and [float(v) for v in layer.clim],
							  'data_range': layer.clim and [float(v) for v in layer.data_range]}
							 for layer in self.viewer._layers]
		self._view.update()

	def set_visible(self, name, visible):
		"""Show or hide a layer, like ``Viewer.set_layer(name, visible=...)``"""
		self.viewer.set_layer(name, visible=visible)

	def set_opacity(self, name, opacity):
		"""Set the opacity of a layer, like ``Viewer.set_layer(name, opacity=...)``"""
		self.viewer.set_layer(name, opacity=opacity)

	def set_clim(self, name, end, value):
		"""Set one end of the colour limits of a layer (``end`` 0 for min, 1 for max), like
		``Viewer.set_layer(name, clim=...)``. A value that does not work is explained in the panel"""
		clim = list(self.viewer._layers[self.viewer.layers.index(name)].clim)
		try:
			clim[end] = float(value)
			self.viewer.set_layer(name, clim=clim)
			self.state.message = ''
		except ValueError as error:
			self.state.message = str(error)
			self.refresh()
