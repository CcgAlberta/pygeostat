# !/usr/bin/env python
#  -*- coding: utf-8 -*-
"""view3d shows drill holes, block models and surfaces together in an interactive 3D scene"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


class _Source:
	"""The whole mesh of a layer, with what is needed to cut a part of it, so the part shown can
	change after the layer is drawn"""

	def __init__(self, mesh, var, griddef):
		self.mesh = mesh
		self.var = var
		self.griddef = griddef

	def cut(self, threshold, section, plane=None):
		"""The part of the mesh in the threshold (min, max of the variable) and in the slice
		(a dictionary of the slice arguments of view3d), whole if either is None. With a plane
		(a dictionary with its ``origin`` and ``normal``), only the cut through that plane"""
		mesh = self.mesh
		if threshold is not None:
			if self.var is None:
				raise ValueError("`threshold` needs a variable `var`")
			mesh = mesh.threshold(threshold, scalars=self.var)
		if section is not None:
			mesh = _slice(mesh, self.griddef, **section)
		if plane is not None:
			mesh = mesh.slice(normal=plane['normal'], origin=plane['origin'])
		return mesh


def _actor_name(layer_name):
	"""The name a layer's actor has in the plotter, so that a layer can be drawn again in its place"""
	return 'layer ' + layer_name


class _Layer:
	"""A mesh in the scene with the options it was drawn with, which the panel changes later"""

	def __init__(self, name, mesh, kwargs, actor, data_range=None, source=None, threshold=None,
				 section=None, kind='blocks'):
		self.name = name
		# 'blocks' (or any mesh) is cut by the slicer's plane, 'holes' are drill hole intervals and
		# 'points' are samples, which the slicer does not cut (see set_slicer)
		self.kind = kind
		self.mesh = mesh
		self.kwargs = kwargs
		self.actor = actor
		self.visible = True
		# The drill hole intervals crossing the slicer's plane, with their actor (see set_slicer)
		self.crossing = None
		self.crossing_actor = None
		self.source = source
		self.threshold = threshold
		self.section = section
		# The range of the variable over the whole layer, not only the part shown after a
		# threshold or slice
		if data_range is None and 'scalars' in kwargs:
			data_range = mesh.get_data_range(kwargs['scalars'])
		self.data_range = data_range

	@property
	def key(self):
		return _actor_name(self.name)

	@property
	def clim(self):
		"""The colour limits (min, max) of a layer coloured by a continuous variable, else None"""
		if 'scalars' not in self.kwargs or 'annotations' in self.kwargs:
			return None
		return tuple(self.kwargs['clim'] or self.data_range)


class Viewer:
	"""
	A 3D scene that layers are added to with :func:`gs.view3d() <pygeostat.plotting3d.viewer.view3d>`.

	Requires PyVista (``pip install pyvista``).

	Parameters:
		off_screen (bool): render without opening a window, e.g. for tests or saved images

	Attributes:
		plotter (pyvista.Plotter): the PyVista plotter behind the scene, for anything not covered
			here
	"""

	def __init__(self, off_screen=None):
		try:
			import pyvista as pv
		except ImportError:
			raise ImportError("Viewer requires PyVista: pip install pyvista") from None
		self.plotter = pv.Plotter(off_screen=off_screen)
		# The layers, so show_static can rebuild the scene
		self._layers = []
		# The extent of each layer before any threshold, which the ruler wraps
		self._extents = []
		self._ruler = None
		self._axis_origin = None
		self._slicer = {'enabled': False, 'origin': None, 'normal': (0.0, 0.0, 1.0), 'tolerance': 10.0,
						'show_holes': False}
		self._axes_from_zero = False
		# The place of the colour bar of each title, in the order they were first needed
		self._bar_slots = {}
		self._vertical_exaggeration = 1.0
		# The panel made by panel(), which shows the layers
		self._panel = None

	@property
	def layers(self):
		"""Names of the layers, in the order they were added (list of str)"""
		return [layer.name for layer in self._layers]

	def add(self, mesh, extent=None, name=None, data_range=None, source=None, threshold=None,
			section=None, kind='blocks', **kwargs):
		"""
		Add a mesh to the scene. Keyword arguments go to ``pyvista.Plotter.add_mesh``.

		Parameters:
			mesh (pyvista.DataSet): the mesh
			extent (tuple): bounds (xmin, xmax, ymin, ymax, zmin, zmax) for the ruler to wrap,
				e.g. the whole model when ``mesh`` is part of it. If None, the bounds of ``mesh``
			name (str): name of the layer, to change it with :meth:`set_layer`. A name that is
				already taken gets a number added. If None, ``'Layer 1'``, ``'Layer 2'``, ...
			data_range (tuple): range (min, max) of the variable the mesh is coloured by over the
				whole layer, when ``mesh`` is only part of it. If None, the range of ``mesh``
			source (_Source): the whole mesh and what is needed to cut it, when ``mesh`` is a part
				of it that :meth:`set_layer` can change. ``threshold`` and ``section`` are how
				``mesh`` was cut from it
			kind (str): ``'holes'`` for drill hole intervals, ``'points'`` for samples, otherwise
				``'blocks'``; the slicer cuts a ``'blocks'`` layer, hides a ``'holes'`` layer and
				leaves a ``'points'`` layer as it is
		"""
		name = name or 'Layer {}'.format(len(self._layers) + 1)
		if 'scalars' in kwargs:
			kwargs['scalar_bar_args'] = self._placed_bar(kwargs['scalar_bar_args'])
		first, number = name, 1
		while name in self.layers:
			number += 1
			name = '{} ({})'.format(first, number)
		self._layers.append(_Layer(name, mesh, kwargs, self.plotter.add_mesh(mesh, name=_actor_name(name), **kwargs), data_range,
								   source, threshold, section, kind))
		self._extents.append(mesh.bounds if extent is None else extent)
		if self._vertical_exaggeration != 1:
			self.set_vertical_exaggeration(self._vertical_exaggeration)
		elif self._ruler is not None:
			self._draw_ruler(self.plotter)
		self._refresh_panel()

	def set_layer(self, name, visible=None, opacity=None, clim=None, threshold=None, orient=None,
				  slice_number=None, slice_coordinate=None, slice_thickness=None, pointtol=None):
		"""
		Change a layer that is already in the scene. The change shows in :meth:`show` and in
		:meth:`show_static`. The ruler keeps wrapping a hidden layer, so hiding one does not move
		the axes.

		Parameters:
			name (str): name of the layer, one of :attr:`layers`
			visible (bool): show or hide the layer. If None, unchanged
			opacity (float): 0 (invisible) to 1 (opaque). If None, unchanged
			clim (tuple): colour limits (min, max) of a layer coloured by a continuous variable.
				If None, unchanged
			threshold (tuple): show only values of the variable between (min, max), as in
				:func:`gs.view3d() <pygeostat.plotting3d.viewer.view3d>`; ``False`` shows all of
				the layer again. If None, unchanged. Colours keep the limits of the layer
			orient (str): show only a slice, ``'xy'``, ``'xz'`` or ``'yz'``; ``False`` shows all of
				the layer again. If None, the slice stays as it is
			slice_number (int): the slice of a layer on a grid, which stays as it is if None
			slice_coordinate (float): the position of the slice, which stays as it is if None
			slice_thickness (int): extra grid slices each side, which stays as it is if None
			pointtol (float): distance each side for a layer not on a grid, which stays as it is
				if None

		The slice arguments work as in :func:`gs.view3d() <pygeostat.plotting3d.viewer.view3d>`.
		Changing a slice that is already shown needs only the arguments that change. A layer
		made from a ready mesh (not from data) cannot be thresholded or sliced.
		"""
		matches = [layer for layer in self._layers if layer.name == name]
		if not matches:
			raise ValueError("No layer '{}'; the layers are {}".format(name, self.layers))
		layer = matches[0]
		# Everything is checked before anything changes, so a refused call changes nothing
		if opacity is not None and not 0 <= opacity <= 1:
			raise ValueError("`opacity` must be between 0 and 1")
		if clim is not None:
			if layer.clim is None:
				raise ValueError("Layer '{}' is not coloured by a continuous variable".format(name))
			if not clim[0] < clim[1]:
				raise ValueError("`clim` must be (min, max) with min below max")
		slice_args = (orient, slice_number, slice_coordinate, slice_thickness, pointtol)
		cut = threshold is not None or any(v is not None for v in slice_args)
		if cut:
			if layer.source is None:
				raise ValueError("Layer '{}' was added as a mesh, so it cannot be thresholded or sliced".format(name))
			new_threshold = layer.threshold if threshold is None else None if threshold is False else tuple(threshold)
			if new_threshold is not None and not new_threshold[0] <= new_threshold[1]:
				raise ValueError("`threshold` must be (min, max) with min not above max")
			new_section = layer.section
			if any(v is not None for v in slice_args):
				new_section = _new_section(layer.section, *slice_args)
			# The slicer's plane may leave nothing of a layer, which is not an error, but a
			# threshold or slice that leaves nothing is
			if layer.source.cut(new_threshold, new_section).n_cells == 0:
				raise ValueError("Nothing is left of layer '{}' with this threshold and slice".format(name))
			layer.threshold, layer.section = new_threshold, new_section
			self._place(layer)
		if opacity is not None:
			layer.kwargs['opacity'] = opacity
			if layer.actor is not None:
				layer.actor.prop.opacity = opacity
		if clim is not None:
			layer.kwargs['clim'] = tuple(clim)
			if layer.actor is not None:
				layer.actor.mapper.scalar_range = tuple(clim)
		if visible is not None:
			layer.visible = bool(visible)
			if layer.actor is not None:
				layer.actor.visibility = layer.visible
		if layer.kind == 'holes':
			self._place_crossing(layer)
		self._sync_scalar_bars()
		self._refresh_panel()

	def _place(self, layer):
		"""Draw a layer again from its source, with its threshold, slice and, for blocks, the
		slicer's plane. A plane that misses the layer leaves nothing to draw."""
		plane = self._plane() if layer.kind == 'blocks' else None
		layer.mesh = layer.source.cut(layer.threshold, layer.section, plane)
		if not layer.mesh.n_cells:
			if layer.actor is not None:
				self.plotter.remove_actor(layer.actor)
				layer.actor = None
			return
		# Replacing the actor by its name, as the prototype does, swaps the mesh at once and keeps
		# the colour bar: removing the actor first took its bar with it, and a new bar was made
		# on every move of the slicer, so the slice and its legend blinked
		layer.actor = self.plotter.add_mesh(layer.mesh, name=layer.key, **layer.kwargs)
		layer.actor.visibility = layer.visible
		self._keep_clim(layer)

	def _placed_bar(self, bar_args):
		"""The colour bar arguments with a position of its own, where PyVista would put the bar in
		its next free place. PyVista takes a bar down with the layer it was made for only if it
		chose the place itself, so a bar with a position stays while the slicer redraws the layer,
		instead of a new bar being made at every move"""
		import pyvista as pv
		from pyvista.plotting.scalar_bars import MAX_N_COLOR_BARS
		theme = pv.global_theme.colorbar_vertical
		slot = self._bar_slots.setdefault(bar_args['title'], len(self._bar_slots))
		if slot >= MAX_N_COLOR_BARS:
			raise RuntimeError("There is room for {} colour bars, name fewer variables or labels".format(
				MAX_N_COLOR_BARS))
		return {**bar_args, 'vertical': True, 'width': theme.width, 'height': theme.height,
				'position_x': theme.position_x - slot * (theme.width + 0.2 * theme.width),
				'position_y': theme.position_y}

	def _keep_clim(self, layer):
		"""Put the colour limits of a layer back on its colour bar. A bar that is reused is widened
		by pyvista to the widest range it has had, so narrowed limits would come back wider"""
		if layer.clim is not None:
			self.plotter.update_scalar_bar_range(layer.clim, name=layer.kwargs['scalar_bar_args']['title'])

	def add_ruler(self, unit=None):
		"""
		Draw axes with grid lines on the back walls around everything in the scene, titled like
		the 2D plots (e.g. "Easting (m)").

		Parameters:
			unit (str): unit in the axis titles. If None, based on gs.Parameters['plotting.unit'];
				use ``''`` for no unit
		"""
		from ..pygeostat_parameters import Parameters
		if unit is None:
			unit = Parameters['plotting.unit']
		unit = ' ({})'.format(unit) if unit else ''
		self._ruler = {'xtitle': Parameters['plotting.xname'] + unit,
					   'ytitle': Parameters['plotting.yname'] + unit,
					   'ztitle': Parameters['plotting.zname'] + unit}
		self._draw_ruler(self.plotter)
		self._refresh_panel()

	def set_slicer(self, enabled=None, origin=None, normal=None, tolerance=None, show_holes=None):
		"""
		Cut the scene with a plane. Block layers show only the cut through the plane; the plane
		can be tilted. (Drill holes are described at ``show_holes``.)

		Parameters:
			enabled (bool): switch the slicer on or off. If None, unchanged
			origin (tuple): a point (x, y, z) on the plane. If None, unchanged; at first the
				centre of the scene
			normal (tuple): the direction (nx, ny, nz) the plane faces, e.g. (0, 0, 1) for a plan
				view. If None, unchanged; at first (0, 0, 1). (0, 0, 0) is taken as (0, 0, 1)
			tolerance (float): thickness in data units of the band of drill holes shown around the
				plane (see ``show_holes``). If None, unchanged; at first 10
			show_holes (bool): show the drill hole intervals that cross the plane. If None,
				unchanged
		"""
		if origin is not None and len(origin) != 3:
			raise ValueError("`origin` must be (x, y, z)")
		if normal is not None and len(normal) != 3:
			raise ValueError("`normal` must be (nx, ny, nz)")
		if tolerance is not None and not tolerance >= 0:
			raise ValueError("`tolerance` must be 0 or more")
		slicer = self._slicer
		was_enabled = slicer['enabled']
		if origin is not None:
			slicer['origin'] = tuple(float(v) for v in origin)
		if normal is not None:
			normal = tuple(float(v) for v in normal)
			slicer['normal'] = normal if any(normal) else (0.0, 0.0, 1.0)
		if tolerance is not None:
			slicer['tolerance'] = float(tolerance)
		if show_holes is not None:
			slicer['show_holes'] = bool(show_holes)
		if enabled is not None:
			# Like the prototype, switching the slicer on hides the drill holes, which would stand
			# in front of the cut, and switching it off shows them again
			if bool(enabled) != slicer['enabled']:
				for layer in self._layers:
					if layer.kind == 'holes':
						layer.visible = not enabled
						if layer.actor is not None:
							layer.actor.visibility = layer.visible
			slicer['enabled'] = bool(enabled)
		if slicer['enabled'] and slicer['origin'] is None:
			lower, upper = _corners(self._extents)
			slicer['origin'] = tuple(float(v) for v in (lower + upper) / 2)
		# While the slicer is off (and was off), moving its plane changes nothing that is drawn
		if slicer['enabled'] or was_enabled:
			for layer in self._layers:
				if layer.kind == 'blocks' and layer.source is not None:
					self._place(layer)
				elif layer.kind == 'holes':
					self._place_crossing(layer)
			self._sync_scalar_bars()
		self._refresh_panel()

	def _place_crossing(self, layer):
		"""Draw the drill hole intervals of a layer that cross the plane, when the slicer is on and
		``show_holes`` is set. As in the prototype: only holes whose tube reaches the plane, and of
		those the intervals within the tolerance band, clipped to the band. They follow the layer's
		threshold and colours, and show even when the layer itself is hidden."""
		import numpy as np
		from .meshes import SEGMENT, HOLE, UPPER, LOWER, RADIUS
		slicer = self._slicer
		if not (slicer['enabled'] and slicer['show_holes']):
			self._clear_crossing(layer)
			return
		mesh = layer.source.mesh
		upper, lower = mesh.field_data[UPPER], mesh.field_data[LOWER]
		radius = float(mesh.field_data[RADIUS][0])
		origin = np.array(slicer['origin'], dtype=float)
		normal = np.array(slicer['normal'], dtype=float)
		normal = normal / np.linalg.norm(normal)
		segment, hole = mesh.cell_data[SEGMENT], mesh.cell_data[HOLE]
		hole_of_segment = np.empty(len(upper), dtype=int)
		hole_of_segment[segment] = hole
		touching = _in_band(upper, lower, origin, normal, 2 * radius)
		crossing = (_in_band(upper, lower, origin, normal, slicer['tolerance'])
					& np.isin(hole_of_segment, hole_of_segment[touching]))
		tubes = mesh.extract_cells(np.flatnonzero(np.isin(segment, np.flatnonzero(crossing))))
		if tubes.n_cells:
			tubes = _Source(tubes, layer.source.var, layer.source.griddef).cut(layer.threshold, layer.section)
		if tubes.n_cells:
			half = slicer['tolerance'] / 2
			tubes = tubes.clip(normal=normal, origin=origin + half * normal, invert=True)
			tubes = tubes.clip(normal=-normal, origin=origin - half * normal, invert=True)
		if not tubes.n_cells:
			self._clear_crossing(layer)
			return
		layer.crossing = tubes
		layer.crossing_actor = self.plotter.add_mesh(tubes, name=layer.key + ' on slice', nan_color='lightgray',
													 show_scalar_bar=False, **layer.kwargs)

	def _clear_crossing(self, layer):
		"""Take the holes on the slice off the scene"""
		if layer.crossing_actor is not None:
			self.plotter.remove_actor(layer.crossing_actor)
		layer.crossing, layer.crossing_actor = None, None

	def _sync_scalar_bars(self):
		"""A colour bar is shown while a layer with its title is shown or has holes on the slice"""
		titles = {}
		for layer in self._layers:
			if 'scalars' in layer.kwargs:
				title = layer.kwargs['scalar_bar_args']['title']
				titles[title] = titles.get(title, False) or (
					layer.actor is not None and layer.visible) or layer.crossing is not None
		for title, shown in titles.items():
			if title in self.plotter.scalar_bars.keys():
				self.plotter.scalar_bars[title].SetVisibility(shown)

	def _plane(self):
		"""The plane of the slicer, or None when it is off"""
		slicer = self._slicer
		if not slicer['enabled']:
			return None
		return {'origin': slicer['origin'], 'normal': slicer['normal']}

	def set_axis_origin(self, origin=None, from_zero=False):
		"""
		Count the labels of the ruler from a point instead of from the coordinates, e.g. to
		compare two models whose coordinates are far apart. The titles say where the labels start.

		Parameters:
			origin (tuple): (easting, northing, elevation) that the label 0 stands for. If None,
				the coordinates are shown as they are
			from_zero (bool): start every axis at 0 on the lowest edge of the scene, instead of
				``origin``
		"""
		if origin is not None and len(origin) != 3:
			raise ValueError("`origin` must be (easting, northing, elevation)")
		self._axis_origin = None if origin is None else tuple(float(v) for v in origin)
		self._axes_from_zero = bool(from_zero)
		if self._ruler is not None:
			self._draw_ruler(self.plotter)
		self._refresh_panel()

	def _origin_values(self):
		"""The point the ruler counts from, as three numbers"""
		import numpy as np
		if self._axes_from_zero and self._extents:
			return tuple(np.array(self._extents, dtype=float).reshape(-1, 3, 2)[:, :, 0].min(axis=0))
		return self._axis_origin or (0.0, 0.0, 0.0)

	def _draw_ruler(self, plotter):
		_draw_ruler(plotter, self._extents, self._ruler, self._vertical_exaggeration, self._origin_values())

	def _refresh_panel(self):
		if self._panel is not None:
			self._panel.refresh()

	def set_vertical_exaggeration(self, factor):
		"""
		Stretch the scene vertically, e.g. to see a thin deposit spread over a large area. The
		ruler keeps showing the real elevations.

		Parameters:
			factor (float): vertical exaggeration, 1 for none
		"""
		self._vertical_exaggeration = factor
		self.plotter.set_scale(zscale=factor)
		# set_scale also stretches the ruler, so it is drawn again
		if self._ruler is not None:
			self._draw_ruler(self.plotter)

	def panel(self, views_folder='views'):
		"""
		The scene with a panel of controls beside it: switch layers on and off, change their
		opacity, colour limits, threshold and slice, move the origin of the ruler, and save the
		view or load a saved one. Requires ``pip install pygeostat[view3d-panel]``. In a notebook,
		make the result the last line of a cell to show it.

		Parameters:
			views_folder (str): folder for the saved views, one json file each (see
				:meth:`get_view`). It is made when the first view is saved

		Returns:
			panel (pygeostat.plotting3d.panel.Panel): the panel
		"""
		from .panel import Panel
		self._panel = Panel(self, views_folder)
		return self._panel

	def show(self, **kwargs):
		"""Show the scene. Keyword arguments go to ``pyvista.Plotter.show``."""
		return self.plotter.show(**kwargs)

	def get_camera(self):
		"""
		The current view, as a dictionary of plain numbers that can be saved (e.g. with ``json``)
		and given back to :meth:`set_camera` to reproduce the view exactly.

		Returns:
			camera (dict): ``position``, ``focal_point`` and ``view_up`` (x, y, z lists),
			``parallel_projection`` (bool), and the zoom: ``view_angle`` (degrees, used by a
			perspective view) and ``parallel_scale`` (used by a parallel view)
		"""
		camera = self.plotter.camera
		return {'position': [float(v) for v in camera.position],
				'focal_point': [float(v) for v in camera.focal_point],
				'view_up': [float(v) for v in camera.up],
				'parallel_projection': bool(camera.parallel_projection),
				'view_angle': float(camera.view_angle),
				'parallel_scale': float(camera.parallel_scale)}

	def set_camera(self, camera):
		"""
		Set the view from a dictionary made by :meth:`get_camera`.

		Parameters:
			camera (dict): the view, as returned by :meth:`get_camera`
		"""
		_apply_camera(self.plotter, camera)

	def get_view(self):
		"""
		Everything about how the scene is shown, as a dictionary of plain numbers and strings that
		can be saved (e.g. with ``json``) and given back to :meth:`set_view` on a scene with the
		same layers: the camera, the vertical exaggeration, the origin of the ruler, the slicer and, for each
		layer, whether it is shown, its opacity, colour limits, threshold and slice.

		Returns:
			view (dict): ``camera`` (see :meth:`get_camera`), ``vertical_exaggeration``,
			``axis_origin`` (list or None), ``axes_from_zero`` (bool), ``slicer`` (``enabled``,
			``origin``, ``normal``, ``tolerance`` and ``show_holes``) and ``layers`` (a list with
			one dictionary per layer, with its ``name``)
		"""
		layers = []
		for layer in self._layers:
			layers.append({
				'name': layer.name, 'visible': layer.visible,
				'opacity': float(layer.kwargs.get('opacity', 1.0)),
				'clim': None if layer.clim is None else [float(v) for v in layer.clim],
				'threshold': None if layer.threshold is None else [float(v) for v in layer.threshold],
				'slice': None if layer.section is None else {
					k: None if v is None else v if k == 'orient' else float(v)
					for k, v in layer.section.items()}})
		slicer = dict(self._slicer)
		slicer['origin'] = slicer['origin'] and list(slicer['origin'])
		slicer['normal'] = list(slicer['normal'])
		return {'camera': self.get_camera(), 'slicer': slicer,
				'vertical_exaggeration': float(self._vertical_exaggeration),
				'axis_origin': self._axis_origin and list(self._axis_origin),
				'axes_from_zero': self._axes_from_zero, 'layers': layers}

	def set_view(self, view):
		"""
		Show the scene as saved by :meth:`get_view`. The layers are matched by name, and the scene
		may have more layers than the view, which stay as they are. A view that cannot be shown
		(it names a layer the scene does not have, or a threshold leaves nothing) is refused and the
		scene stays as it was.

		Parameters:
			view (dict): the view, as returned by :meth:`get_view`
		"""
		missing = [saved['name'] for saved in view['layers'] if saved['name'] not in self.layers]
		if missing:
			raise ValueError("The view has layers {} that the scene does not; the layers are {}".format(
				missing, self.layers))
		before = self.get_view()
		try:
			self._apply_view(view)
		except ValueError:
			# e.g. a threshold that leaves nothing in this data: the scene goes back as it was
			self._apply_view(before)
			raise

	def _apply_view(self, view):
		# The slicer first: switching it on hides the drill holes, then the layers are set as saved
		if view.get('slicer') is not None:
			self.set_slicer(**view['slicer'])
		for saved in view['layers']:
			cut = dict(saved['slice'] or {'orient': False})
			if cut.get('slice_number') is not None:
				cut['slice_number'] = int(cut['slice_number'])
			if cut.get('slice_thickness') is not None:
				cut['slice_thickness'] = int(cut['slice_thickness'])
			self.set_layer(saved['name'], visible=saved['visible'], opacity=saved['opacity'],
						   clim=saved['clim'], threshold=saved['threshold'] or False, **cut)
		self.set_axis_origin(view['axis_origin'], view['axes_from_zero'])
		self.set_vertical_exaggeration(view['vertical_exaggeration'])
		self.set_camera(view['camera'])

	def show_static(self, figsize=None, dpi=None, title=None, cbar=True, ax=None):
		"""
		Draw the scene as a picture on a matplotlib axis, so it can sit next to 2D plots in one
		figure. The scene is rendered off screen, so this works before or after ``show``.

		Parameters:
			figsize (tuple): figure size (width, height) in inches, if ``ax`` is None
			dpi (float): resolution of the picture. If None, the figure's own dpi
			title (str): title of the axis
			cbar (bool): draw a matplotlib colour bar for each layer coloured by a variable
			ax (mpl.axis): axis to draw on. If None, a new figure is made

		Returns:
			fig (mpl.figure.Figure): the figure
			ax (mpl.axis): the axis with the picture
		"""
		import pyvista as pv
		import matplotlib.pyplot as plt
		from matplotlib.cm import ScalarMappable
		from matplotlib.colors import Normalize

		if ax is None:
			fig, ax = plt.subplots(figsize=figsize)
		else:
			fig = ax.figure
		if dpi is None:
			dpi = fig.dpi

		# PyVista's colour bars overlap the picture at small sizes, so matplotlib draws them instead.
		# They go on first because they shrink the axis, and the picture is rendered at its size.
		shown = [layer for layer in self._layers if layer.visible and layer.actor is not None]
		# A layer with holes on the slice keeps its colour bar, though it is hidden
		cbars = {}
		for layer in shown + [l for l in self._layers if l.crossing is not None and l not in shown]:
			kwargs = layer.kwargs
			if 'scalars' in kwargs:
				label = kwargs['scalar_bar_args']['title']
				clim = kwargs['clim'] or layer.mesh.get_data_range(kwargs['scalars'])
				ticks = kwargs.get('annotations')
				cbars.setdefault(label, (kwargs['cmap'], clim, ticks))
		if cbar:
			for label, (cmap, clim, ticks) in cbars.items():
				bar = fig.colorbar(ScalarMappable(Normalize(*clim), cmap), ax=ax, label=label)
				if ticks is not None:
					bar.set_ticks(list(ticks), labels=list(ticks.values()))
		bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		window_size = [int(round(bbox.width * dpi)), int(round(bbox.height * dpi))]

		static = pv.Plotter(off_screen=True, window_size=window_size)
		for layer in shown:
			static.add_mesh(layer.mesh, show_scalar_bar=False, **layer.kwargs)
		for layer in self._layers:
			if layer.crossing is not None:
				static.add_mesh(layer.crossing, show_scalar_bar=False, nan_color='lightgray', **layer.kwargs)
		static.set_scale(zscale=self._vertical_exaggeration)
		if self._ruler is not None:
			self._draw_ruler(static)
		if self.plotter.camera_set:
			_apply_camera(static, self.get_camera())
		image = static.screenshot(return_img=True)
		static.close()

		ax.imshow(image)
		ax.axis('off')
		if title is not None:
			ax.set_title(title)
		return fig, ax


def _in_band(upper, lower, origin, normal, thickness):
	"""Whether each segment (from ``upper`` to ``lower``, one row each) has a point within half of
	``thickness`` of the plane. Distance to a plane is linear along a segment, so its two ends
	give the nearest and the farthest point. With thickness 0 this is a crossing of the plane."""
	import numpy as np
	at_upper, at_lower = (upper - origin) @ normal, (lower - origin) @ normal
	half = thickness / 2
	return (np.minimum(at_upper, at_lower) <= half) & (np.maximum(at_upper, at_lower) >= -half)


def _corners(extents):
	"""The lowest and the highest corner (x, y, z) of a list of bounds"""
	import numpy as np
	bounds = np.array(extents, dtype=float).reshape(-1, 3, 2)
	return bounds[:, :, 0].min(axis=0), bounds[:, :, 1].max(axis=0)


def _draw_ruler(plotter, extents, titles, vertical_exaggeration=1.0, origin=(0.0, 0.0, 0.0)):
	"""Draw the ruler around the layers' extents. They come from the layers, not the plotter,
	which also holds the ruler itself, so redrawing it never grows the box. With vertical
	exaggeration the box is stretched like the scene but labelled with the real elevations. The
	labels count from ``origin``: the box is rounded in those counts, then moved back."""
	import numpy as np
	bounds = np.array(extents, dtype=float).reshape(-1, 3, 2)
	origin = np.array(origin, dtype=float)
	lower, upper = bounds[:, :, 0].min(axis=0) - origin, bounds[:, :, 1].max(axis=0) - origin
	# VTK spaces the labels evenly over the box, so for round labels the box ends on multiples of
	# a round step (1, 2 or 5 times a power of 10). The longest axis on screen gets about 5 labels
	# and a shorter one proportionally fewer, so labels on a short axis don't run together
	spans = np.where(upper > lower, upper - lower, 1.0)
	on_screen = spans * [1, 1, vertical_exaggeration]
	box, nlabels, decimals = [], [], []
	for lo, hi, span, shown in zip(lower, upper, spans, on_screen):
		intervals = max(2, round(4 * shown / on_screen.max()))
		power = 10 ** np.floor(np.log10(span / intervals))
		step = next(m * power for m in (1, 2, 5, 10) if m * power >= span / intervals)
		lo, hi = np.floor(lo / step) * step, np.ceil(hi / step) * step
		if hi == lo:
			hi = lo + step
		box += [lo, hi]
		nlabels.append(int(round((hi - lo) / step)) + 1)
		decimals.append(max(0, int(-np.floor(np.log10(step)))))
	placed = [v + origin[i // 2] for i, v in enumerate(box)]
	stretched = placed[:4] + [placed[4] * vertical_exaggeration, placed[5] * vertical_exaggeration]
	titles = {key: title + (' from {:g}'.format(origin[i]) if origin[i] else '')
			  for i, (key, title) in enumerate(titles.items())}
	axes = plotter.show_bounds(bounds=stretched, axes_ranges=box, grid='back', location='outer', ticks='both',
							   minor_ticks=False, font_size=9, n_xlabels=nlabels[0],
							   n_ylabels=nlabels[1], n_zlabels=nlabels[2], color='black', **titles)
	axes.x_label_format = '%.{}f'.format(decimals[0])
	axes.y_label_format = '%.{}f'.format(decimals[1])
	axes.z_label_format = '%.{}f'.format(decimals[2])


# Cell array with each cell's category as 0, 1, 2, ... in code order, which is what gets coloured
_CATEGORY = 'category index'


def _categories(mesh, var, cmap, catdict):
	"""Colour setup for a categorical variable, following location_plot: the categories present,
	in code order, each with the colour at its place in catdict (or in code order without one).
	Returns a copy of the mesh with the category index added, the colormap, and the labels."""
	import numpy as np
	from matplotlib.colors import ListedColormap
	from ..plotting.utils import _get_cmap
	values = np.asarray(mesh.cell_data[var], dtype=float)
	known = np.isfinite(values)
	codes = np.unique(values[known]).astype(int)
	if isinstance(catdict, dict):
		missing = sorted(int(code) for code in set(codes) - set(catdict))
		if missing:
			raise ValueError("Codes {} of '{}' are not in catdict".format(missing, var))
		colours = _get_cmap(cmap, True, len(catdict)).colors
		colours = [colours[list(catdict).index(code)] for code in codes]
		names = [catdict[code] for code in codes]
	else:
		colours = _get_cmap(cmap, True, len(codes)).colors
		names = [str(code) for code in codes]
	index = np.full(len(values), np.nan)
	index[known] = np.searchsorted(codes, values[known].astype(int))
	mesh = mesh.copy(deep=False)
	mesh.cell_data[_CATEGORY] = index
	return mesh, ListedColormap(colours), dict(enumerate(names))


def _slice(mesh, griddef, orient, slice_number, slice_coordinate, slice_thickness, pointtol):
	"""The part of a mesh in a slice band. Blocks are kept whole (by their centre); drill holes,
	points and surfaces are cut at the band's edges."""
	import numpy as np
	import pyvista as pv
	axes = {'xy': 2, 'xz': 1, 'yz': 0}
	if orient not in axes:
		raise ValueError("`orient` must be 'xy', 'xz' or 'yz'")
	axis = axes[orient]
	if slice_coordinate is None:
		if slice_number is None or griddef is None:
			raise ValueError("Pass `slice_coordinate`, or `slice_number` for data on a grid")
		slice_coordinate = griddef.get_slice_coordinate(orient, slice_number)
	if griddef is not None:
		size = (griddef.xsiz, griddef.ysiz, griddef.zsiz)[axis]
		# Half a block past the outermost centre, so exactly 2 * slice_thickness + 1 layers
		tolerance = ((slice_thickness or 0) + 0.5) * size
	elif pointtol is None:
		raise ValueError("Pass `pointtol`, the distance either side of the slice to show")
	else:
		tolerance = pointtol
	if isinstance(mesh, pv.PolyData):
		band = np.array(mesh.bounds, dtype=float)
		band[2 * axis:2 * axis + 2] = slice_coordinate - tolerance, slice_coordinate + tolerance
		return mesh.clip_box(band.tolist(), invert=False)
	distance = np.abs(mesh.cell_centers().points[:, axis] - slice_coordinate)
	return mesh.extract_cells(np.flatnonzero(distance <= tolerance))


def _new_section(section, orient, slice_number, slice_coordinate, slice_thickness, pointtol):
	"""The slice arguments of a layer after set_layer: the current ones with the given ones
	changed, or None (the whole layer) for orient=False"""
	if orient is False:
		return None
	if section is None and orient is None:
		raise ValueError("Pass `orient` to slice a layer that is not sliced")
	given = dict(orient=orient, slice_number=slice_number, slice_coordinate=slice_coordinate,
				 slice_thickness=slice_thickness, pointtol=pointtol)
	section = {**dict.fromkeys(given), **(section or {}), **{k: v for k, v in given.items() if v is not None}}
	# The position is the number or the coordinate, whichever was given last
	if slice_number is not None:
		section['slice_coordinate'] = None
	if slice_coordinate is not None:
		section['slice_number'] = None
	return section


def _apply_camera(plotter, camera):
	"""Set a plotter's camera from a get_camera dictionary. The position alone misses the zoom and
	the projection."""
	plotter.camera.parallel_projection = camera['parallel_projection']
	plotter.camera_position = (camera['position'], camera['focal_point'], camera['view_up'])
	plotter.camera.view_angle = camera['view_angle']
	plotter.camera.parallel_scale = camera['parallel_scale']


def view3d(data, var=None, viewer=None, cmap=None, clim=None, opacity=1.0, color=None, label=None,
		   threshold=None, orient=None, slice_number=None, slice_coordinate=None,
		   slice_thickness=None, pointtol=None, catdict=None, **kwargs):
	"""
	Add drill holes, a block model or a mesh to a 3D scene. Calls share one scene through
	``viewer``, the same way 2D plots share ``ax``.

	Drill holes and block models are converted with
	:func:`gs.drillhole_mesh() <pygeostat.plotting3d.meshes.drillhole_mesh>` and
	:func:`gs.block_mesh() <pygeostat.plotting3d.meshes.block_mesh>`. ``data`` is a block model if
	``griddef`` or ``size`` is passed, or if it is a DataFile with a griddef; otherwise it is drill
	hole data.

	Requires PyVista (``pip install pyvista``).

	Parameters:
		data (pd.DataFrame, gs.DataFile or pyvista.DataSet): data to show, or a mesh such as a
			surface
		var (str): variable to colour by. If None and data is a DataFile with a single variable,
			that variable. If there is none, the layer is drawn in ``color``
		viewer (gs.Viewer): scene to add to. If None, a new one is made
		cmap (str or mpl.colors.Colormap): colormap. If None, based on
			gs.Parameters['plotting.cmap']
		clim (tuple): colour limits (min, max). If None, the range of ``var``
		opacity (float): 0 (invisible) to 1 (opaque)
		color: colour of a layer that is not coloured by a variable
		label (str): title of the colour bar and name of the layer (see ``Viewer.layers``). If None,
			``var``
		catdict (dict): names of the categories of ``var``, e.g. ``{1: 'Granite', 2: 'Diorite'}``.
			``var`` is coloured by category if this is passed, or if data is a DataFile and
			``var`` is its ``cat`` column (then DataFile.catdict is used). Categories get the
			same colours as in :func:`gs.location_plot() <pygeostat.plotting.location_plot.location_plot>`
		threshold (tuple): (min, max) of ``var`` to show; blocks, intervals or points outside it,
			or without a value, are hidden. Colours keep the full range of ``var`` unless ``clim``
			is passed
		orient (str): show only a slice: ``'xy'`` (plan), ``'xz'`` or ``'yz'`` (sections), as in
			:func:`gs.slice_plot() <pygeostat.plotting.slice_plot.slice_plot>`
		slice_number (int): grid index of the slice along the axis not in ``orient``, for data
			on a grid
		slice_coordinate (float): position of the slice along the axis not in ``orient``, instead
			of ``slice_number`` (e.g. for drill holes or sub-blocks)
		slice_thickness (int): for data on a grid, the number of extra grid slices shown on each
			side of the slice. Default 0, a single layer of blocks
		pointtol (float): for data that is not on a grid, the distance either side of the slice
			to show; drill holes are cut at that distance. Like ``threshold``, a slice keeps the
			colours of the full range of ``var`` unless ``clim`` is passed
		**kwargs: passed to ``drillhole_mesh`` or ``block_mesh`` (e.g. ``start``, ``end``,
			``trace``, ``radius``, ``griddef``, ``size``)

	Returns:
		viewer (gs.Viewer): the scene, to add more layers to or to show

	Examples:

	Drill holes on their own, then a block model and points in one scene:

	.. code-block:: python

		import pygeostat as gs
		data = gs.ExampleData('oilsands')
		gs.view3d(data, 'Bitumen', trace=True).show()

		blocks = gs.ExampleData('3d_grid')
		viewer = gs.view3d(blocks, 'value', griddef=gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1]),
						   opacity=0.3)
		gs.view3d(gs.ExampleData('3d_estimate'), 'Estimate', viewer=viewer)
		viewer.show()

	"""
	from ..data.data import DataFile
	from .meshes import drillhole_mesh, block_mesh, SEGMENT
	from ..plotting.utils import _get_cmap
	try:
		import pyvista as pv
	except ImportError:
		raise ImportError("view3d requires PyVista: pip install pyvista") from None

	if var is None and isinstance(data, DataFile) and isinstance(data.variables, str):
		var = data.variables
	griddef = None
	kind = 'blocks'
	if isinstance(data, pv.DataSet):
		mesh = data
	elif 'griddef' in kwargs or 'size' in kwargs or getattr(data, 'griddef', None) is not None:
		mesh = block_mesh(data, var=var, **kwargs)
		if 'size' not in kwargs:
			griddef = kwargs['griddef'] if 'griddef' in kwargs else data.griddef
	else:
		mesh = drillhole_mesh(data, var=var, **kwargs)
		kind = 'holes' if SEGMENT in mesh.cell_data else 'points'
	categorical = var is not None and (catdict is not None or (
		isinstance(data, DataFile) and var == data.cat))
	if categorical:
		if catdict is None:
			catdict = data.catdict
		mesh, cmap, annotations = _categories(mesh, var, cmap, catdict)
	# When only part of the layer is shown, the ruler keeps wrapping the whole layer and the
	# colours keep its full range, so moving a slice or cutoff doesn't change them
	extent = None
	data_range = None
	if threshold is not None or orient is not None:
		extent = mesh.bounds
		if var is not None:
			data_range = mesh.get_data_range(var)
			if clim is None:
				clim = data_range
	source = _Source(mesh, var, griddef)
	section = None if orient is None else dict(
		orient=orient, slice_number=slice_number, slice_coordinate=slice_coordinate,
		slice_thickness=slice_thickness, pointtol=pointtol)
	mesh = source.cut(threshold, section)

	if viewer is None:
		viewer = Viewer()
	name = label or var
	if var is None:
		viewer.add(mesh, extent=extent, name=name, source=source, threshold=threshold, section=section, kind=kind,
				   color=color, opacity=opacity)
	elif categorical:
		viewer.add(mesh, extent=extent, name=name, source=source, threshold=threshold, section=section, kind=kind,
				   scalars=_CATEGORY, cmap=cmap, opacity=opacity,
				   clim=(-0.5, len(annotations) - 0.5), n_colors=len(annotations),
				   annotations=annotations, scalar_bar_args={'title': label or var})
	else:
		viewer.add(mesh, extent=extent, name=name, source=source, threshold=threshold, section=section, kind=kind,
				   data_range=data_range, scalars=var, cmap=_get_cmap(cmap, False, None), clim=clim,
				   opacity=opacity, scalar_bar_args={'title': label or var})
	return viewer
