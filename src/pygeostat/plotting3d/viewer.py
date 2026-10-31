# !/usr/bin/env python
#  -*- coding: utf-8 -*-
"""view3d shows drill holes, block models and surfaces together in an interactive 3D scene"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


class _Layer:
	"""A mesh in the scene with the options it was drawn with, which the panel changes later"""

	def __init__(self, name, mesh, kwargs, actor):
		self.name = name
		self.mesh = mesh
		self.kwargs = kwargs
		self.actor = actor
		self.visible = True


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
		self._vertical_exaggeration = 1.0
		# The panel made by panel(), which shows the layers
		self._panel = None

	@property
	def layers(self):
		"""Names of the layers, in the order they were added (list of str)"""
		return [layer.name for layer in self._layers]

	def add(self, mesh, extent=None, name=None, **kwargs):
		"""
		Add a mesh to the scene. Keyword arguments go to ``pyvista.Plotter.add_mesh``.

		Parameters:
			mesh (pyvista.DataSet): the mesh
			extent (tuple): bounds (xmin, xmax, ymin, ymax, zmin, zmax) for the ruler to wrap,
				e.g. the whole model when ``mesh`` is part of it. If None, the bounds of ``mesh``
			name (str): name of the layer, to change it with :meth:`set_layer`. A name that is
				already taken gets a number added. If None, ``'Layer 1'``, ``'Layer 2'``, ...
		"""
		name = name or 'Layer {}'.format(len(self._layers) + 1)
		first, number = name, 1
		while name in self.layers:
			number += 1
			name = '{} ({})'.format(first, number)
		self._layers.append(_Layer(name, mesh, kwargs, self.plotter.add_mesh(mesh, **kwargs)))
		self._extents.append(mesh.bounds if extent is None else extent)
		if self._vertical_exaggeration != 1:
			self.set_vertical_exaggeration(self._vertical_exaggeration)
		elif self._ruler is not None:
			_draw_ruler(self.plotter, self._extents, self._ruler)
		if self._panel is not None:
			self._panel.refresh()

	def set_layer(self, name, visible=None, opacity=None):
		"""
		Change a layer that is already in the scene. The change shows in :meth:`show` and in
		:meth:`show_static`. The ruler keeps wrapping a hidden layer, so hiding one does not move
		the axes.

		Parameters:
			name (str): name of the layer, one of :attr:`layers`
			visible (bool): show or hide the layer. If None, unchanged
			opacity (float): 0 (invisible) to 1 (opaque). If None, unchanged
		"""
		matches = [layer for layer in self._layers if layer.name == name]
		if not matches:
			raise ValueError("No layer '{}'; the layers are {}".format(name, self.layers))
		layer = matches[0]
		if opacity is not None:
			if not 0 <= opacity <= 1:
				raise ValueError("`opacity` must be between 0 and 1")
			layer.kwargs['opacity'] = opacity
			layer.actor.prop.opacity = opacity
		if visible is not None:
			layer.visible = bool(visible)
			layer.actor.visibility = layer.visible
		if self._panel is not None:
			self._panel.refresh()

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
		_draw_ruler(self.plotter, self._extents, self._ruler, self._vertical_exaggeration)

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
			_draw_ruler(self.plotter, self._extents, self._ruler, factor)

	def panel(self):
		"""
		The scene with a panel of controls beside it, to switch layers on and off and change their
		opacity. Requires ``pip install pygeostat[view3d-panel]``. In a notebook, make the result the
		last line of a cell to show it.

		Returns:
			panel (pygeostat.plotting3d.panel.Panel): the panel
		"""
		from .panel import Panel
		self._panel = Panel(self)
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
		shown = [layer for layer in self._layers if layer.visible]
		cbars = {}
		for layer in shown:
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
		static.set_scale(zscale=self._vertical_exaggeration)
		if self._ruler is not None:
			_draw_ruler(static, self._extents, self._ruler, self._vertical_exaggeration)
		if self.plotter.camera_set:
			_apply_camera(static, self.get_camera())
		image = static.screenshot(return_img=True)
		static.close()

		ax.imshow(image)
		ax.axis('off')
		if title is not None:
			ax.set_title(title)
		return fig, ax


def _draw_ruler(plotter, extents, titles, vertical_exaggeration=1.0):
	"""Draw the ruler around the layers' extents. They come from the layers, not the plotter,
	which also holds the ruler itself, so redrawing it never grows the box. With vertical
	exaggeration the box is stretched like the scene but labelled with the real elevations."""
	import numpy as np
	bounds = np.array(extents, dtype=float).reshape(-1, 3, 2)
	lower, upper = bounds[:, :, 0].min(axis=0), bounds[:, :, 1].max(axis=0)
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
	stretched = box[:4] + [box[4] * vertical_exaggeration, box[5] * vertical_exaggeration]
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
	from .meshes import drillhole_mesh, block_mesh
	from ..plotting.utils import _get_cmap
	try:
		import pyvista as pv
	except ImportError:
		raise ImportError("view3d requires PyVista: pip install pyvista") from None

	if var is None and isinstance(data, DataFile) and isinstance(data.variables, str):
		var = data.variables
	griddef = None
	if isinstance(data, pv.DataSet):
		mesh = data
	elif 'griddef' in kwargs or 'size' in kwargs or getattr(data, 'griddef', None) is not None:
		mesh = block_mesh(data, var=var, **kwargs)
		if 'size' not in kwargs:
			griddef = kwargs['griddef'] if 'griddef' in kwargs else data.griddef
	else:
		mesh = drillhole_mesh(data, var=var, **kwargs)
	categorical = var is not None and (catdict is not None or (
		isinstance(data, DataFile) and var == data.cat))
	if categorical:
		if catdict is None:
			catdict = data.catdict
		mesh, cmap, annotations = _categories(mesh, var, cmap, catdict)
	# When only part of the layer is shown, the ruler keeps wrapping the whole layer and the
	# colours keep its full range, so moving a slice or cutoff doesn't change them
	extent = None
	if threshold is not None or orient is not None:
		extent = mesh.bounds
		if clim is None and var is not None:
			clim = mesh.get_data_range(var)
	if threshold is not None:
		if var is None:
			raise ValueError("`threshold` needs a variable `var`")
		mesh = mesh.threshold(threshold, scalars=var)
	if orient is not None:
		mesh = _slice(mesh, griddef, orient, slice_number, slice_coordinate, slice_thickness,
					  pointtol)

	if viewer is None:
		viewer = Viewer()
	name = label or var
	if var is None:
		viewer.add(mesh, extent=extent, name=name, color=color, opacity=opacity)
	elif categorical:
		viewer.add(mesh, extent=extent, name=name, scalars=_CATEGORY, cmap=cmap, opacity=opacity,
				   clim=(-0.5, len(annotations) - 0.5), n_colors=len(annotations),
				   annotations=annotations, scalar_bar_args={'title': label or var})
	else:
		viewer.add(mesh, extent=extent, name=name, scalars=var, cmap=_get_cmap(cmap, False, None), clim=clim,
				   opacity=opacity, scalar_bar_args={'title': label or var})
	return viewer
