# !/usr/bin/env python
#  -*- coding: utf-8 -*-
"""view3d shows drill holes, block models and surfaces together in an interactive 3D scene"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


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
		# Each layer's mesh and drawing options, so show_static can rebuild the scene
		self._layers = []

	def add(self, mesh, **kwargs):
		"""Add a mesh to the scene. Keyword arguments go to ``pyvista.Plotter.add_mesh``."""
		self._layers.append((mesh, kwargs))
		self.plotter.add_mesh(mesh, **kwargs)

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
		cbars = {}
		for mesh, kwargs in self._layers:
			if 'scalars' in kwargs:
				label = kwargs['scalar_bar_args']['title']
				clim = kwargs['clim'] or mesh.get_data_range(kwargs['scalars'])
				cbars.setdefault(label, (kwargs['cmap'], clim))
		if cbar:
			for label, (cmap, clim) in cbars.items():
				fig.colorbar(ScalarMappable(Normalize(*clim), cmap), ax=ax, label=label)
		bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		window_size = [int(round(bbox.width * dpi)), int(round(bbox.height * dpi))]

		static = pv.Plotter(off_screen=True, window_size=window_size)
		for mesh, kwargs in self._layers:
			static.add_mesh(mesh, show_scalar_bar=False, **kwargs)
		if self.plotter.camera_set:
			_apply_camera(static, self.get_camera())
		image = static.screenshot(return_img=True)
		static.close()

		ax.imshow(image)
		ax.axis('off')
		if title is not None:
			ax.set_title(title)
		return fig, ax


def _apply_camera(plotter, camera):
	"""Set a plotter's camera from a get_camera dictionary. The position alone misses the zoom and
	the projection."""
	plotter.camera.parallel_projection = camera['parallel_projection']
	plotter.camera_position = (camera['position'], camera['focal_point'], camera['view_up'])
	plotter.camera.view_angle = camera['view_angle']
	plotter.camera.parallel_scale = camera['parallel_scale']


def view3d(data, var=None, viewer=None, cmap=None, clim=None, opacity=1.0, color=None, label=None,
		   **kwargs):
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
		label (str): title of the colour bar. If None, ``var``
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
	if isinstance(data, pv.DataSet):
		mesh = data
	elif 'griddef' in kwargs or 'size' in kwargs or getattr(data, 'griddef', None) is not None:
		mesh = block_mesh(data, var=var, **kwargs)
	else:
		mesh = drillhole_mesh(data, var=var, **kwargs)

	if viewer is None:
		viewer = Viewer()
	if var is None:
		viewer.add(mesh, color=color, opacity=opacity)
	else:
		viewer.add(mesh, scalars=var, cmap=_get_cmap(cmap, False, None), clim=clim, opacity=opacity,
				   scalar_bar_args={'title': label or var})
	return viewer
