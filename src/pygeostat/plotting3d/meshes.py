# !/usr/bin/env python
#  -*- coding: utf-8 -*-
"""meshes converts pygeostat data into PyVista meshes for 3D viewing"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


def drillhole_mesh(data, var=None, x=None, y=None, z=None, start=None, end=None, radius=2.0,
				   n_sides=20):
	"""
	Convert drill hole data to a PyVista mesh. Each interval with ``start`` and ``end``
	coordinate columns becomes a cylinder from its start to its end. Without them, each row
	becomes a point at ``x``, ``y``, ``z``. The values of ``var`` are stored as cell data, one
	value per interval or point.

	Requires PyVista (``pip install pyvista``).

	Parameters:
		data (pd.DataFrame or gs.DataFile): drill hole data
		var (str or list): variable(s) to store on the mesh. If None and data is a DataFile,
			based on DataFile.variables
		x (str): x coordinate column. If None and data is a DataFile, based on DataFile.x
		y (str): y coordinate column. If None and data is a DataFile, based on DataFile.y
		z (str): z coordinate column. If None and data is a DataFile, based on DataFile.z
		start (tuple): x, y, z columns of the start of each interval
		end (tuple): x, y, z columns of the end of each interval
		radius (float): cylinder radius, in the units of the coordinates
		n_sides (int): number of sides of each cylinder

	Returns:
		mesh (pyvista.PolyData): cylinders, or points if ``start`` and ``end`` are None

	Examples:

	Points from a DataFile that knows its coordinate columns:

	.. code-block:: python

		import pygeostat as gs
		data = gs.ExampleData('oilsands')
		mesh = gs.drillhole_mesh(data, var='Bitumen')
		mesh.plot()

	"""
	import numpy as np
	from ..data.data import DataFile
	try:
		import pyvista as pv
	except ImportError:
		raise ImportError("drillhole_mesh requires PyVista: pip install pyvista") from None

	if isinstance(data, DataFile):
		if var is None:
			var = data.variables
		x = data.x if x is None else x
		y = data.y if y is None else y
		z = data.z if z is None else z
		data = data.data
	if var is None:
		var = []
	elif isinstance(var, str):
		var = [var]
	if (start is None) != (end is None):
		raise ValueError("Pass both `start` and `end`, or neither")

	if start is None:
		if any(col is None for col in (x, y, z)):
			raise ValueError("The `x`, `y` and `z` columns have not been set and could not be"
							 " automatically retrieved")
		mesh = pv.PolyData(data[[x, y, z]].to_numpy(dtype=float))
	else:
		n = len(data)
		points = np.empty((2 * n, 3))
		points[0::2] = data[list(start)].to_numpy(dtype=float)
		points[1::2] = data[list(end)].to_numpy(dtype=float)
		lines = np.column_stack([np.full(n, 2), np.arange(0, 2 * n, 2),
								 np.arange(1, 2 * n, 2)]).ravel()
		mesh = pv.PolyData(points, lines=lines)
	for v in var:
		mesh.cell_data[v] = data[v].to_numpy(dtype=float)
	if start is not None:
		mesh = mesh.tube(radius=radius, n_sides=n_sides)
	return mesh
