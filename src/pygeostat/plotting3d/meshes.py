# !/usr/bin/env python
#  -*- coding: utf-8 -*-
"""meshes converts pygeostat data into PyVista meshes for 3D viewing"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


def drillhole_mesh(data, var=None, x=None, y=None, z=None, start=None, end=None, dh=None,
				   trace=False, max_gap=None, radius=2.0, n_sides=20):
	"""
	Convert drill hole data to a PyVista mesh. Each interval with ``start`` and ``end``
	coordinate columns becomes a cylinder from its start to its end. With ``trace=True``, the
	samples of each drill hole are joined from top to bottom by cylinders, each carrying the value
	of its upper sample. Otherwise each row becomes a point at ``x``, ``y``, ``z``. The values of
	``var`` are stored as cell data, one value per cylinder or point.

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
		dh (str): drill hole id column, used by ``trace``. If None and data is a DataFile, based
			on DataFile.dh
		trace (bool): join the samples of each drill hole, for data without ``start`` and ``end``
		max_gap (float): with ``trace``, consecutive samples further apart than this are not
			joined, so no cylinder is drawn across a stretch with no samples. If None, three times
			the median distance between consecutive samples
		radius (float): cylinder radius, in the units of the coordinates
		n_sides (int): number of sides of each cylinder

	Returns:
		mesh (pyvista.PolyData): cylinders, or points without ``start``, ``end`` or ``trace``

	Examples:

	Points from a DataFile that knows its coordinate columns:

	.. code-block:: python

		import pygeostat as gs
		data = gs.ExampleData('oilsands')
		mesh = gs.drillhole_mesh(data, var='Bitumen')
		mesh.plot()

	"""
	import numpy as np
	import pandas as pd
	from ..data.data import DataFile
	from ..utility.logging import printerr
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
		dh = data.dh if dh is None else dh
		data = data.data
	if var is None:
		var = []
	elif isinstance(var, str):
		var = [var]
	if (start is None) != (end is None):
		raise ValueError("Pass both `start` and `end`, or neither")
	for v in var:
		if not pd.api.types.is_numeric_dtype(data[v]):
			raise ValueError("The variable '{}' is not numeric. Store categories as integer codes"
							 " (see DataFile.catdict)".format(v))

	if start is None and any(col is None for col in (x, y, z)):
		raise ValueError("The `x`, `y` and `z` columns have not been set and could not be"
						 " automatically retrieved")
	if start is None and not trace:
		mesh = pv.PolyData(data[[x, y, z]].to_numpy(dtype=float))
		for v in var:
			mesh.cell_data[v] = data[v].to_numpy(dtype=float)
		return mesh

	if start is not None:
		rows = np.arange(len(data))
		upper = data[list(start)].to_numpy(dtype=float)
		lower = data[list(end)].to_numpy(dtype=float)
	else:
		if dh is None:
			raise ValueError("`trace` needs the drill hole id column `dh`")
		xyz = data[[x, y, z]].to_numpy(dtype=float)
		holes = pd.factorize(data[dh])[0]
		order = np.lexsort((-xyz[:, 2], holes))
		same_hole = holes[order][1:] == holes[order][:-1]
		rows, below = order[:-1][same_hole], order[1:][same_hole]
		upper, lower = xyz[rows], xyz[below]
		if max_gap is None and len(rows):
			max_gap = 3 * np.median(np.linalg.norm(upper - lower, axis=1))
		keep = np.linalg.norm(upper - lower, axis=1) <= max_gap
		rows, upper, lower = rows[keep], upper[keep], lower[keep]
	n = len(rows)
	points = np.empty((2 * n, 3))
	points[0::2] = upper
	points[1::2] = lower
	lines = np.column_stack([np.full(n, 2), np.arange(0, 2 * n, 2),
							 np.arange(1, 2 * n, 2)]).ravel()
	mesh = pv.PolyData(points, lines=lines)
	for v in var:
		mesh.cell_data[v] = data[v].to_numpy(dtype=float)[rows]
	# The tube filter drops cylinders with no length, so they would vanish silently
	nzero = int((np.linalg.norm(upper - lower, axis=1) == 0).sum())
	if nzero:
		printerr("{} interval(s) have the same start and end and are not drawn".format(nzero),
				 errtype='warning')
	return mesh.tube(radius=radius, n_sides=n_sides)
