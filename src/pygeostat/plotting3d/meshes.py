# !/usr/bin/env python
#  -*- coding: utf-8 -*-
"""meshes converts pygeostat data into PyVista meshes for 3D viewing"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


# Arrays that a drill hole mesh carries for slicing, named so they do not clash with a variable
SEGMENT = 'segment index'
HOLE = 'hole index'
UPPER = 'segment upper end'
LOWER = 'segment lower end'
RADIUS = 'tube radius'
SIDES = 'tube sides'
HOLE_NAMES = 'hole names'


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
		# A sample with no location can't be joined, and would make the median spacing NaN
		located = ~np.isnan(xyz).any(axis=1)
		if not located.all():
			printerr("{} sample(s) with missing coordinates are not drawn".format(
				int((~located).sum())), errtype='warning')
		order = np.flatnonzero(located)[np.lexsort((-xyz[located, 2], holes[located]))]
		same_hole = holes[order][1:] == holes[order][:-1]
		rows, below = order[:-1][same_hole], order[1:][same_hole]
		upper, lower = xyz[rows], xyz[below]
		if max_gap is None and len(rows):
			max_gap = 3 * np.median(np.linalg.norm(upper - lower, axis=1))
		keep = np.linalg.norm(upper - lower, axis=1) <= max_gap
		rows, upper, lower = rows[keep], upper[keep], lower[keep]
	n = len(rows)
	mesh = _segments(upper, lower)
	for v in var:
		mesh.cell_data[v] = data[v].to_numpy(dtype=float)[rows]
	# What a plane through the holes needs to find the intervals it crosses: each cell's interval
	# and drill hole, and the intervals' ends and the radius
	mesh.cell_data[SEGMENT] = np.arange(n)
	if dh is None:
		mesh.cell_data[HOLE] = np.arange(n)
	else:
		codes, holes = pd.factorize(data[dh])
		mesh.cell_data[HOLE] = codes[rows]
		# The name of each drill hole, by its index, for telling which hole an interval is of
		mesh.field_data[HOLE_NAMES] = np.array([
			'{:.0f}'.format(h) if isinstance(h, (float, np.floating)) and float(h).is_integer() else str(h)
			for h in holes])
	mesh.field_data[UPPER] = upper
	mesh.field_data[LOWER] = lower
	mesh.field_data[RADIUS] = [radius]
	mesh.field_data[SIDES] = [n_sides]
	# The tube filter drops cylinders with no length, so they would vanish silently
	nzero = int((np.linalg.norm(upper - lower, axis=1) == 0).sum())
	if nzero:
		printerr("{} interval(s) have the same start and end and are not drawn".format(nzero),
				 errtype='warning')
	return mesh.tube(radius=radius, n_sides=n_sides)


def _segments(upper, lower):
	"""A line from each row of ``upper`` to the same row of ``lower``, as one cell each"""
	import numpy as np
	import pyvista as pv
	n = len(upper)
	points = np.empty((2 * n, 3))
	points[0::2] = upper
	points[1::2] = lower
	lines = np.column_stack([np.full(n, 2), np.arange(0, 2 * n, 2), np.arange(1, 2 * n, 2)]).ravel()
	return pv.PolyData(points, lines=lines)


def _with_radius(tubes, radius):
	"""The tubes of :func:`drillhole_mesh` drawn again with another radius. The cells are the same
	and in the same order, so every array of the mesh is carried across, including the ones added
	to it since"""
	new = _segments(tubes.field_data[UPPER], tubes.field_data[LOWER]).tube(
		radius=radius, n_sides=int(tubes.field_data[SIDES][0]))
	if new.n_cells != tubes.n_cells:
		raise RuntimeError("The tubes of radius {:g} do not have the cells of the others".format(radius))
	for name in tubes.cell_data.keys():
		new.cell_data[name] = tubes.cell_data[name]
	for name in tubes.field_data.keys():
		new.field_data[name] = tubes.field_data[name]
	new.field_data[RADIUS] = [radius]
	return new


def block_mesh(data, var=None, griddef=None, x=None, y=None, z=None, size=None):
	"""
	Convert a block model to a PyVista mesh, one cell per block, with the values of ``var`` stored
	as cell data. A regular grid needs ``griddef``. Sub-blocks of any size need ``size``: each row
	is a box centred at ``x``, ``y``, ``z`` with the sizes in the three ``size`` columns.

	Requires PyVista (``pip install pyvista``).

	Parameters:
		data (pd.DataFrame or gs.DataFile): block model values, in grid order (x fastest, then y,
			then z)
		var (str or list): variable(s) to store on the mesh. If None and data is a DataFile,
			based on DataFile.variables
		griddef (gs.GridDef): grid definition. If None and data is a DataFile, based on
			DataFile.griddef
		x (str): x column of the sub-block centres. If None and data is a DataFile, based on
			DataFile.x
		y (str): y column of the sub-block centres. If None and data is a DataFile, based on
			DataFile.y
		z (str): z column of the sub-block centres. If None and data is a DataFile, based on
			DataFile.z
		size (tuple): x, y, z size columns of the sub-blocks

	Returns:
		mesh (pyvista.ImageData or pyvista.UnstructuredGrid): one cell per block

	Examples:

	.. code-block:: python

		import pygeostat as gs
		data = gs.ExampleData('3d_grid')
		griddef = gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1])
		mesh = gs.block_mesh(data, var='value', griddef=griddef)
		mesh.plot()

	Sub-blocks, one row per block with its centre and its sizes:

	.. code-block:: python

		mesh = gs.block_mesh(subblocks, var='Au', x='x', y='y', z='z',
							 size=('xsize', 'ysize', 'zsize'))

	"""
	import numpy as np
	from ..data.data import DataFile
	try:
		import pyvista as pv
	except ImportError:
		raise ImportError("block_mesh requires PyVista: pip install pyvista") from None

	if griddef is not None and size is not None:
		raise ValueError("Pass `griddef` for a regular grid or `size` for sub-blocks, not both")
	if isinstance(data, DataFile):
		if var is None:
			var = data.variables
		if griddef is None and size is None:
			griddef = data.griddef
		x = data.x if x is None else x
		y = data.y if y is None else y
		z = data.z if z is None else z
		data = data.data
	if var is None:
		var = []
	elif isinstance(var, str):
		var = [var]

	if size is not None:
		if any(col is None for col in (x, y, z)):
			raise ValueError("The `x`, `y` and `z` columns have not been set and could not be"
							 " automatically retrieved")
		n = len(data)
		centres = data[[x, y, z]].to_numpy(dtype=float)
		half = data[list(size)].to_numpy(dtype=float) / 2
		# Corners in VTK hexahedron order: bottom face counterclockwise, then top face
		corners = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
							[-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)
		points = (centres[:, None, :] + corners[None, :, :] * half[:, None, :]).reshape(-1, 3)
		cells = np.hstack([np.full((n, 1), 8), np.arange(8 * n).reshape(n, 8)]).ravel()
		mesh = pv.UnstructuredGrid(cells, np.full(n, pv.CellType.HEXAHEDRON), points)
		for v in var:
			mesh.cell_data[v] = data[v].to_numpy(dtype=float)
		return mesh

	if griddef is None:
		raise ValueError("`griddef` (regular grid) or `size` (sub-blocks) has not been set and"
						 " could not be automatically retrieved")
	if len(data) != griddef.count():
		raise ValueError("The data has {} rows but the grid has {} blocks".format(
			len(data), griddef.count()))
	mesh = pv.ImageData(dimensions=(griddef.nx + 1, griddef.ny + 1, griddef.nz + 1),
						spacing=(griddef.xsiz, griddef.ysiz, griddef.zsiz),
						origin=(griddef.xmn - griddef.xsiz / 2, griddef.ymn - griddef.ysiz / 2,
								griddef.zmn - griddef.zsiz / 2))
	for v in var:
		mesh.cell_data[v] = data[v].to_numpy(dtype=float)
	return mesh
