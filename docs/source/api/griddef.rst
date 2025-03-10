.. _griddef:

Grid Definitions
################



GridDef Class
=============
.. autoclass:: pygeostat.data.grid_definition.GridDef

Copy
++++
.. automethod:: pygeostat.data.grid_definition.GridDef.copy

Block Volume
++++++++++++
.. autoattribute:: pygeostat.data.grid_definition.GridDef.block_volume

Cell Count
++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.count


Convert To 2D
+++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.convert_to_2d

Convert Index: 3D to 1D
+++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.index3d_to_index1d

Convert Index: 1D to 3D
+++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.index1d_to_index3d

Coordinate(s) to 1D Index
+++++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.get_index

Coordinate(s) to 3D Index
+++++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.get_index3d

Change Block Size
+++++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.change_blocksize

Change Block Number
+++++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.change_blocknum

Pad Grid
++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.pad_grid

Extents
+++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.extents


Get Coordinates based on Index/Indices
++++++++++++++++++++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.get_coordinates

Get Slice Index
+++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.get_slice_index


Grid Array
++++++++++
.. autoattribute:: pygeostat.data.GridDef.grid_array

Origin
++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.origin

Outline Points
++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.outline_points


Sample Random Points from Grid
+++++++++++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.random_points

Sample Random Indices from Grid
++++++++++++++++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.random_indices


Slice Coordinate
++++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.get_slice_coordinate

Slice Index
+++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.get_slice_index

Vertical Indices
++++++++++++++++
.. automethod:: pygeostat.data.grid_definition.GridDef.get_vertical_indices


