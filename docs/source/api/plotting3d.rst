.. _plotting3d:

3D Plotting Functions
#####################

Functions to view drill holes, block models and sub-blocks in 3D. They are built on
`PyVista <https://docs.pyvista.org/>`_, an optional dependency installed with
``pip install pygeostat[view3d]``. See the :ref:`3D viewer user guide <view3d>` for examples.

View3D
******

.. autofunction:: pygeostat.plotting3d.view3d

Viewer
******

.. autoclass:: pygeostat.plotting3d.Viewer
   :members: layers, add, set_layer, add_ruler, set_vertical_exaggeration, get_camera, set_camera, show,
              show_static

Meshes
******

These convert pygeostat data to PyVista meshes. :func:`view3d() <pygeostat.plotting3d.view3d>` calls them,
and they can be used alone to work with the mesh in PyVista.

.. autofunction:: pygeostat.plotting3d.drillhole_mesh

.. autofunction:: pygeostat.plotting3d.block_mesh
