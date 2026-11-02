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
   :members: layers, add, set_layer, set_slicer, add_ruler, set_axis_origin, set_vertical_exaggeration, get_camera,
             set_camera, get_view, set_view, panel, show, show_static

Panel
*****

The interactive panel made by :meth:`Viewer.panel() <pygeostat.plotting3d.Viewer.panel>`. Everything
it does is also a method of the viewer, so a view set up by clicking can be reproduced by a script.

.. autoclass:: pygeostat.plotting3d.panel.Panel
   :members: refresh

Meshes
******

These convert pygeostat data to PyVista meshes. :func:`view3d() <pygeostat.plotting3d.view3d>` calls them,
and they can be used alone to work with the mesh in PyVista.

.. autofunction:: pygeostat.plotting3d.drillhole_mesh

.. autofunction:: pygeostat.plotting3d.block_mesh
