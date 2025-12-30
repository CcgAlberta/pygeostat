.. _datautils:

Data Utilities
##############

The :mod:`pygeostat.datautils` module includes data utility functions for desurveying drillhole data,
compositing, labeling, and other data manipulation tasks.

.. currentmodule:: pygeostat.datautils

Desurveying
===========

Drillhole Class
---------------
.. autoclass:: Drillhole
   :members:

set_desurvey
------------
.. autofunction:: set_desurvey

get_desurvey
------------
.. autofunction:: get_desurvey

Compositing
===========

set_comps
---------
.. autofunction:: set_comps

get_comps
---------
.. autofunction:: get_comps

fast_comps
----------
.. autofunction:: fast_comps

Labeling
========

insert_real_idx
---------------
.. autofunction:: insert_real_idx

make_labels
-----------
.. autofunction:: make_labels

Utility Functions
=================

round_sigfig
------------
.. autofunction:: round_sigfig

fileheader
----------
.. autofunction:: fileheader

corrmatstr
----------
.. autofunction:: corrmatstr

slice_grid
----------
.. autofunction:: slice_grid

slicescatter
------------
.. autofunction:: slicescatter

fixpath
-------
.. autofunction:: fixpath

is_numeric
----------
.. autofunction:: is_numeric

ensure_dir
----------
.. autofunction:: ensure_dir

ensure_path
-----------
.. autofunction:: ensure_path

nearest_eucdist
---------------
.. autofunction:: nearest_eucdist

VTK Export
==========

write_vtp
---------
.. autofunction:: write_vtp
