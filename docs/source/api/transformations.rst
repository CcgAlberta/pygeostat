.. _transformations:

Transformations
###############

The :mod:`pygeostat.transformations` module provides functions for calculating GSLIB rotations,
coordinate transformations, and 3D geometric operations.

.. currentmodule:: pygeostat.transformations

GSLIB Rotations
===============

get_rotation_matrix
-------------------
.. autofunction:: get_rotation_matrix

principalvectors
----------------
.. autofunction:: principalvectors

azmdip
------
.. autofunction:: azmdip

principaldirs
-------------
.. autofunction:: principaldirs

Visualization
=============

Arrow3D Class
-------------
.. autoclass:: Arrow3D
   :members:

plotprincipalvectors
--------------------
.. autofunction:: plotprincipalvectors

drawgsvectorwidget
------------------
.. autofunction:: drawgsvectorwidget

drawellipsoid
-------------
.. autofunction:: drawellipsoid

drawgsaniswidget
----------------
.. autofunction:: drawgsaniswidget

Matrix Operations
=================

identity_matrix
---------------
.. autofunction:: identity_matrix

inverse_matrix
--------------
.. autofunction:: inverse_matrix

concatenate_matrices
--------------------
.. autofunction:: concatenate_matrices

decompose_matrix
----------------
.. autofunction:: decompose_matrix

compose_matrix
--------------
.. autofunction:: compose_matrix

is_same_transform
-----------------
.. autofunction:: is_same_transform

Translation
===========

translation_matrix
------------------
.. autofunction:: translation_matrix

translation_from_matrix
-----------------------
.. autofunction:: translation_from_matrix

Rotation
========

rotation_matrix
---------------
.. autofunction:: rotation_matrix

rotation_from_matrix
--------------------
.. autofunction:: rotation_from_matrix

random_rotation_matrix
----------------------
.. autofunction:: random_rotation_matrix

Reflection
==========

reflection_matrix
-----------------
.. autofunction:: reflection_matrix

reflection_from_matrix
----------------------
.. autofunction:: reflection_from_matrix

Scaling
=======

scale_matrix
------------
.. autofunction:: scale_matrix

scale_from_matrix
-----------------
.. autofunction:: scale_from_matrix

Projection
==========

projection_matrix
-----------------
.. autofunction:: projection_matrix

projection_from_matrix
----------------------
.. autofunction:: projection_from_matrix

clip_matrix
-----------
.. autofunction:: clip_matrix

Shear
=====

shear_matrix
------------
.. autofunction:: shear_matrix

shear_from_matrix
-----------------
.. autofunction:: shear_from_matrix

Affine Transformations
======================

affine_matrix_from_points
-------------------------
.. autofunction:: affine_matrix_from_points

superimposition_matrix
----------------------
.. autofunction:: superimposition_matrix

orthogonalization_matrix
------------------------
.. autofunction:: orthogonalization_matrix

Euler Angles
============

euler_matrix
------------
.. autofunction:: euler_matrix

euler_from_matrix
-----------------
.. autofunction:: euler_from_matrix

euler_from_quaternion
---------------------
.. autofunction:: euler_from_quaternion

Quaternions
===========

quaternion_from_euler
---------------------
.. autofunction:: quaternion_from_euler

quaternion_about_axis
---------------------
.. autofunction:: quaternion_about_axis

quaternion_matrix
-----------------
.. autofunction:: quaternion_matrix

quaternion_from_matrix
----------------------
.. autofunction:: quaternion_from_matrix

quaternion_multiply
-------------------
.. autofunction:: quaternion_multiply

quaternion_conjugate
--------------------
.. autofunction:: quaternion_conjugate

quaternion_inverse
------------------
.. autofunction:: quaternion_inverse

quaternion_real
---------------
.. autofunction:: quaternion_real

quaternion_imag
---------------
.. autofunction:: quaternion_imag

quaternion_slerp
----------------
.. autofunction:: quaternion_slerp

random_quaternion
-----------------
.. autofunction:: random_quaternion

Arcball
=======

Arcball Class
-------------
.. autoclass:: Arcball
   :members:

arcball_map_to_sphere
---------------------
.. autofunction:: arcball_map_to_sphere

arcball_constrain_to_axis
-------------------------
.. autofunction:: arcball_constrain_to_axis

arcball_nearest_axis
--------------------
.. autofunction:: arcball_nearest_axis

Vector Operations
=================

vector_norm
-----------
.. autofunction:: vector_norm

unit_vector
-----------
.. autofunction:: unit_vector

random_vector
-------------
.. autofunction:: random_vector

vector_product
--------------
.. autofunction:: vector_product

angle_between_vectors
---------------------
.. autofunction:: angle_between_vectors
