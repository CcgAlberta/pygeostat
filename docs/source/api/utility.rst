.. _utility:

Utility Functions
#################

The :mod:`pygeostat.utility` module provides utility functions for file management,
logging, and deprecation warnings.

.. currentmodule:: pygeostat.utility

File Management
===============

mkdir
-----
.. autofunction:: mkdir

rmdir
-----
.. autofunction:: rmdir

rmfile
------
.. autofunction:: rmfile

get_executable
--------------
.. autofunction:: get_executable

list_executable
---------------
.. autofunction:: list_executable

Logging
=======

printerr
--------
.. autofunction:: printerr

log_progress
------------
.. autofunction:: log_progress

Deprecation
===========

PygeostatDeprecationWarning
---------------------------
.. autoclass:: PygeostatDeprecationWarning
   :members:

warn_deprecated
---------------
.. autofunction:: warn_deprecated

deprecated
----------
.. autofunction:: deprecated
