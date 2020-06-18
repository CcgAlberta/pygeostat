.. _data:

Data Files
##########

The core class in pygeostat is the ``DataFile`` class which contains a Pandas DataFrame with the
data values and column names in addition to metadata, such as the name of the x, y and z coordinates
or grid definition.


.. _datafile:

DataFile Class
===============
.. autoclass:: pygeostat.data.data.DataFile


DataFile Attributes
+++++++++++++++++++
Attributes of a ``datafile`` object are accessed with ``datafile.<attribute>``.

Columns
-------
Access the ``columns`` of the datafile. Wrapper for ``datafile.data.columns``.

Num Variables
-------------
Access the ``nvar`` of the datafile. e.g., the ``len(datafile.variables)``

Locations
---------
Access the ``locations`` stored in the datafile. Wrapper for ``datafile[datafile.xyz]``.

Example:
    
    >>> datafile = gs.DataFile("somefile.out")  # this file has an x, y[, z] attribute that is found
    >>> datafile.locations
    ... dataframe of x, y, z locations

Shape
-----
Access the ``shape`` of the data stored in the datafile. Wrapper for ``datafile.data.shape``

Example:

    >>> datafile = gs.DataFile("somefile.out")
    >>> datafile.shape
    ... shape of datafile.data



Head
++++
.. automethod:: pygeostat.data.data.DataFile.head

Rename Columns
++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.rename

Drop Columns
++++++++++++
.. automethod:: pygeostat.data.data.DataFile.drop

Check for Duplicate Columns
+++++++++++++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.check_for_duplicate_cols

Set Columns
+++++++++++
.. automethod:: pygeostat.data.data.DataFile.setcol

Set Variable Columns
++++++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.setvarcols

Set Categorical Dictionary
++++++++++++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.setcatdict

Check DataFile
++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.check_datafile

Add Coord
+++++++++
.. automethod:: pygeostat.data.data.DataFile.addcoord

Apply Dictionary
++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.applydict

Describe DataFile
+++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.describe

Infer Grid Definition
++++++++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.infergriddef

File Name String
++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.__str__

Generate Dictionary
+++++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.gendict

GSLIB Column
++++++++++++
.. automethod:: pygeostat.data.data.DataFile.gscol

Truncate NaN's
++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.truncatenans

Unique Categories
+++++++++++++++++
.. automethod:: pygeostat.data.data.DataFile.unique_cats

Write file
++++++++++++
.. automethod:: pygeostat.data.data.DataFile.write_file

Data Spacing
++++++++++++
.. automethod:: pygeostat.data.data.DataFile.spacing


Example Data
============
.. autofunction:: pygeostat.data.data.ExampleData



.. _iotools:

Input/Ouput Tools
=================
.. automodule:: pygeostat.data.iotools

Read File
+++++++++
.. autofunction:: pygeostat.data.iotools.read_file

Read CSV
++++++++
.. autofunction:: pygeostat.data.iotools.read_csv


Read GSLIB Python
++++++++++++++++++++++++++
.. autofunction:: pygeostat.data.iotools.read_gslib



Write GSLIB Python
+++++++++++++++++++++++++++
.. autofunction:: pygeostat.data.iotools.write_gslib

Write CSV
+++++++++
.. autofunction:: pygeostat.data.iotools.write_csv


Write VTK
+++++++++
.. autofunction:: pygeostat.data.iotools.write_vtk

Write HDF5 VTK
++++++++++++++
.. autofunction:: pygeostat.data.iotools.write_hvtk

Count Lines in File
+++++++++++++++++++
.. autofunction:: pygeostat.data.iotools.file_nlines

Write CCG GMM
+++++++++++++
.. autofunction:: pygeostat.data.iotools.writeout_gslib_gmm




.. _hdf5:

HDF5 I/O
========

Write HDF5
++++++++++
.. autofunction:: pygeostat.data.h5_io.write_h5

Read HDF5
+++++++++
.. autofunction:: pygeostat.data.h5_io.read_h5

Is HDF5
++++++++
.. autofunction:: pygeostat.data.h5_io.ish5dataset

Combine Datasets from Multiple Paths
++++++++++++++++++++++++++++++++++++
.. autofunction:: pygeostat.data.h5_io.h5_combine_data


Pygeostat HDF5 Class
++++++++++++++++++++
.. autoclass:: pygeostat.data.h5_io.H5Store

Write Data
----------
.. automethod:: pygeostat.data.h5_io.H5Store.__setitem__

Read Data
---------
.. automethod:: pygeostat.data.h5_io.H5Store.__getitem__

Print Contents of HDF5 File
---------------------------
.. automethod:: pygeostat.data.h5_io.H5Store.__str__

Close the HDF5 File
-------------------
.. automethod:: pygeostat.data.h5_io.H5Store.close

Datasets in H5 Store
--------------------
.. automethod:: pygeostat.data.H5Store.datasets

Generate Iterator
-----------------
.. automethod:: pygeostat.data.H5Store.iteritems



DictFile Class
==============
.. autoclass:: pygeostat.data.data.DictFile

Read Dictionary
+++++++++++++++
.. automethod:: pygeostat.data.data.DictFile.read_dict

Write Dictionary
++++++++++++++++
.. automethod:: pygeostat.data.data.DictFile.write_dict
