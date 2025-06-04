#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
This module Contains basic data class to read data in different formats including GeoEAS, csv and hdf5
'''

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import csv
import os
import numpy as np
import pandas as pd

from . import iotools as iotools
from ..utility.logging import printerr
from .. pygeostat_parameters import Parameters

# Private dictionary of checks that are used in setcol for the specialized column attributes
_CHECKS = {'dh': ['dh', 'drillhole', 'drill hole', 'dhid', 'drillholeid', 'drill hole id', 'Drillhole Number',
				  'well', 'wellid', 'well id', 'holeid', 'hole id', 'id', 'ID'],
		   'ifrom':  ['ifrom', 'from', 'from depth', 'ifrom(m)', 'from(m)', 'from depth(m)'],
		   'ito': ['ito', 'to', 'to depth', 'ito(m)', 'to(m)',
				   'To(m)', 'TO(m)', 'to depth(m)', 'To Depth(m)', 'TO DEPTH(m)'],
		   'x': ['x', 'east', 'easting', 'x(m)', 'east(m)', 'easting(m)', 'xlocation'],
		   'y': ['y', 'north', 'northing', 'y(m)', 'north(m)', 'northing(m)', 'ylocation'],
		   'z': ['z', 'elev', 'elevation', 'z(m)', 'elev(m)', 'elevation(m)', 'zlocation'],
		   'cat': ['cat', 'category', 'categories', 'rt', 'rocktype', 'rock type', 'fac',
				   'facies', 'lithofacies', 'lithology', 'Facies Code'],
		   'weights': ['wt', 'wts', 'weights', 'weight', 'declusteringweight', 'declustering weight']}


class DataFile:
	"""
	This class stores geostatistical data values and metadata.

	DataFile classes may be created on initialization, or generated
	using pygeostat functions. This is the primary class for pygeostat
	and is used for reading and writing GSLIB, CSV, VTK, and HDF5 file formats.

	Parameters:
		flname (str): Path (or name) of file to read
		readfl (bool): True if the data file should be read on class initialization
		fltype (str): Type of data file: either ``csv``, ``gslib`` or ``hdf5`` or ``gsb``
		dftype (str): Data file type as either 'point' or 'grid' used for writing out VTK files
			for visualization
		data (pandas.DataFrame): Pandas dataframe containing array of data values
		dicts (List[dict] or dict): List of dictionaries or dictionary for converting alphanumeric
			to numeric data
		null (float): Null value for missing values in the data file
		title (str): Title, or name, of the data file
		griddef (pygeostat.GridDef): Grid definition for a gridded data file
		dh (str): Name of drill hole variable
		x (str): Name of X coordinate column
		y (str): Name of Y coordinate column
		z (str): Name of Z coordinate column
		ifrom (str): Name of 'from' columns
		ito (str): Name of 'to' columns
		weights (str or list): Name of declustering weight column(s)
		cat (str): Name of categorical (e.g., rock type or facies) column
		catdict (dict): Set a dictionary for the categories, which should be formatted as:
			``catdict = {catcode:catname}``
		variables (str or list): Name of continuous variable(s), which if unspecified, are the
			columns not assigned to the above attributes (via kwargs or inference)
		notvariables (str or list): Name of column(s) to exclude from variables
		delimiter (str): Delimiter used in data file (ie: comma or space)
		headeronly (bool): True to just read header + 1 line of data file This is useful for
			getting column numbers of large files OR if reading hdf5 files will only read in the
			hdf5 store information
		h5path (str): Forward slash (/) delimited path through the group hierarchy you wish to
			read the dataset(s) specified by the argument ``datasets`` from. The dataset name
			cannot be passed using this argument, it is interpreted as a group name only. A
			value of ``None`` places the dataset into the root directory of the HDF5 file. A value
			of ``False`` loads a blank pd.DataFrame().
		h5datasets (str or list): Name of the dataset(s) to read from the group specified by
			``h5path``. Does nothing if ``h5path`` points to a dataset.
		columns (list): List of column labels to use for the resulting ``data`` pd.DataFrame
		nreals (int): number of realizations to read in. -1 will read all
		tmin (float): If a number is provided, values less than this number (e.g., trimmed
					  or null values) are convernted to NaN. May be useful since NaN's are
					  more easily handled within python, matplotlib and pandas. Set to None
					  to disable.

	Examples:

		Quickly reading in a GeoEAS data file:

		.. code-block:: python

			data_file = gs.DataFile(flname='../data/oilsands.dat')


		To read in a GeoEAS datafile and assign attributes

		.. code-block:: python

			# Point Data Example
			data_file = gs.DataFile(flname='../data/oilsands.dat',readfl=True,dh='Drillhole Number', x='East',y='North',z='Elevation')

		.. code-block:: python

			# Gridded Data Example
			griddef = gs.GridDef('''10 0.5 1
			10 0.5 1
			10 0.5 1''')
			data_file = gs.DataFile(flname='../data/3DDecor.dat', griddef=griddef)

			# To view grid definition string
			print(data_file.griddef)
			# Access some Grid Deffinition attributes
			data_file.griddef.count() # returns number of blocks in grid
			data_file.griddef.extents() # returns an array of the extents for all directions
			data_file.griddef.nx # returns nunmber of blocks in x direction

	**HDF5**

	Using the HDF5 file format has its own positive features. For one it reads and writes much
	faster then using the ASCII format. Attributes (like the grid definition) can also be saved
	within the file. All files for a single project can also be saved in the same file. Please
	refer to the :ref:`introduction on HDF5 files <hdf5>` for more information

	This class currently only searches for and loads a grid definition.

	Examples:
		
		HDF5 file simple read example:

		.. code-block:: python

			data_file = gs.DataFile(flname='../data/oilsands_out.hdf5')

		To view the HDF5 header information (tables stored in the file):

		.. code-block:: python

			data_file.store

		If you have a HDF5 file with multiple tables and you just want to read in the file
		information to view what tables are in the file and any attributes saved to the file you
		can do a header style only read:

		.. code-block:: python

			data_file = gs.DataFile(flname='../data/oilsands_out.hdf5', dftype='hdf5', headeronly=True)

		Then to see what tables are written in the hdf5 file:

		.. code-block:: python

			data_file.store
	"""

	def __init__(self, flname=None, readfl=None, fltype=None, dftype=None, data=None, columns=None,
				 null=None, title='data', griddef=None, dh=None, x=None, y=None, z=None, ifrom=None,
				 ito=None, weights=None, cat=None, catdict=None, variables=None, notvariables=None,
				 delimiter=r"\s+", headeronly=False, h5path=None, h5datasets=None, nreals=-1,
				 tmin=None):
		self.flname = flname
		self.readfl = readfl
		# File types include 'GSLIB' and 'CSV' and 'HDF5'
		self.fltype = fltype
		self.headeronly = headeronly
		# Assign a default delimiter
		self.delimiter = delimiter
		self.data = data
		if null is None:
			self.null = Parameters['data.null']
		else:
			self.null = null
		self.title = title
		# Note that the default grid definition will be checked against Parameters after
		# initializing DataFile.data
		self.griddef = griddef
		self.store = None
		self.h5fl = False
		self.nreals = nreals
		# Read the data in if ``data`` is None but a file is passed
		if (data is None) and (flname is not None) and (readfl is None):
			self.readfl = True
		if (data is None) and (readfl is None):
			self.data = pd.DataFrame()
		# Data types include 'point', 'grid' - not that grid is assigned later after
		# default griddefs are considered... which may override the point assignment below
		if dftype is not None:
			self.dftype = dftype
		else:
			self.dftype = 'point'
		# Read the file
		if self.readfl:
			_, extention = os.path.splitext(flname)
			extention = extention.lower()
			if extention in ['.h5', '.hdf5', '.hd5']:
				self.h5fl = True
			from .iotools import read_file
			# Load the h5fl as a H5Store class as required
			if self.h5fl:
				from .h5_io import H5Store
				self.store = H5Store(flname)
			# read the file and import data and/or header
			self.data = read_file(self.flname, fltype=fltype, headeronly=headeronly,
								 delimiter=delimiter, h5path=h5path, h5datasets=h5datasets,
								 columns=columns, griddef=self.griddef, ireal=nreals, tmin=tmin)
		# Try and coerce data into a pandas dataframe if it is not already
		if (data is not None) and (not isinstance(data, pd.DataFrame)):
			try:
				self.data = pd.DataFrame(data=data, columns=columns)
			except:
				raise ValueError("Please ensure the passed `data` can be coerced into a pandas"
								 " dataframe (i.e., data = pd.DataFrame(data=data)")
		# Initialize the specialized columns
		self.dh = dh
		self.ifrom = ifrom
		self.ito = ito
		self.x = x
		self.y = y
		self.z = z
		self.weights = weights
		self.cat = cat
		self.variables = variables
		# Will set/validate specialized columns by iterating through this dictionary
		_specialattr = self._get_specialattr()
		if self.data is not None:
			# Check for duplicate column names
			self.check_for_duplicate_cols()
			# Set specialized columns, validating their existence DataFile.data if user provided,
			# or assigning based on common names if present in DataFile.data
			for colattr, colval in _specialattr.items():
				self.setcol(colattr, colval)
				# Get variables
				self.setvarcols(variables, notvariables)
		else:  # data is None, so `self.variables` has no meaning
			if self.variables is not None:
				# Simply assign the names, where the only validation is that a string is provided
				for colattr, colval in _specialattr.items():
					if not isinstance(colval, str) and colval is not None:
						raise TypeError(colattr + ' should be a list or tuple or string!')
					setattr(self, colattr, colval)
				if isinstance(self.variables, str):
					self.variables = [variables]
				elif isinstance(self.variables, tuple):
					self.variables = list(variables)
				elif not isinstance(self.variables, list):
					raise TypeError('variables should be a list or tuple or string!')

		# If a hdf5 file was read in, check to see if a griddef can be found
		if self.readfl and self.h5fl:
			# Handle some input parameters
			if h5path in [None, '', ' ', '/']:
				h5path = '/'
			if isinstance(h5datasets, str):
				h5datasets = [h5datasets]
			if h5datasets is None:
				h5datasets = ['']
			gridstrs = []
			# If multiple datasets are loaded, check them all for a griddef
			for dataset in h5datasets:
				attrs = self.store.h5data[h5path + '/' + dataset].attrs
				if 'griddef' in attrs:
					gridstrs.append(attrs['griddef'])
			# Set the griddef if only one or multiple griddefs were found that are all the same;
			# otherwise, don't load anything
			if len(gridstrs) > 0:
				from .grid_definition import GridDef
				from ..utility.logging import printerr
				if len(np.unique(gridstrs)) > 1:
					printerr("Multiple grid definitions are in the HDF5 file. Nothing has been"
							 " set.", errtype='error')
				else:
					gridstr = gridstrs[0].replace("b'", "'").replace("'", "").replace("\\n", "\n")
					if self.griddef is not None:
						printerr("Grid definition is being rewritten by griddef saved in h5 file",
								 errtype='warning')

					self.griddef = GridDef(gridstr)
		# If griddef is None, check if a default griddef is initialized - assign it to this
		# DataFile if its length matches. Convert the fltype if the user hasn't specified
		# a fltype
		if self.data is not None and griddef is None and Parameters['data.griddef'] is not None:
			if Parameters['data.griddef'].count() == self.data.shape[0]:
				self.griddef = Parameters['data.griddef']
		if dftype is None and self.griddef is not None:
			if self.griddef.count() == self.data.shape[0]:
				self.dftype = 'grid'
		# Set the catdict if appropriate
		self.catdict = None
		if self.cat is not None:
			if catdict is None:
				catdict = Parameters['data.catdict']
			if catdict is not None:
				try:
					self.setcatdict(catdict)
				except:
					# The catdict will not be set if catdict.keys() doesn't match self.cat
					pass

	def __str__(self):
		"""
		Return the name of the data file if asked to 'print' the data file... or use the datafile
		in a string!
		"""

		return str(self.flname)

	def __len__(self):
		"""
		Return the length of the dataframe this DataFile object contains
		"""
		if self.data is not None:
			return len(self.data)
		else:
			return 0

	def __getitem__(self, col):
		"""
		Access the column of data corresponding to `col` in the self.data `DataFrame`

		"""
		from copy import deepcopy
		assert(isinstance(self.data, pd.DataFrame)), "`obj.data` is not a `pd.DataFrame`"
		if isinstance(col, (bool, pd.core.series.Series)):
			new_datafile = deepcopy(self)
			try:
				new_datafile.data = self.data[col]
			except Exception as exception_message:
				raise TypeError(str(exception_message))
			return new_datafile

		assert(isinstance(col, (str, list))), \
			"`DataFile` can be indexed by `str` or `list` of columns"
		if isinstance(col, (list, tuple)):
			for c in col:
				assert(c in self.data.columns), "`{}` is not in this file!".format(col)
		else:
			assert(col in self.data.columns), "`{}` is not in this file!".format(col)
		return self.data[col]

	def __setitem__(self, col, value):
		if self.data is None:
			self.data = pd.DataFrame()
		self.data.loc[:, col] = value


	def _get_specialattr(self):
		'''Return a dictionary of the current special attribute values'''
		specialattr = {'dh': self.dh, 'ifrom': self.ifrom, 'ito': self.ito, 'x': self.x,
					   'y': self.y, 'z': self.z, 'cat': self.cat, 'weights': self.weights}
		return specialattr

	def _get_shape(self):
		assert isinstance(self.data, pd.DataFrame), "DataFile.data must be a pandas.DataFrame!"
		return self.data.shape
	shape = property(_get_shape)

	def head(self, n=5):
		'''Return the first n rows of the data, accessing self.data.head()'''
		assert isinstance(self.data, pd.DataFrame), "DataFile.data must be a pandas.DataFrame!"
		return self.data.head(n)

	def tail(self, n=5):
		'''Return the last n rows of the data, accessing self.data.tail()'''
		assert isinstance(self.data, pd.DataFrame), "DataFile.data must be a pandas.DataFrame!"
		return self.data.tail(n)

	def _get_columns(self):
		""" Get the columns of the pandas data frame"""
		return self.data.columns

	columns = property(_get_columns)

	def _get_locations(self):
		return self[self.xyz]
	locations = property(_get_locations)

	def _get_info(self):
		"""
		Print a summary of the special attributes found in this datafile

		Parameters:
			verbose (bool): if True prints the datafile.data.info() showing the types and counts
				of each column in the dataframe

		"""
		printstr = 'DataFile: {}\n'.format(self.flname)
		specialattr = self._get_specialattr()
		if not all([v is None for v in specialattr.values()]):
			printstr += 'Attributes:\n'
			for attr, value in specialattr.items():
				if value is not None:
					printstr += "{}: '{}',  ".format(attr, value)
			printstr += '\n'
		else:
			printstr += 'No Special Attributes Found \n'
		if self.data is not None:
			if self.variables is not None:
				printstr += 'Variables:\n'
				if isinstance(self.variables, str):
					printstr += self.variables
				else:
					printstr += ', '.join(["'{}'".format(var) for var in self.variables])
		if self.griddef is not None:
			printstr += '\nGrid Definitions:\n'
			printstr += str(self.griddef)
		return printstr

	info = property(_get_info)

	def rename(self, columns):
		'''
		Applies a dictionary to alter self.DataFrame column names. This applies the DataFrame.rename
		function, but updates any special attributes (dh, x, y, etc.) with the new name, if previously
		set to the old name. Users should consider using the self.columns property if changing all
		column names.

		Parameters:
			columns(dict): formatted as {oldname1: newname1, oldname2:newname2}, etc, where the old and
				new names are strings. The old names must be present in data.columns.

		'''
		# if this column already exists... drop it before renaming to avoid dups?
		for old, new in columns.items():
			if new in self.data.columns:
				self.data.drop(new, axis=1, inplace=True)
		self.data.rename(index=str, columns=columns, inplace=True)
		# Update the special attributes
		specialattr = self._get_specialattr()
		for old, new in columns.items():
			found = False
			for colattr, colval in specialattr.items():
				if isinstance(colval, str):
					colval = [colval]
				if isinstance(colval, list):
					for i in range(len(colval)):
						if colval[i] == old:
							colval[i] = new
							found = True
							break
					if found:
						if len(colval) < 1:
							colval = None
						elif len(colval) == 1:
							colval = colval[0]
						self.setcol(colattr, colval)
						break
		# Update the variables
		if isinstance(self.variables, list):
			variables = self.variables
		elif isinstance(self.variables, str):
			variables = [self.variables]
		else:
			variables = None
		if isinstance(variables, list):
			for old, new in columns.items():
				if old in variables:
					i = variables.index(old)
					variables[i] = new
			if len(variables) == 1:
				variables = variables[0]
			self.variables = variables

	def drop(self, columns):
		'''
		This applies the DataFrame.drop function, where axis=1, inplace=True and columns is
		used in place of the labels. It also updates any special attributes (dh, x, y, etc.),
		setting them to None if dropped. Similarly, if any variables are dropped, they are
		removed from self.variables.

		Parameters:
			columns(str or list): column names to drop

		'''
		if isinstance(columns, str):
			columns = [columns]
		self.data.drop(axis=1, labels=columns, inplace=True)
		# Set the special attributes to None if in columns
		specialattr = self._get_specialattr()
		for key, val in specialattr.items():
			if isinstance(val, str):
				val = [val]
			if isinstance(val, list):
				temp = []
				for v in val:
					if v not in columns:
						temp.append(v)
				if len(temp) < 1:
					temp = None
				elif len(temp) == 1:
					temp = temp[0]
				setattr(self, key, temp)
		# Update the variables
		if isinstance(self.variables, list):
			variables = self.variables
		elif isinstance(self.variables, str):
			variables = [self.variables]
		else:
			variables = None
		if isinstance(variables, list):
			temp = []
			for var in variables:
				if var not in columns:
					temp.append(var)
			if len(temp) < 1:
				temp = None
			elif len(temp) == 1:
				temp = temp[0]
			self.variables = temp

	def _get_xyz(self):
		assert(self.x is not None), 'No x column has been assigned for the current instance'
		assert(self.y is not None), 'No y column has been assigned for the current instance'
		if self.z is None:
			return [self.x, self.y, None]
		else:
			return [self.x, self.y, self.z]
	xyz = property(_get_xyz)

	def check_for_duplicate_cols(self):
		"""
		Run a quick check on the column names to see if any of them are duplicated. If they
		are duplicated then print a Warning and rename the columns

		"""

		cols=pd.Series(self.data.columns)

		# Check for duplicates by comparing unique list to full column list
		if any(cols.duplicated()):
			from ..utility.logging import printerr
			printerr("Duplicate column names found when reading in data!!",
					 errtype='warning')

			for dup in cols[cols.duplicated()].unique():
				cols[cols[cols == dup].index.values.tolist()] = [dup + '_' + str(i) if i != 0 else dup for i in range(sum(cols == dup))]
			self.data.columns=cols
			# print out what we changed the column names to.
			print('New unique column names are: ', list(cols))

	def setcol(self, colattr, colname=None):
		"""
		Set a specialized column attribute (``dh, ifrom, ito, x, y, z, cat or weights``) for the
		``DataFile``, where ``DataFile.data`` must be initialized. If colname is None, then the
		attribute will be set if a common name for it is detected in ``DataFile.data``.columns
		(e.g., if ``colattr='dh'`` and ``colname=None``, and ``'DHID'`` is found in
		``DataFile.data``, then ``DataFile.dh='DHID'``. The attribute will be None if none of the
		common names are detected. If colname is not None, then the provided string will be
		assigned to the attribute, e.g. DataFile.colattr=colname. Note, however, that an error
		will be thrown if colname is not None and colname is not in ``DataFile.data``.columns.
		This is used on DataFile initialization, but may also be useful for calling after
		specialized columns are altered.

		Parameters:
			colattr(str) : must match one of: ``'dh'``, ``'ifrom'``, ``'ito'``, ``'x'``, ``'y'``,
				``'z'``, ``'cat'`` or ``'weights'``
			colname(str or list) : if not None, must be the name(s) of a column in ``DataFile.data``.
				List is only valid if ``colattr=weights``

		Examples:

			Set the x attribute (dat.x) based on a specified value:

			>>> dat.setcol('x', 'Easting')

			Set the x attribute (dat.x), where the function checks common names for x:

			>>> dat.setcol('x')

		"""
		if self.data is None:
			raise ValueError('function can only be called after DataFile.data is initialized!')
		isset = False
		columns = []
		for s in list(self.data.columns):
			try:
				columns.append(s.lower())
			except AttributeError:
				columns.append(s)
		# Note that weights is the single attribute that is a list - so it's handled in a non-general
		# way for now. Should probably be cleaned up.
		if colname is None:
			# No attribute name is provided, so the function will check common
			# naming conventions for this attribute
			# First, grab the list of names to check
			checks = _CHECKS[colattr]
			if Parameters['data.'+colattr] is not None:
				# If a default naming convention for this attribute was set for this project,
				# then add it to the start of the checks list
				checks = [Parameters['data.'+colattr].lower()] + checks
			weights = []
			for check in checks:
				i = [i for i, x in enumerate(columns) if x.lower() == check.lower()]
				if len(i) > 0:
					i = i[0]
					colname = self.data.columns[i]
					# Set the attribute since the checked name is in the columns
					if colattr != 'weights':
						setattr(self, colattr, colname)
						isset = True
						break
					else:
						weights.append(colname)
			if len(weights) > 0:
				if len(weights) == 1:
					weights = weights[0]
				setattr(self, colattr, weights)
				isset = True
		else:
			if isinstance(colname, str):
				colnames = [colname]
			elif isinstance(colname, tuple) or isinstance(colname, list):
				if colattr != 'weights':
					raise ValueError('colname must be a str unless weights!')
				colnames = colname
			else:
				raise ValueError(colname+' is invalid!')
			for colname in colnames:
				if colname not in self.data.columns:
					raise ValueError(colname + ' is not a column in the data!')
			if len(colnames) == 1:
				setattr(self, colattr, colnames[0])
			else:
				setattr(self, colattr, colnames)
			isset = True
		# Throw a warning if the attr is a coordinate and has nans
		if isset:
			if any([colattr == coord for coord in ['x', 'y', 'z']]):
				if np.any(np.isnan(self.data[colname].values)):
					print('WARNING: null (NaN) values in the {} coordinate column!'.format(colname))
					print('     You may want to check the tmin setting (also, Parameters[data.tmin])')

	def setvarcols(self, variables=None, notvariables=None):
		"""
		Set the variables for the DataFile. If provided, the function checks that the
		variables are present in the DataFrame. If not provided, the function assigns columns
		that are not specified as the variables (``dh, x, y, z, rt, weights``), as well as a list of
		user specified notvariables.

		This is used on DataFile initialization, but may also be useful for calling after
		variables are added or removed.

		Parameters:
			variables(list or str) : list of strings
			notvariables(list or str) : list of strings

		Examples:

			Set the variables based on a specified list:

			>>> dat.setvarcols(variables=['Au', 'Carbon'])

			Set the variables based on the function excluding specialized columns (dh, x, y, etc.):

			>>> dat.setvarcols()

			Set the variables based on the function excluding specialized columns (dh, x, y, etc.),
			as well as a user specified list of what is not a variable:

			>>> dat.setvarcols(notvariables=['Data Spacing', Keyout'])

		"""
		self.variables = None
		if variables is None:
			# Construct a list of the specialized columns
			if notvariables is None:
				notvariables = []
			elif isinstance(notvariables, str):
				notvariables = [notvariables]
			elif isinstance(notvariables, tuple):
				notvariables = list(notvariables)
			elif not isinstance(notvariables, list):
				raise ValueError('variables should be a string or list of strings')
			# Grab a dictionary of the special attributes
			specialattr = self._get_specialattr()
			for val in specialattr.values():
				if val is not None:
					if isinstance(val, str):
						notvariables.append(val)
					else:
						for v in val:
							notvariables.append(v)
			# Variables are the columns not in notvariables
			self.variables = []
			for column in self.data.columns:
				if column not in notvariables:
					self.variables.append(column)
		else:
			if self.variables is None:
				self.variables = variables
			if not isinstance(self.variables, list) and not isinstance(self.variables, str):
				raise ValueError('variables should be a string or list of strings')
			if isinstance(self.variables, str):
				self.variables = [self.variables]
			for variable in self.variables:
				if variable not in self.data.columns:
					raise ValueError('{} is not data.data.columns!'.format(variable))
		# Convert to a string if the list is length-1
		if isinstance(self.variables, list):
			if len(self.variables) == 1:
				self.variables = self.variables[0]
			elif len(self.variables) == 0:
				self.variables = None
	
	# Read only property number of variables
	def _get_nvar(self):
		if isinstance(self.variables, list):
			nvar = len(self.variables)
		elif isinstance(self.variables, str):
			nvar = 1
		else:
			nvar = 0
		return nvar
	nvar = property(_get_nvar)

	def __repr__(self):
		return self.info

	def setcatdict(self, catdict):
		'''
		Set a dictionary for the categories, which should be formatted as:

			>>> catdict = {catcode:catname}

		Example:

			>>> catdict = {0: "Mudstone", 1: "Sandstone"}
			>>> self.setcatdict(catdict)

		'''
		if self.cat is None:
			raise ValueError('self.cat must be initialized!')
		# Ensure that catdict provides a key for each category in the data
		cats = self[self.cat]
		try:
			cats = np.unique(cats[np.isfinite(cats)])
		except TypeError as exc:
			raise TypeError('Make sure the provided categorical column, "{}", contains numerical values'.format(self.cat)) from exc
		
		for cat in cats:
			if cat not in catdict.keys():
				raise ValueError('{} in self.cat, but not in catdict.keys()!'.format(cat))
		self.catdict = catdict


	def check_datafile(self, flname, variables, sep, fltype):
		"""
		Run some quick checks on the DataFile before writing and grab info if not provided

		"""

		# Default to introspection
		if flname is None:
			flname = self.flname
		# Default separator is comma, unless the file is GSLIB-style
		if sep is None:
			if fltype is not None:
				if fltype.lower() == 'gslib':
					sep = ' '
				elif fltype.lower() == 'csv':
					sep = ','
			else:
				sep = ','
		# Check for single variable only and convert to list of length 1
		if isinstance(variables, str):
			variables = [variables]
		# Get variable names
		if variables is None:
			variables = self.data.columns.tolist()
		return flname, variables, sep, fltype

	def write_file(self, flname, title=None, variables=None, fmt=None, sep=None, fltype=None, data=None, h5path=None, griddef=None,
				null=None, tvar=None, nreals=1):
		"""Writes out a GSLIB-style, VTK, CSV, Excel (XLSX), HDF5 data file.

		Parameters:
			flname (str): Path (or name) of file to write out.

		Keyword Args:
			title (str): Title for output file.
			variables (List(str)): List of variables to write out if only a subset is desired.
			fmt (str): Format to use for floating point numbers.
			sep (str): Delimiter to use for file output, generally don't need to change.
			fltype (str): Type of file to write either ``gslib``, ``vtk``, ``csv``, ``xlsx``,
				or ``hdf5``.
			data (str): Subset of data to write out - cannot be used with variables
				option!
			h5path (str): The h5 group path to write data to (H5 filetype)
			griddef (obj): a gslib griddef object
            tvar (str): Name of variable to use for compression when NaNs exist within it
            nreals (int): number of realizations you are writing out (needed for GSB)
			null (float): If a number is provided, NaN numbers are converted to this value
				prior to writing. May be useful since NaN's are more easily handled
				within python and pandas than null values, but are not valid
				in GSLIB. Set to None to disable (but NaN's must be handled prior to
				this function call if so).

		Note:
			pygeostat.write_file is saved for backwards compatibility or as an overloaded class method.
			Current write functions can be called seperately with the functions listed below:

			>>> import pygeostat as gs
			>>> import pandas as pd
			>>> gs.write_gslib(gs.DataFile or pd.DataFrame)
			>>> gs.write_csv(gs.DataFile or pd.DataFrame)
			>>> gs.write_hdf5(gs.DataFile or pd.DataFrame)
			>>> gs.write_vtk(gs.DataFile or pd.DataFrame)
			>>> gs.write_gsb(gs.DataFile or pd.DataFrame)

		Note: 
		The GSB format is not specifically intended for general users of pygeostat. Some CCG programs use GSB that 
		is a compressed GSLIB-like binary data format that greatly reduces the computational expense.

			The following calls are equivalent:

			>>> data_file.write_file('testgslib.out')
			>>> data_file.write_file('testgsb.gsb')

			is equivalent to:

			>>> gs.write_gslib(data_file, 'testgslib.out')
			>>> gs.write_gsb(data_file, 'testgsb.gsb')

			and similar to:

			>>> gs.write_gslib(data_file.data, 'testgslib.out')
			>>> gs.write_gsb(data_file.data, 'testgsb.gsb')

		"""
		from .grid_definition import GridDef
		# Infer a filetype if none is specified
		# ------------------------------------
		if fltype is None:
			if flname.endswith('.dat') or flname.endswith('.out'):
				fltype = 'gslib'
			elif flname.endswith('.csv'):
				fltype = 'csv'
			elif flname.endswith('.vtk'):
				fltype = 'vtk'
			elif flname.endswith('.hvtk'):
				fltype = 'hvtk'
			elif flname.endswith('.hd5') or flname.endswith('.hdf5') or flname.endswith('.h5'):
				fltype = 'hdf5'
			elif flname.endswith('.gsb'):
				fltype = 'gsb'
			else:
				raise ValueError(" Output file type (fltype) cannot be inferred!")
		# Check the arguments and for any info not passed thats already assigned
		# ------------------------------------
		variables0 = variables  # Storing this for VTK due to conflicting definitions
		flname, variables, sep, fltype = self.check_datafile(flname, variables, sep, fltype)
		fltype = fltype.lower()

		# Check data
		if null is None and fltype != 'vtk':
			null = self.null
		if data is None:
			data = self.data
		# Check GridDef
		if griddef is None and self.griddef is not None:
			griddef = self.griddef
		# For some reason, the following isinstance is returning an erroneous False. Any ideas?
		if griddef is not None:
			if not isinstance(griddef, GridDef):
				printerr("griddef is not an instance of pygeostat.GridDef Class and will not be written to file",
						 errtype='error')
				raise ValueError('griddef must be pygeostat.GridDef instance')
		# Based on the filetype call the correct iotools function
		# ------------------------------------
		if fltype == 'gslib':
			iotools.write_gslib(self.data, flname, title=self.title, variables=variables,
								fmt=fmt, sep=sep, null=null)
		elif fltype == 'csv' or fltype.lower() == 'xlsx' or fltype.lower() == 'excel':
			iotools.write_csv(self.data, flname, variables=variables,
							  fmt=fmt, sep=sep, fltype=fltype, null=null)
		elif fltype == 'hdf5' or fltype == 'h5' or fltype == 'hd5':
			from .h5_io import write_h5
			if isinstance(griddef, GridDef):
				griddef = griddef.__str__().encode()
			else:
				griddef = None
			# As with read and tmin, null is not yet implemented for h5
			write_h5(self.data, flname, h5path=h5path, datasets=variables, gridstr=griddef)
		elif fltype == 'vtk':
			iotools.write_vtk(self, flname, variables=variables0, null=null)
		elif fltype == 'hvtk':
			if griddef is None:
				printerr('ERROR: this DataFile must be a 3D gridded DataFile with a .griddef '
						 'attribute!', errtype='error')
				return
			else:
				griddef = self.griddef
			if flname.endswith('hvtk'):
				flname = flname.replace('.hvtk', '.xdmf')
			iotools.write_hvtk(self.data, flname, self.griddef, variables=variables)
		elif fltype == 'gsb':
			if not isinstance(fmt, list):
				fmt = 0
				iotools.write_gsb(self.data, flname, tvar=tvar, variables=variables,
                              nreals=nreals, fmt=fmt, griddef=griddef)
		else:
			printerr('Unsupported File type. Try with either "gslib"'
					 '"csv", "xlsx", "excel", "hdf5", or "vtk"', errtype='error')

	def gscol(self, variables, string=True):
		"""Returns the GSLIB (1-ordered) column given a (list of) variable(s).

		Parameters:
			variables (str or List(str)): Path, or name, of the data file.

		Keyword Args:
			string (bool): If True returns the columns as a string.

		Returns:
			cols (int or List(int) or string): GSLIB 1-ordered column(s).

		Note:
			None input returns a 0, which may be necessary, for example, with 2-D
			data:
			>>> data.xyz
			... ['East', 'North', None]
			>>> data.gscol(data.xyz)
			... '2 3 0'

		Examples:

			Some simple calls

			>>> data_file.gscol('Bitumen')
			... 5

			>>> data_file.gscol(['Bitumen', 'Fines'])
			... [5, 6]

			>>> data_file.gscol(['Bitumen', 'Fines'], string=True)
			... '5 6'
		"""

		if isinstance(variables, str) or variables is None:
			variables = [variables]
		cols = []
		for var in variables:
			if var is None:
				cols.append(0)
			else:
				cols.append(self.data.columns.get_loc(var) + 1)
		if len(cols) == 1:
			cols = cols[0]
		if not string:
			return cols
		else:
			if not isinstance(cols, list):
				return str(cols)
			else:
				return ' '.join(list(map(str, cols)))


	def describe(self, variables=None):
		"""
		Describe a data set using pandas describe(), but exclude special variables.

		Keyword Args:
			variables (List(str)): List of variables to describe.

		Returns:
			self.data[variables].describe(): Pandas description of variables.

		Examples:

			Describe all none special variables in the DataFrame (will exclued columns set as dh ID,
			coordinate columns, etc.)

			>>> data_file.describe()

			Or describe specific variables

			>>> data_file.describe(['Bitumen', 'Fines'])

		"""

		# Get a list of variables to set to np.nan
		if variables is None:
			if self.variables is None:
				# Should rarely happen... but just in case
				self.setvarcols()
			variables = self.variables
		if variables is None:
			variables = list(self.columns)
		return self.data[variables].describe()

	def copy(self):
		"""
		Copy the DataFile, which provides frequent utility in workflows.
		"""
		from copy import deepcopy
		return deepcopy(self)

	def gendict(self, var, outvar=None):
		"""
		Generates a dictionary with unique IDs from alphanumeric IDs.
		This is particularly useful for alphanumeric drill hole IDs
		which cannot be used in GSLIB software.

		Parameters:
			var (str): Variable to generate a dictionary for

		Keyword Args:
			outvar (str): Variable to generate using generated dictionary.

		Returns:
			newdict (dict): Dictionary of alphanumerics to numeric ids.

		Examples:

			A simple call

			>>> data_file.gendict('Drillhole')

			OR

			>>> dh_dict = data_file.gendict('Drillhole')
		"""

		# Get the unique IDs
		ids = self.unique_cats(var)
		# Build up reasonable numeric IDs given traditional DH labeling
		new_ids = []
		for idx, str_id in enumerate(ids):
			try:
				new_id = int(''.join(c for c in str_id if c.isdigit()))
			except:
				new_id = 0
			# Make it unique!
			while new_id in new_ids:
				new_id = int(''.join([str(new_id), str(idx)]))
			new_ids.append(new_id)
		# Generate the dictionary
		newdict = dict(zip(ids, new_ids))
		# Apply the dictionary if desired
		if outvar is not None:
			self.applydict(var, outvar, newdict)
		# Return the new dictionary
		return newdict

	def applydict(self, origvar, outvar, mydict):
		"""Applies a dictionary to the original variable to get a new variable.

		This is particularly useful for alphanumeric drill hole IDs
		which cannot be used in GSLIB software.

		Parameters:
			origvar (str): Name of original variable.
			outvar (str): Name of output variable.
			mydict (dict): Dictionary of values to apply.

		Examples:

			>>> data_file.applydict('Drillhole', 'Drillhole-mod', mydict)
		"""

		def dictapply(value):
			if value in mydict:
				return mydict[value]
			else:
				return value
		self.data[outvar] = self.data[origvar].apply(dictapply)

	def unique_cats(self, variable, truncatenans=False):
		"""Returns a sorted list of the unique categories given a variable.

		Parameters:
			variable (str): Name of original variable.

		Keyword Args:
			truncatenans (bool): Truncates missing values if True.

		Returns:
			unique_cats (List(object)): Sorted, list of set(object).

		Examples:

			A simple call that

			>>> data_file.unique_cats('Drillhole')

			Or to save the list

			>>> unique_dh_list = data_file.unique_cats('Drillhole')
		"""

		if truncatenans:
			variable = self.truncatenans(variable)
		if isinstance(variable, str):
			return sorted(list(set(self.data[variable])))
		else:
			try:
				return sorted(list(set(variable)))
			except:
				raise Exception('Could not get unique categories!')

	def truncatenans(self, variable):
		"""Returns a truncated list with nans removed for a variable.

		Parameters:
			variable (str): Name of original variable.

		Returns:
			truncated (values): Truncated values.

		Examples:

			A simple call that will return the list

			>>> data_file.truncatenans('Bitumen')

		"""

		if isinstance(variable, str):
			return [x for x in self.data[variable] if not np.isnan(x)]
		else:
			try:
				return [x for x in variable if not np.isnan(x)]
			except:
				raise Exception('Could not truncate NaNs!')

	def addcoord(self):
		"""
		Only use on DataFile classes containing GSLIB style gridded data.

		If x, y, or z coordinate column(s) do not exist they are created. If the created or current
		columns only have null values, they are populated based on the GridDef class pass to the
		DataFile class.

		Note:

			A griddef must be assigned to the DataFile class either at read in like here

			>>> data_file = gs.DataFile(flname='test.out', griddef=grid)

			Or later such it can be manually assigned such as here

			>>> data_file.griddef = gs.GridDef(gridstr=my_grid_str)

		"""

		# Make sure the griddef in the class is set
		if self.griddef is None:
			print('ERROR: GridDef is not defined within the DataFile class')
			return
		# Make sure the number of data makes sense with the griddef
		if len(self.data) % self.griddef.count() != 0:
			print('ERROR: The number of data found does not work with the number of cells')
			return
		# Check to make sure coordinate columns don't exist and if they do, force them to be loaded
		# to the DataFile class
		columns = list(self.data.columns)
		if 'x' in columns and self.x is None:
			print("ERROR: The column 'x' already exists and isn't saved to the DataFile class")
			return
		if 'y' in columns and self.y is None:
			print("ERROR: The column 'y' already exists and isn't saved to the DataFile class")
			return
		if 'z' in columns and self.z is None:
			print("ERROR: The column 'z' already exists and isn't saved to the DataFile class")
			return
		# Add coordinate columns as needed
		if self.x is None:
			self.data['x'] = float("nan")
			self.x = 'x'
		if self.y is None:
			self.data['y'] = float("nan")
			self.y = 'y'
		if self.z is None:
			self.data['z'] = float("nan")
			self.z = 'z'
		# Move the new coord columns to the start
		xyzcols = [self.x, self.y, self.z]
		columns = [columns[idx] for idx in range(len(columns)) if columns[idx] not in xyzcols]
		xyzcols.extend(list(columns))
		self.data = self.data[xyzcols]
		# Figure out what columns need to be populated
		reqcols = []
		for column in xyzcols:
			if self.data[column].isnull().values.sum() == len(self.data[column]):
				reqcols.append(True)
			else:
				reqcols.append(False)
		# Add the coordinates
		nreal = len(self.data) // self.griddef.count()
		x, y, z = self.griddef.get_coordinates()
		gridxyz = np.stack((x, y, z), axis=1)
		xyzreals = gridxyz.copy()
		for ireal in range(1, nreal):
			xyzreals = np.concatenate((xyzreals, gridxyz), axis=0)
		x = xyzreals[:, 0]
		y = xyzreals[:, 1]
		z = xyzreals[:, 2]
		if reqcols[0]:
			self.data.loc[:, self.x] = x
		if reqcols[1]:
			self.data.loc[:, self.y] = y
		if reqcols[2]:
			self.data.loc[:, self.z] = z

	def infergriddef(self, blksize=None, databuffer=5, nblk=None):
		"""
		Infer a grid definition with the specified dimensions to cover the set of data values. The
		function operates with two primary options:

			1. Provide a block size (node spacing), the function infers the required number
			   of blocks (grid nodes) to cover the data
			2. Provide the number of blocks, the function infers the required block size

		A data buffer may be used for expanding the grid beyond the data extents. Basic integer
		rounding is also used for attempting to provide a 'nice' grid in terms of the origin
		alignment.

		Parameters:
			blksize(float or 3-tuple): provides (xsiz, ysiz, zsiz). If blksize is not None,
				nblk must be None. Set zsiz None if the grid is 2-D. A float may also be provided,
				where xsiz = ysiz = zsiz = float is assumed.
			databuffer (float or 3-tuple): buffer between the data and the edge of the model,
				optionally for each direction
			nblk (int or 3-tuple): provides (nx, ny, nz). If blksize is not None,
				nblk must be None. Set nz to None or 1 if the grid is 2-D. An int may also be provided,
				where nx = ny = nz = int is assumed.

		Returns:
			griddef (GridDef): this function returns the grid definition object as well as
			assigns the griddef to the current `gs.DataFile`

		Note:
			this function assumes things are either 3D or 2D along the xy plane. If nx == 1
			or ny == 1, nonsense will result!

		Usage:

			First, import a datafile using gs.DataFile(), make sure to assign the correct columns
			to x, y and z:

			>>> datfl = gs.DataFile('test.dat',x='x',y='y',z='z')

			Now create the griddef from the data contained within the dataframe:

			>>> blksize = (100, 50, 1)
			>>> databuffer = (10, 25, 0) # buffer in the x, y and z directions
			>>> griddef = datfl.infergriddef(blksize, databuffer)

			Check by printing out the resulting griddef:

			>>> print(griddef)

		Examples:

			For 3D data, infergriddef() returns a 3D grid definition even if zsiz is given as None or 0 or 1:

			.. code-block:: python

				df3d = gs.ExampleData("point3d_ind_mv")
				a = df3d.infergriddef(blksize = [50,60,1])
				b = df3d.infergriddef(blksize = [50,60,None])
				c = df3d.infergriddef(blksize = [50,60,0])
				#a,b,c are returned as Pygeostat GridDef:
				#					20 135.0 50.0 
				#					19 1230.0 60.0 
				#					82 310.5 1.0

			For 3D data, nz given as None or 0 or 1 returns a 2D grid that covers the vertical extent of the 3D data:

			.. code-block:: python				

				d = df3d.infergriddef(nblk = [50,60,1])
				e = df3d.infergriddef(nblk = [50,60,None])
				f = df3d.infergriddef(nblk = [50,60,0])
				#d,e,f are returned as Pygeostat GridDef:
				#					50 119.8 19.6 
				#					60 1209.1 18.2 
				#					1 350.85 81.7

			Where xsiz = ysiz = zsiz, a float can also be provided, or where nx = ny = nz, an int can also be provided:

			.. code-block:: python

				df3d.infergriddef(blksize = 75)
				df3d.infergriddef(blksize = [75,75,75])#returns the same as its above line	

				df3d.infergriddef(nblk = 60)
				df3d.infergriddef(nblk = [60,60,60])#returns the same as its above line	

       		If data is 2-D, zsiz or nz must be provided as None. Otherwise it raise exception:

	   		.. code-block:: python

				df2d = gs.ExampleData("point2d_ind")
				df2d.infergriddef(nblk = [60, 60, None])
				df2d.infergriddef(blksize = [50,60,None])
		"""
		
		from .grid_definition import GridDef
		from ..datautils.utils import round_sigfig
		import math
		# Check the parameter inputs
		twod = False
		if nblk is None and blksize is None:
			raise ValueError(("ERROR: both nblk and blksize are specified."
							  "One is specified and the other is inferred!"))
		elif nblk is not None and blksize is not None:
			raise ValueError(("ERROR: both nblk and blksize are not specified."
							  "One is specified and the other is inferred!"))
		elif nblk is not None:
			# nblk is the constant, so check its input
			if isinstance(nblk, tuple) or isinstance(nblk, list):
				if len(nblk) != 3:
					raise ValueError(("ERROR: nblk should be an integer or a length 3 tuple!"))
				else:
					nx, ny, nz = nblk[0], nblk[1],  nblk[2]
			else:
				nx = nblk
				ny = nblk
				nz = nblk
			nx = int(nx); ny = int(ny) 
			if nz is not None: nz = int(nz)
			if nz is None or nz == 0:
				nz = 1
				twod = True
		else:
			# blksize is the constant, so check its input
			if isinstance(blksize, tuple) or isinstance(blksize, list):
				if len(blksize) != 3:
					raise ValueError(("ERROR: blksize should be a float or a length 3 tuple!"))
				else:
					xsiz, ysiz, zsiz = blksize[0], blksize[1],  blksize[2]
			else:
				xsiz = blksize
				ysiz = blksize
				zsiz = blksize
			if zsiz is None or zsiz == 0:
				twod = True
				zsiz = 1.0
		# Make sure either 1 or 3 values passed to the buffer.
		if isinstance(databuffer, tuple) or isinstance(databuffer, list):
			if len(databuffer) != 3:
				raise ValueError(("ERROR: Ensure a 3 buffers are passed if specifying"
								  "a buffer for each direction"))
			else:
				bx, by, bz = databuffer
		else:
			bx = databuffer
			by = databuffer
			bz = databuffer
		# Very basic check on the data inputs
		if self.x is None:
			raise ValueError("ERROR: the x column must be saved to the DataFile class")
		if self.y is None:
			raise ValueError("ERROR: the y column must be saved to the DataFile class")
		if not twod:
			if self.z is None:
				raise ValueError("ERROR: the z column must be saved to the DataFile class unless 2-D")
		# get the min's and max's
		xmin, xmax = (min(self.data[self.x]) - bx), (max(self.data[self.x]) + bx)
		ymin, ymax = (min(self.data[self.y]) - by), (max(self.data[self.y]) + by)
		if not twod or self.z is not None:
			zmin, zmax = (min(self.data[self.z]) - bz), (max(self.data[self.z]) + bz)
		else:
			zmin, zmax = (0.0, 1.0)
		# round mins down/up to the nearest integer, likely providing a nicer
		# origin/division
		xmin -= xmin % 1.0
		ymin -= ymin % 1.0
		zmin -= zmin % 1.0
		xmax += xmax % 1.0
		ymax += ymax % 1.0
		zmax += zmax % 1.0
		if nblk is not None:
			# Infer the block sizes - note that this should be divided by nx, not nx-1
			# Since it's the distance of the grid extents (not the min/max centroid)
			xsiz = (xmax - xmin) / nx
			ysiz = (ymax - ymin) / ny
			if not twod or self.z is not None:
				zsiz = (zmax - zmin) / nz
			else:
				zsiz = 1.0
			sizes = [xsiz, ysiz, zsiz]
			xsiz, ysiz, zsiz = [float(round_sigfig(v, 3)) for v in sizes]
		else:
			# Infer the number of blocks
			nx = math.ceil((xmax - xmin) / xsiz)
			ny = math.ceil((ymax - ymin) / ysiz)
			if not twod or self.z is not None:
				nz = math.ceil((zmax - zmin) / zsiz)
			else:
				nz = 1

		griddef = GridDef(grid_arr=[nx, xmin + 0.5 * xsiz, xsiz,
								   ny, ymin + 0.5 * ysiz, ysiz,
								   nz, zmin + 0.5 * zsiz, zsiz])

		self.griddef = griddef
		return griddef

	def spacing(self, n_nearest, var=None, inplace=True, dh=None, x=None, y=None):
		'''
		Calculates data spacing in the xy plane, based on the average distance to the nearest
		n_nearest neighbours. The x, y coordinates of 3-D data may be provided in combination
		with a dh (drill hole or well), in which case the mean x, y of each dh is calculated
		before performing the calculation. If a dh is not provided in combination with 3-D xy's,
		then calculation is applied to all data and may create memory issues if greater than
		~5000-10000 records are provided. A var specifier allows for the calculation to only
		applied where var is not NaN.

		If ``inplace==True``:

			The output is concatenated as a 'Data Spacing ({Parameters['plotting.unit']})' column
			if ``inplace=False`` (or 'Data Spacing' if Parameters['plotting.unit'] is None). If var is
			used, then the calculation is only performed where DataFile[var] is not NaN, and
			the output is concatenated as '{var} Data Spacing ({Parameters['plotting.unit']})'.

		If ``inplace==False``:

			The funciton returns dspace as a numpy array if dspace.shape[0] is equal to
			DataFile.shape[0], meaning that dh and var functionality was not used, or did
			not lead to differences in the length of dspace and DataFile (so that the x and
			y in DataFile can be used for plotting dspace in map view).
			The function returns a tuple of the form (dspace, dh, x, y), if dh is not None and
			dspace.shape[0] is not equal to DataFile.shape[0]. The function returns a tuple of
			the form (dspace, x, y) if dh is None and and var is not None and dspace.shape[0]
			is not equal to DataFile.shape[0].

		Parameters:
			n_nearest (int): number of nearest neighbours to consider in data spacing
				calculation
			var (str): variable for calculating data spacing, where the calculation is only
				applied to locations where var is not NaN.  If None, the calculation is
				to all locations.
			inplace (bool): if True, the output data spacing is concatenated
			dh (str): dh name, which can override self.dh
			x (str): x coordinate name, which can override self.x
			y (str): y coordinate name, which can override self.y

		Examples:

			Calculate data spacing without consideration of underlying variables, based on the
			nearest 8 neighbours.

			>>> dat.spacing(8)

			Output as a numpy array rather than concatenating a column:

			>>> dspace = dat.spacing(8, inplace=False):

			Only consider values where Au is non NaN for the calculation:

			>>> (dspace, x, y) = dat.spacing(8, inplace=False, var=Au)

		'''
		# Check the data input
		if self.data is None:
			raise ValueError('DataFile.data must be initialized!')
		if self.data.ndim < 2:
			raise ValueError(('DataFile.data.ndim must be greater than 1, given than an x and y'
							  'columns are required!'))
		# Check the coordinate inputs
		if x is None:
			if self.x is None:
				raise Exception('x must be provided either as DataFile.x or the kwarg!')
			xc = self[self.x].values
		else:
			xc = self[x].values
		if y is None:
			if self.y is None:
				raise Exception('y must be provided either as DataFile.y or the kwarg!')
			yc = self[self.y].values
		else:
			yc = self[y].values
		# Check the drill hole
		if dh is None:
			if self.dh is None:
				dhc = np.arange(self.shape[0])
			else:
				dhc = self[self.dh].values
		else:
			dhc = self[dh].values
		# Check the var
		if var is not None:
			if var not in self.columns:
				raise KeyError('{} is not in DataFile.data!'.format(var))
			vidx = np.logical_not(np.isnan(self[var].values))
			xc = xc[vidx]
			yc = yc[vidx]
			dhc = dhc[vidx]
		# Record this dht for output indexing
		dhr = dhc
		# Calculate the x and y as the average of each dh
		dhu = np.unique(dhc)
		if dhu.shape[0] == dhc.shape[0]:
			xu = xc
			yu = yc
			dhu = dhc
		else:
			xu = np.zeros(dhu.shape[0])
			yu = np.zeros(dhu.shape[0])
			for i, dhut in enumerate(dhu):
				idx = dhc == dhut
				xu[i] = np.mean(xc[idx])
				yu[i] = np.mean(yc[idx])
		# Record the c's for future calculations
		xc = xu
		yc = yu
		dhc = dhu
		# Check the n_nearest now that we've reduced down to the calculated data
		n = xc.shape[0]
		if not isinstance(n_nearest, int):
			raise ValueError('n_nearest must be an integer!')
		elif n_nearest < 1:
			raise ValueError('n_nearest must be larger than 0!')
		elif n_nearest > n:
			raise ValueError(('n_nearest must be less than or equal to the '
							  'number of values (after dh/var filtering)!'))
		if n > 5000:
			print(('WARNING: current implementation of function is likely too memory intensive'
				   'for greater than 5000 data'))
		# Repeat the vector n times
		xc = np.tile(xc, (n, 1))
		yc = np.tile(yc, (n, 1))
		# Calculate the squared distance between all points
		# Distance matrix results (0s on the diagonal)
		dspace = np.square(np.subtract(xc, xc.T))
		dspace = np.add(dspace, np.square(np.subtract(yc, yc.T)))
		dspace = np.sort(dspace, axis=0)
		# Calcate the average distance to the n_nearest data
		dspace = np.sqrt(dspace[1:n_nearest+1, :])
		dspace = np.mean(dspace, axis=0)
		if inplace:
			# Build the output array
			dspace1 = np.zeros(dhr.shape[0])
			for i, dhut in enumerate(dhu):
				idx = dhr == dhut
				dspace1[idx] = dspace[i]
			if var is not None:
				dspace = np.zeros(self.shape[0])
				dspace[np.where(vidx)] = dspace1
				dspace[np.where(np.logical_not(vidx))] = np.nan
			else:
				dspace = dspace1
			# Figure out the data spacing name and assign
			if var is not None:
				name = var+' '
			else:
				name = ''
			name = name+'Data Spacing'
			if Parameters['plotting.unit'] is not None:
				name = name+' ({})'.format(Parameters['plotting.unit'])
			self[name] = dspace
		else:
			# Determine the form of the output
			if self.shape[0] == n:
				return dspace
			elif dh is None:
				return (dspace, xc, yc)
			else:
				return (dspace, dhc, xc, yc)


class DictFile:

	'Class containing dictionary file information'

	def __init__(self, flname=None, readfl=False, dictionary={}):
		self.flname = flname
		self.dict = dictionary
		if readfl:
			self.read_dict()

	def read_dict(self):
		'Read dictionary information from file'
		reader = csv.reader(open(self.flname, 'r'), delimiter=',')
		for item in reader:
			if item:
				self.dict[item[0]] = item[1]

	def write_dict(self):
		'Write dictionary information to csv style dictionary'
		writerfile = open(self.flname, 'w')
		writer = csv.writer(writerfile, delimiter=',', lineterminator='\n')
		for key, value in self.dict.items():
			writer.writerow([key, value])
		writerfile.close()



def ExampleData(testfile, griddef=None, **kwargs):
	"""
	Get an example pygeostat DataFile

	Parameters:
		testfile (str): one of the available pygeostat test files, listed below

	Test files available in pygeostat include:

		* "point2d_ind": 2d indicator dataset
		* "point2d_surf": 2d point dataset sampling a surface
		* "grid2d_surf": 'Thickness' from 'point2d_surf' interpolated on the grid
		* "point3d_ind_mv": 3d multivariate and indicator dataset
		* "oilsands": 3D Oil sands data set
		* "accuracy_plot": Simulated realizations to test accuracy plot
		* "location_plot": 2D data set to test location plot
		* "3d_grid": 3D gridded data set
		* "point2d_mv" : 2D multivariate data set
		* "cluster": GSLIB datafile (data with declustering weights)
		* "97data": GSLIB datafile (the first 97 rows of cluster datafile)
		* "data": GSLIB datafile (2D data set of primary and secondary variable)
		* "parta": GSLIB datafile (small 2D dataset part A)
		* "partb": GSLIB datafile (small 2D dataset part B)
		* "partc": GSLIB datafile (small 2D dataset part C)
		* "true": GSLIB datafile (Primary secondary data pairs)
		* "ydata": GSLIB datafile (2D spatial seondary data with some primary data)

	"""

	_testfiles = {
		'point2d_ind': "point2d_ind.dat",
		'point2d_surf': "point2d_surf.dat",
		'grid2d_surf': "grid2d_surf.dat",
		'point3d_ind_mv': "point3d_ind_mv.dat",
		'3d_grid' : '3d_grid.out',
		'accuracy_plot': 'accuracy_plot.dat',
		'location_plot': 'location_plot.dat',
		'oilsands': 'oilsands.dat',
		'point2d_mv': 'point2d_mv.dat',
		'3d_correlation': '3d_correlation.dat',
		'3d_estimate': '3d_estimate.out',
		'experimental_variogram': 'varcalc.out',
		'variogram_model': 'varmodel.out',
		'reservoir_boundary': 'reservoir_boundary.dat',
		'reservoir_data': 'reservoir_data.dat',
		'reservoir_surface': 'reservoir_surface.dat',
		'cluster': 'cluster.dat',
		'97data': '97data.dat',
		'data': 'data.dat',
		'parta': 'parta.dat',
		'partb': 'partb.dat',
		'partc': 'partc.dat',
		'true': 'true.dat',
		'ydata': 'ydata.dat',
	}
	if testfile not in _testfiles:
		raise ValueError("Invalida file name. Choose one of {}".format(list(_testfiles.keys())))

	data_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), r'example_data'))

	return DataFile(os.path.join(data_dir, _testfiles[testfile]), griddef=griddef, **kwargs)


