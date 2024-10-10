#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function, unicode_literals

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import six
import os
import sys

if sys.version_info >= (3, 10):
	from collections.abc import MutableMapping
else:
	from collections import MutableMapping
	
import warnings

# Validation functions for parameters
#----------------------------------------------------------------------------------------------------------
def _validate_string(s, accept_none = False):
	"""
	 A validation method to convert input s to string or raise error if it is not convertable
	"""
	if s is None and accept_none :
		return None
	try:
		if isinstance(s,list):
			return [six.text_type(item) for item in s]
		elif isinstance(s,dict):
			return dict((six.text_type(key), six.text_type(value)) for key, value in s.items())
		else:
			return six.text_type(s)
	except ValueError:
		raise ValueError('Could not convert "%s" to string' % s)

def _validate_string_or_None(s):
	return _validate_string(s, accept_none=True)

def _validate_float(s, accept_none = False):
	"""
	A validation method to convert input s to float or raise error if it is not convertable
	"""

	if s is None and accept_none :
		return None
	try:
		return float(s)
	except ValueError:
		raise ValueError('Could not convert "%s" to float' % s)

def _validate_float_or_None(s):
	return _validate_float(s, accept_none=True)

def _validate_int(s, accept_none = False):
	"""
	A validation method to convert input s to int or raise error if it is not convertable
	"""
	if s is None and accept_none :
		return None
	try:
		return int(s)
	except ValueError:
		raise ValueError('Could not convert "%s" to int' % s)

def _validate_int_or_None(s):
	return _validate_int(s, accept_none=True)

def _validate_color(s):

	"""
	A validation method to check s is a valid matplotlib color
	"""

	from matplotlib.colors import is_color_like

	try:
		if s.lower() == 'none':
			return 'none'
	except AttributeError:
		pass

	if isinstance(s, six.string_types):
		if len(s) == 6 or len(s) == 8:
			stmp = '#' + s
			if is_color_like(stmp):
				return stmp

	if is_color_like(s):
		return s

def _validate_list(s, accept_none=False):
	"""
	A validation method to convert input s to a list or raise error if it is not convertable
	"""
	if s is None and accept_none:
		return None
	try:
		return list(s)
	except ValueError:
		raise ValueError('Could not convert input to list')

def _validate_list_or_None(s):
	return _validate_list(s, accept_none=True)

def _validate_list_int(s, accept_none=False):
	"""
	A validation method to convert input s to int or list and raise error if it is not convertable
	"""
	if s is None and accept_none:
		return None
	try:
		return int(s)
	except (ValueError, TypeError):
		try:
			return list(s)
		except ValueError:
			raise ValueError('Could not convert input to int or list')

def _validate_list_int_or_None(s):
	return _validate_list_int(s, accept_none=True)

def _validate_bool(s):
	"""
	A validation method to convert input s to boolean or raise error if it is not convertable
	"""
	try:
		return bool(s)
	except ValueError:
		raise ValueError('Could not convert input to boolean')

def _validate_dict_or_string(s):
	"""
	A validation method to convert input s to string or dictionary and raise error if it is not convertable
	"""
	if isinstance(s, str) or isinstance(s, dict):
		return s
	try:
		return str(s)
	except ValueError:
		raise ValueError('Could not convert input to dict or string')

def _validate_dict(s, accept_none=False):
	"""
	A validation method to check if the input s is a dictionary otherwise raise error if it is not convertable
	"""
	if s is None and accept_none:
		return None
	if isinstance(s, dict):
		return s
	else:
		raise ValueError('{} is not a dictionary!'.format(s))

def _validate_dict_or_None(s):
	return _validate_dict(s, accept_none=True)

def _validate_kde_or_color(s):
	'''
	A validation method to return a valid color arg if not KDE
	'''

	from matplotlib.colors import is_color_like

	try:
		if s.lower() == 'kde':
			return 'kde'
	except AttributeError:
		pass
	try:
		if s.lower() == 'none':
			return 'none'
	except AttributeError:
		pass

	if isinstance(s, six.string_types):
		if len(s) == 6 or len(s) == 8:
			stmp = '#' + s
			if is_color_like(stmp):
				return stmp

	if is_color_like(s):
		return s

def _validate_griddef_or_None(s):
	"""
	A validation method to convert input s to pygeostat grid def and raise an error if not convertable.
	"""
	import pygeostat as gs
	if s is None:
		return None
	elif isinstance(s, str):
		try:
			return gs.GridDef(s)
		except AssertionError:
			raise ValueError('{} should be a pygeostat.GridDef object or griddef string!'.format(s))
	elif isinstance(s, gs.GridDef):
		return s
	else:
		raise ValueError('{} should be a pygeostat.GridDef object or griddef string!'.format(s))

#----------------------------------------------------------------------------------------------------------

# Pygeostat configurations file
_params_cfg_file = os.path.join(os.path.expanduser('~'),r".Pygeostat\Parameters.json")

# Writing/loading parameters to a json file
#----------------------------------------------------------------------------------------------------------
def _write_pardict(filename, dictionary, section, description=None):
	import json
	from .data.grid_definition import GridDef

	if os.path.isfile(filename):
		with open(filename, 'r') as cfg_file:
			dict_main = dict(json.load(cfg_file))
	else:
		dict_main = {}
	if section in dict_main.keys():
		dict_main.pop(section, None)
	
	section_dict={}
	for key in list(sorted(dictionary.keys())):
		if 'griddef' in str(key).lower() and dictionary[key] is not None:
			value = str(dictionary[key])
		else:
			value = dictionary[key]
		section_dict.update({key: {'Value': value, 'Description':None if description is None else description[key]}})
	
	dict_main.update({section:section_dict})

	with open(filename, 'w') as cfg_file:
		json.dump(dict_main, cfg_file)

def _parse_pardict(filename, section):
	import json
	if os.path.isfile(filename):
		try:
			with open(filename, 'r') as fp:
				dict_main = dict(json.load(fp))
		except:
			print('Unable to open the provided parameter file: {}'.format(filename))
			return {}
	else:
		print('Unable to load the system defaults as it does no exist! Use gs.Programs.set_systems_defaults to create one.')
		return {}
	# Allowing no values
	if section not in dict_main.keys():
		return {}
		
	new_dict = {}
	input_dict = dict_main[section]
	for key in input_dict:
		new_dict.update({key: input_dict[key]['Value']})

	return new_dict

# Parameter list
#---------------------------------------------------------------------------------------------------
default_parameters = {
	'config.verbose': [True, _validate_bool],
	'config.autoload.parameters': [True, _validate_bool],
	'config.autoload.plot_style': [False, _validate_bool],
	'config.use_shipped_executable': [True, _validate_bool],
	'config.nprocess': [4, _validate_int_or_None],
	'config.getpar': [False, _validate_bool],
	'config.ignore_mpl_warnings': [True, _validate_bool],
	'plotting.unit': ['m', _validate_string_or_None],
	'plotting.xname': ['Easting', _validate_string_or_None],
	'plotting.yname': ['Northing', _validate_string_or_None],
	'plotting.zname': ['Elevation', _validate_string_or_None],
	'plotting.xabbrev': ['E', _validate_string_or_None],
	'plotting.yabbrev': ['N', _validate_string_or_None],
	'plotting.zabbrev': ['Elev', _validate_string_or_None],
	'plotting.gammasize': [3.0, _validate_float],
	'plotting.lagname': ['Lag Distance', _validate_string_or_None],
	'plotting.rotateticks': [[0, 0], _validate_list_or_None],
	'plotting.nticks': [None, _validate_list_int_or_None],
	'plotting.grid': [False, _validate_bool],
	'plotting.axis_xy': [False, _validate_bool],
	'plotting.axis_xy_spatial': [False, _validate_bool],
	'plotting.stat_fontsize': [None, _validate_float_or_None],
	'plotting.stat_ha': ['right', _validate_string],
	'plotting.stat_linespacing': [1.0, _validate_float],
	'plotting.roundstats': [True, _validate_bool],
	'plotting.sigfigs': [2, _validate_int],
	'plotting.cmap': ['viridis', _validate_string],
	'plotting.cmap_cat': ['tab20', _validate_dict_or_string],
	'plotting.assumecat': [11, _validate_int],
	'plotting.log_lowerval': [0.0001, _validate_float],
	'plotting.histogram_plot.facecolor': ['.9', _validate_color],
	'plotting.histogram_plot.edgecolor': ['k', _validate_color],
	'plotting.histogram_plot.edgeweight': [None, _validate_float_or_None],
	'plotting.histogram_plot.cdfcolor': ['.5', _validate_color],
	'plotting.histogram_plot.histbins': [15, _validate_int_or_None],
	'plotting.histogram_plot.stat_blk': ['all', _validate_string],
	'plotting.histogram_plot.stat_xy': [(0.95, 0.95), _validate_list],
	'plotting.histogram_plot.stat_xy_cdf': [(0.95, 0.05), _validate_list],
	'plotting.histogram_plot_simulation.refclr': ['C3', _validate_color],
	'plotting.histogram_plot_simulation.simclr': ['0.2', _validate_color],
	'plotting.histogram_plot_simulation.alpha': ['0.5', _validate_float_or_None],
	'plotting.location_plot.s': [20, _validate_float_or_None],
	'plotting.location_plot.c': ['.4', _validate_color],
	'plotting.location_plot.lw': ['5', _validate_float_or_None],
	'plotting.scatter_plot.s': [None, _validate_float_or_None],
	'plotting.scatter_plot.c': ['kde', _validate_kde_or_color],
	'plotting.scatter_plot.cmap': ['viridis', _validate_string],
	'plotting.scatter_plot.alpha': [1.0, _validate_float_or_None],
	'plotting.scatter_plot.stat_blk': ['pearson', _validate_string],
	'plotting.scatter_plot.stat_xy': [(0.95, 0.05), _validate_list],
	'plotting.variogram_plot.color': ['.5', _validate_color],
	'plotting.variogram_plot.ms': [6.0, _validate_float],
	'plotting.variogram_plot_simulation.refclr': ['C3', _validate_color],
	'plotting.variogram_plot_simulation.simclr': ['0.2', _validate_color],
	'plotting.variogram_plot_simulation.alpha': ['0.5', _validate_float_or_None],
	'plotting.vplot.colors': [['C0', 'C1', 'C2'], _validate_string],
	'data.tmin': [-998.0, _validate_float_or_None],
	'data.null': [-999.0, _validate_float_or_None],
	'data.legacy_null': [(-999, -998, -99, -98), _validate_list_or_None],
	'data.fix_legacy_null': [False, _validate_bool],
	'data.null_vtk': [False, _validate_float_or_None],
	'data.griddef': [None, _validate_griddef_or_None],
	'data.nreal': [None, _validate_int_or_None],
	'data.dh': [None, _validate_string_or_None],
	'data.ifrom': [None, _validate_string_or_None],
	'data.ito': [None, _validate_string_or_None],
	'data.x': [None, _validate_string_or_None],
	'data.y': [None, _validate_string_or_None],
	'data.z': [None, _validate_string_or_None],
	'data.cat': [None, _validate_string_or_None],
	'data.catdict': [None, _validate_dict_or_None],
	'data.weights': [None, _validate_string_or_None],
	'data.io.pandas_engine': ['python', _validate_string],
	'data.write.python_floatfmt': ['%.6f', _validate_string],
	'data.write_vtk.vdtype': ['float64', _validate_string],
	'data.write_vtk.cdtype': ['float64', _validate_string],
}

# Description of the parameters
default_parameter_descriptions = {
	'config.verbose': ('If `true`, mentions when the `gsParams` are updated automatically '
					   'when pygeostat is loaded. Useful for interactive sessions'),
	'config.autoload.parameters': ('If `True` the `gsParams` configurations found in\n'
								 '``%USER%/.gsParams`` are parsed and loaded when pygeostat loads'),
	'config.autoload.plot_style': ('If `True` the `gsPlotStyle` configurations found in\n'
									'``%USER%/.gsParams`` are parsed and loaded when pygeostat loads'),
	'config.use_shipped_executable': ('If `True` adds the shipped executable files to the path environmental variable'),
	'config.nprocess': ('The number of parallel processes to run functions that provide that '
						'functionality'),
	'config.getpar': ('If True, getpar=True when calling pygeostat.Program() with getpar=None'),
	'config.ignore_mpl_warnings': 'If True, matplotlib warnings will be ignored!',
	'plotting.unit': ('The unit (str) that appears in spatial/map/variogram plots '
					  'e.g., m or ft. Use None if no unit should be displayed.'),
	'plotting.xname': ('Name (str) that is used for the x coordinate label for '
					   'map/cross-section plots, e.g., Easting or X'),
	'plotting.yname': ('Name (str) that is used for the y coordinate label for '
					   'map/cross-section plots, e.g., Northing or Y'),
	'plotting.zname': ('Name (str) that is used for the z coordinate label for '
					   'map/cross-section plots, e.g., Elevation or Z'),
	'plotting.xabbrev': ('Abbreviated name (str) that is used for the x coordinate '
						 'label for map/cross-section plots, e.g., Easting or X'),
	'plotting.yabbrev': ('Abbreviated name (str) that is used for the y coordinate '
						 'label for map/cross-section plots, e.g., Northing or Y'),
	'plotting.zabbrev': ('Abbreviated name (str) that is used for the z coordinate '
						 'label for map/cross-section plots, e.g., Elevation or Z'),
	'plotting.gammasize': ('The size of gamma symbols in variogram plots are '
						   'Matplotlib.rcParams["font.size"] multiplied with this value.'),
	'plotting.lagname': ('Name (str) to use on the x-axis of variogram plots, '
						 'e.g., Lag or h'),
	'plotting.rotateticks': ('If None, the pygeostat tickoverlap attempts to optimize tick, '
							 'label angles to minimize overlap in the absence of a kwarg., '
							 'If [xangle, yangle] (e.g., [0, 0], then these angles are used '
							 'in the absence of a kwarg.'),
	'plotting.nticks': ('Specify the target number of ticks for plotting, as an integer, or '
						'tuple/list of integers specifying the number of ticks on the x and y '
						'axis, respectively'),
	'plotting.grid': ('If True, a grid is plotted in pygeostat ploting function unless, '
					  'overridden by the associated kwarg'),
	'plotting.axis_xy': ('Converts plotting axes to GSLIB-style visibility (only left and '
						 'bottom visible) if axis_xy is True. See plotting.axis_xy_spatial '
						 'for the default setting of functions such as slice_plot and location_plot.'),
	'plotting.axis_xy_spatial': ('axis_xy setting that is specific to spatial plotting functions '
								 'such as location_plot and slice_plot. Provided since many users '
								 'prefer axis_xy to be applied to all plots (via '
								 'plotting.axis_xy) other than spatial plots. '),
	'plotting.stat_fontsize': ('Font size of statistics in pygeosat plots. If None, then '
							   'the font size matches rcParams[font.size]. If a fraction, '
							   'then the font size is rcParams[font.size]*stat_fontsize. '
							   'If greater than or equal to 1, then stat_fontsize is the '
							   'font size.'),
	'plotting.stat_linespacing': ('The line spacing to use for statistics blocks in plots. '
								  'Default is 1.0'),
	'plotting.stat_ha': ('The horizontal alignment of statistics in plots, which should be one of '
						 '"right", "left" or "center"'),
	'plotting.roundstats': ('Display a set number of significant figures (False) or digits (True) '
							'for statistics in plots.'),
	'plotting.sigfigs': ('Significant figures (roundstats=False) or digits (roundstats=True) '
						 'for statistics in plots.'),
	'plotting.cmap': ('Matplotlib colormap object or a registered Matplotlib colormap, which is '
					  'used for continuous variables'),
	'plotting.cmap_cat': ('Matplotlib colormap object or a registered Matplotlib colormap, which '
						  'is used for categorical variables'),
	'plotting.assumecat': ('When executing plotting functions, data is assumed to be categorical '
						   'for the purposes of selecting colormaps, if less than this number '
						   'of unique values are found'),
	'plotting.log_lowerval': ('For a log axis, specify the default lower value that replaces a '
							  '0.0 or negative valued data. Applies to both axes.'),
	'plotting.histogram_plot.facecolor': 'Color of the histogram faces for the histogram_plot function.',
	'plotting.histogram_plot.edgecolor': 'Color of the histogram edges for the histogram_plot function.',
	'plotting.histogram_plot.edgeweight': '`lw` of the histogram edges for the histogram_plot function.',
	'plotting.histogram_plot.cdfcolor': 'Color of the CDF for the histogram_plot function.',
	'plotting.histogram_plot.histbins': 'Number of bins for a histogram (not CDF) calculation.',
	'plotting.histogram_plot.stat_blk': ('Default stat_blk setting, which is either a string code, such '
								  'as all or minimal, or a list of valid statistic names'),
	'plotting.histogram_plot.stat_xy': ('Location of the histogram stats location as a tuple '
								 'of (xloc, yloc), where the locations should be between 0 and 1.'),
	'plotting.histogram_plot.stat_xy_cdf': ('Location of the CDF stats location as a tuple of  '
									 '(xloc, yloc), where the locations should be 0 to 1.'),
	'plotting.histogram_plot_simulation.refclr': 'Color of the reference CDF',
	'plotting.histogram_plot_simulation.simclr': 'Color of the realization CDFs',
	'plotting.histogram_plot_simulation.alpha': 'Transparency of the realizations CDFs',
	'plotting.location_plot.s': ('Size of scatter in location_plot (and slice_plot), which is based on '
						  'Matplotlib.rcParams if None'),
	'plotting.location_plot.c': ('Color of scatter in location_plot (and slice_plot), which is a valid Matplotlib '
						  'color specifier'),
	'plotting.location_plot.lw': ('Line width for location plot if the orient is xz or yz and the drill hole id is provided/inferred'),
	'plotting.scatter_plot.s': ('Size of scatter in scatter_plot, which is based on Matplotlib.rcParams if '
						   'None'),
	'plotting.scatter_plot.c': ('Color of scatter in scatter_plot, which is either a valid Matplotlib '
						   'color specifier, or the "KDE" string. KDE leads to calculation of '
						   'the kernel density estimate at each scatter location'),
	'plotting.scatter_plot.cmap': ('Matplotlib colormap object or a registered Matplotlib colormap '
							  'which overrides plotting.cmap'),
	'plotting.scatter_plot.alpha': ('Alpha transparency of scatter in scatter_plot, which should be '
							   'between 0 and 1.'),
	'plotting.scatter_plot.stat_blk': ('Statistics that are plotted in scatter_plot, which should be '
								  'either "all" or a list that may contains the strings: '
								  '["count", "pearson", "spearmanr"]'),
	'plotting.scatter_plot.stat_xy': 'A 2-list that provides the x/y location of the statistics block',
	'plotting.variogram_plot.color': 'Color for variogram_plot ',
	'plotting.variogram_plot.ms': 'Marker size for variogram_plot (e.g. experimental dot size) ',
	'plotting.variogram_plot_simulation.refclr': 'Color of the reference variogram',
	'plotting.variogram_plot_simulation.simclr': 'Color of the realization variograms',
	'plotting.variogram_plot_simulation.alpha': 'Transparency of the realizations variograms',
	'plotting.vplot.colors': ('Colors of the major, minor and vertical direction variograms, '
							  'which is used by the Variogram.plot function. May be a 3-list, '
							  'or a color pallete in gs.avail_palettes'),
	'data.tmin': ('Values less than this number (float or int) are assigned NaN on '
				  'import by all io read functions, e.g., -98 or -998. Using -1.0e21 '
				  'or False disables this functionality.'),
	'data.null': ('NaN values are assigned this number (float or int) by io write '
				  'functions (excluding write_vtk), e.g., -99, -999, or -999.0. Using '
				  'None disables the functionality, but will likely lead to issues '
				  'in the case of subsequent GSLIB operations.'),
	'data.null_vtk': ('NaN values are assigned this number (float or int) by the '
					  'write_vtk function, e.g., -99, -999, or -999.0. Using None '
					  'disables the functionality, which is recommended since NaN '
					  'values are handled by Paraview.'),
	'data.legacy_null': ('Nulls used by default in some GSLIB programs that can be fixed if '
						 'encountered'),
	'data.fix_legacy_null': ('Automatically fix legacy null values if encountered and '
							 'trimming null values'),
	'data.griddef': ('When initializing a DataFile, this will be used as '
					 'DataFile.GridDef if GridDef.count() matches DataFile.shape[0]. '
					 'A pygeostat.GridDef object or valid gridstr/gridarr may be used '
					 'for intitialization.'),
	'data.nreal': ('The number of realizations for modeling (once implemented) '
				   'and model checking.'),
	'data.dh': ('When initializing a DataFile, this str will be used as '
				'DataFile.dh if it is in DataFile.columns, e.g., WellID'),
	'data.ifrom': ('When initializing a DataFile, this str will be used as '
				   'DataFile.ifrom if it is in DataFile.columns, e.g., From'),
	'data.ito': ('When initializing a DataFile, this str will be used as '
				 'DataFile.ito if it is in DataFile.columns, e.g., To'),
	'data.x': ('When initializing a DataFile, this st will be used as '
			   'DataFile.x if it is in DataFile.columns, e.g., Easting'),
	'data.y': ('When initializing a DataFile, this str will be used as '
			   'DataFile.y if it is in DataFile.columns, e.g., Northing'),
	'data.z': ('When initializing a DataFile, this str will be used as '
			   'DataFile.z if it is in DataFile.columns, e.g., Elevation'),
	'data.cat': ('When initializing a DataFile, this str will be used as '
				 'DataFile.cat if it is in DataFile.columns, e.g., Facies'),
	'data.catdict': ('When initializing a DataFile, this dictionary will be used as '
					 'DataFile.catdict (if catdict.keys() matches DataFile.cat codes). This '
					 'dictionary is formaatted as {catcodes: catnames}.'),
	'data.weights': ('When initializing a DataFile, this str or list will be used as '
				 'DataFile.wts if wts is in DataFile.columns, e.g., DeclusteringWeight'),
	'data.io.pandas_engine': 'Default engine to use for pandas read, either `python` or `c`',
	'data.write.python_floatfmt': ('A format string to use for `gslib` fltypes (written with '
								   'pandas)'),
	'data.write_vtk.vdtype': ('Precision of variables when writing to VTK. Must be a valid '
							  'numpy specifier, such as "float32", "float64", "int8", etc.'),
	'data.write_vtk.cdtype': ('Precision of coordinates when writing to VTK. Must be a valid '
							  'numpy specifier, such as "float32", "float64", etc.'),
}

#----------------------------------------------------------------------------------------------------------

# Plotstyles
#----------------------------------------------------------------------------------------------------------
import matplotlib as mpl

nochanges = {}

darkcontent = {'axes.axisbelow': True, 'axes.edgecolor': 'white', 'axes.facecolor': 'darkgray',
			   'axes.grid': False, 'axes.labelcolor': 'white', 'axes.labelsize': 8,
			   'axes.linewidth': 0.5, 'axes.titlesize': 12,
			   'figure.figsize': (5, 5), 'figure.subplot.bottom': 0.125,
			   'figure.facecolor': 'black', 'figure.edgecolor': 'black',
			   'font.family': 'Calibri', 'font.size': 8, 'font.weight': 400, 'text.color': 'white',
			   'grid.color': 'lightgray', 'grid.linestyle': '-', 'grid.linewidth': 0.5,
			   'legend.fontsize': 8.0, 'legend.frameon': False, 'legend.numpoints': 1,
			   'legend.scatterpoints': 1,
			   'lines.linewidth': 1.5, 'lines.markeredgewidth': 0, 'lines.markersize': 7,
			   'lines.solid_capstyle': 'round', 'patch.linewidth': 0.3,
			   'pdf.fonttype': 3, 'ps.fonttype': 3, 'ps.useafm': True,
			   'xtick.color': 'white', 'xtick.direction': 'out', 'xtick.labelsize': 8,
			   'xtick.major.pad': 3, 'xtick.major.width': 1, 'xtick.major.size': 0,
			   'xtick.minor.width': 0.5, 'xtick.minor.size': 0, 'ytick.color': 'white',
			   'ytick.direction': 'out', 'ytick.labelsize': 8, 'ytick.major.pad': 3,
			   'ytick.major.width': 1, 'ytick.major.size': 0, 'ytick.minor.width': 0.5,
			   'ytick.minor.size': 0}

ccgpaper = {'axes.axisbelow': True, 'axes.edgecolor': 'black', 'axes.facecolor': 'white',
			'axes.grid': False, 'axes.labelcolor': 'black', 'axes.labelsize': 8,
			'axes.linewidth': 0.5, 'axes.titlesize': 10,
			'figure.figsize': (8, 8), 'figure.subplot.bottom': 0.125, 'figure.facecolor': 'white',
			'figure.edgecolor': 'black',
			'font.family': 'Calibri', 'font.size': 12, 'font.weight': 400, 'text.color': 'black',
			'grid.color': 'lightgray', 'grid.linestyle': '-', 'grid.linewidth': 0.5,
			'legend.fontsize': 8.0, 'legend.frameon': False, 'legend.numpoints': 1,
			'legend.scatterpoints': 1,
			'lines.linewidth': 1.5, 'lines.markeredgewidth': 0, 'lines.markersize': 7,
			'lines.solid_capstyle': 'round', 'patch.linewidth': 0.3,
			'pdf.fonttype': 3, 'ps.fonttype': 3, 'ps.useafm': True,
			'xtick.color': 'black', 'xtick.direction': 'out', 'xtick.labelsize': 11,
			'xtick.major.pad': 3, 'xtick.major.width': 1, 'xtick.major.size': 0,
			'xtick.minor.width': 0.5, 'xtick.minor.size': 0, 'ytick.color': 'black',
			'ytick.direction': 'out', 'ytick.labelsize': 11, 'ytick.major.pad': 3,
			'ytick.major.width': 1, 'ytick.major.size': 0, 'ytick.minor.width': 0.5,
			'ytick.minor.size': 0}

mpldefault = {} # Avoid loading all rcParams to minimize deprectation warnings
mpldefault.update((key,mpl.rcParams[key]) for key in ccgpaper.keys())

presentation = {'axes.axisbelow': True, 'axes.edgecolor': 'black', 'axes.facecolor': 'white',
				'axes.grid': False, 'axes.labelcolor': 'black', 'axes.labelsize': 8,
				'axes.linewidth': 0.5, 'axes.titlesize': 8, 'figure.figsize': (8, 8),
				'figure.subplot.bottom': 0.125, 'figure.facecolor': 'white',
				'figure.edgecolor': 'black',
				'font.family': 'Calibri', 'font.size': 12, 'font.weight': 400,
				'text.color': 'black', 'grid.color': 'lightgray', 'grid.linestyle': '-',
				'grid.linewidth': 0.5, 'legend.fontsize': 8.0, 'legend.frameon': False,
				'legend.numpoints': 1, 'legend.scatterpoints': 1, 'lines.linewidth': 1.5,
				'lines.markeredgewidth': 0, 'lines.markersize': 5, 'lines.solid_capstyle': 'round',
				'patch.linewidth': 0.3, 'pdf.fonttype': 3, 'ps.fonttype': 3, 'ps.useafm': True,
				'xtick.color': 'black', 'xtick.direction': 'out', 'xtick.labelsize': 8,
				'xtick.major.pad': 3, 'xtick.major.width': 1, 'xtick.major.size': 0,
				'xtick.minor.width': 0.5, 'xtick.minor.size': 0, 'ytick.color': 'black',
				'ytick.direction': 'out', 'ytick.labelsize': 8, 'ytick.major.pad': 3,
				'ytick.major.width': 1, 'ytick.major.size': 0, 'ytick.minor.width': 0.5,
				'ytick.minor.size': 0}

jdplots = {'axes.axisbelow': True, 'axes.edgecolor': 'k', 'axes.facecolor': 'white',
		   'axes.grid': False, 'axes.labelcolor': 'k', 'axes.labelsize': 'medium',
		   'axes.linewidth': 1.0, 'axes.titlesize': 'large',
		   'figure.figsize': (7, 5), 'figure.subplot.bottom': 0.125, 'figure.facecolor': 'white',
		   'figure.edgecolor': 'white',
		   'font.family': 'sans-serif', 'font.size': 12.0, 'font.weight': 'medium',
		   'grid.color': 'lightgray', 'grid.linestyle': '-', 'grid.linewidth': 0.5,
		   'legend.fontsize': 'medium', 'legend.frameon': False, 'legend.numpoints': 1,
		   'legend.scatterpoints': 1,
		   'lines.linewidth': 1.0, 'lines.markeredgewidth': 0.25, 'lines.markersize': 6,
		   'lines.solid_capstyle': 'round', 'patch.linewidth': 0.3,
		   'pdf.fonttype': 3, 'ps.fonttype': 3, 'ps.useafm': False,
		   'text.color': 'k', 'xtick.color': 'k', 'xtick.direction': 'in',
		   'xtick.labelsize': 'medium', 'xtick.major.pad': 4, 'xtick.major.width': 0.5,
		   'xtick.major.size': 4, 'xtick.minor.width': 0.5, 'xtick.minor.size': 0,
		   'ytick.color': 'k', 'ytick.direction': 'in', 'ytick.labelsize': 'medium',
		   'ytick.major.pad': 4, 'ytick.major.width': 0.5, 'ytick.major.size': 4,
		   'ytick.minor.width': 0.5, 'ytick.minor.size': 0}

# Create alternate styles of the default for different font sizes
fontparms = ['axes.labelsize', 'axes.titlesize', 'legend.fontsize', 'font.size', 'xtick.labelsize',
			 'ytick.labelsize']

# Define a dictionary pointing to the right data
plot_styles = {'nochanges': nochanges,
				'darkcontent': darkcontent,
				'ccgpaper': ccgpaper,
				'presentation': presentation,
				'mpldefault': mpldefault,
				'jdplots': jdplots}


class Parameters(MutableMapping, dict):

	"""
	A dictionary object for defining pygeostat default parameters, which is intended
	to add convenience to modeling workflows where common parameters must otherwise be
	repeated in function calls. An instance of this object is initialized as Parameters when
	pygeostat is imported. Default parameters are implemented as thoughtfully
	as possible, but may be modified by users once to suit their requirements. Keyword
	arguments for functions override any Parameters defaults for the function call itself,
	without modifying the defaults themselves.

	Examples of modifying plotting parameters include:

		>>> gs.Parameters['plotting.xname'] = 'Easting'
		>>> gs.Parameters['plotting.lagname'] = 'h'
		>>> gs.Parameters['plotting.unit'] = 'm'

	Examples of data parameters include:

		>>> gs.Parameters['data.dh'] = 'WellID'
		>>> gs.Parameters['data.griddef'] = pygeostat.GridDef()
		>>> gs.Parameters['data.x'] = 'Easting'

	The object performs validation on any modificaitons to the dictionary. If users are confused
	about the application of any defaults, a describe function allows for individual parameters
	(e.g., Parameters.describe('plotting.xname') or the entire dictionary (e.g., Parameters.describe()
	to be printed to the console.

	Attributation to Matplotlib:

		The template for the Parameters class is inspired by adapted from the
		Matplotlib rcParams class and related functionality, as it provides a template for the
		approach to setting defaults to a Python instance.
	"""
	# get validation items
	validate = dict((key, converter) for key, (_, converter) in six.iteritems(default_parameters))
	

	def __init__(self, *args, **kwargs):
		self.update(*args, **kwargs)
		self._section_name = 'PygeostatParameters'
		self._default_dict = dict((key, parameter) for key, (parameter, _) in six.iteritems(default_parameters))
		self._default_descriptions = default_parameter_descriptions


	def __setitem__(self, key, val):
		try:
			cval = self.validate[key](val)
			dict.__setitem__(self, key, cval)
		except ValueError as ve:
				raise ValueError("Key %s: %s" % (key, str(ve)))
		except:
			warnings.warn('%s is not a valid parameter! See Parameters.keys() for a '
						   'list of valid parameters.' % (key,))
	
	def __getitem__(self, key):
		val = dict.__getitem__(self, key)
		return val

	def __iter__(self):
		"""
		Yield sorted list of keys.
		"""
		for k in sorted(dict.__iter__(self)):
			yield k

	def __repr__(self):
		import pprint
		class_name = self.__class__.__name__
		indent = len(class_name) + 1
		repr_split = pprint.pformat(dict(self), indent=1,
									width=80 - indent).split('\n')
		repr_indented = ('\n' + ' ' * indent).join(repr_split)
		return '{0}({1})'.format(class_name, repr_indented)

	def __str__(self):
		return '\n'.join('{}: {}'.format(k, v) for k, v in sorted(self.items()))

	def __delitem__(self, key):
		raise KeyError('Not allowed to remove a paramater!')

	def find_all(self, pattern):
		"""
		Return the subset of this Parameters dictionary whose keys match,
		using :func:`re.search`, the given ``pattern``.

		.. note::

			Changes to the returned dictionary are *not* propagated to
			the parent Parameter dictionary.
		"""
		import re
		pattern_re = re.compile(pattern)
		return dict([(key, value) for key, value in self.items() if pattern_re.search(key)])

	def restore_defaults(self):
		"""
		Restore to the original dictionary of defaults values for the current instance
		"""
		self.update([(key, default_value) for key, default_value in self._default_dict.items()])


	def describe(self, key=None, verbose=True):
		"""
		Print a description of an individual parameter in the
		pygeostat parameter collection, as specified by its dictionary key. If key
		is None, all parameter descriptions are printed.

		Parameters:
			key(str) : dictionary key within Parameters

		"""
		if key is None:
			# Print the entire dictionary description
			for iterkey in iter(self):
				print(iterkey + ':\n' + self._default_descriptions[iterkey] + '\n')
		else:
			if key not in iter(self):
				raise KeyError('%s is not a valid parameter. See pygeostat.Parameters.keys() for a '
							   'list of valid parameters.' % (key,))
			if verbose:
				print(key + ':\n' + self._default_descriptions[key])

	def save(self, filename):
		"""
		Save the current pygeostat configurations to the given file
		"""
		_write_pardict(filename, self, self._section_name, description=self._default_descriptions)

	def load(self, filename):
		"""
		Load the pygeostat configurations from a given file
		"""
		import os
		assert os.path.isfile(filename), "ERROR: `{}` does not exist".format(filename)
		newdict = _parse_pardict(filename, self._section_name)
		self.update(newdict)

	def reset_systemdefault(self):
		"""
		Overwrite the system defaults stored in ``%USER%/.Pygeostat`` with the default parameters of
		pygeostat
		"""
		_write_pardict(_params_cfg_file, self._default_dict, self._section_name, description=self._default_descriptions)


	def set_systemdefault(self):
		"""
		Write the current configuration from the current instance i.e. ``self`` to the system defaults in ``%USER%/.Pygeostat``
		"""
		print("Writing current `{}` to system default at: {}".format(self._section_name, _params_cfg_file))
		main_dir = os.path.dirname(_params_cfg_file)
		if not os.path.exists(main_dir):
			os.makedirs(main_dir)
		self.save(_params_cfg_file)


	def get_systemdefault(self, confirm_autoload=False):
		"""
		Get the Parameters stored in the system defaults in ``%USER%/.Pygeostat``. ``confirm_autoload``
		is passed only on the very first call on module init.
		"""
		if os.path.isfile(_params_cfg_file):
			newdict = _parse_pardict(_params_cfg_file, section=self._section_name)
			
			if newdict: # If dictionary is not empty
				if confirm_autoload:
					if newdict.get("config.autoload.parameters", False):
						try:
							self.update(newdict)
							if newdict.get("config.verbose", False):
								print("Loading default Pygeostat Parameters from {}".format(_params_cfg_file))
						except:
							print('Unable to load the system defaults due to version issues.\nThe current system defaults setting was removed. Please try Pygeostat.Parameters.set_systemdefault()')
							os.remove(_params_cfg_file)
				else:
					self.update(newdict)
					if newdict.get("config.verbose", False):
						print("Loading default Pygeostat Parameters from {}".format(_params_cfg_file))


class PlotStyle(Parameters):
	"""
	Manages a set of pygeostat default plotting styles VIA ``mpl.rcParams``.

	Parameters:
		plot_style (str): One of: "darkcontent", "ccgpaper", "presentation", "mpldefault", "jdplots"
		custom_dict (dict): Dictionary of custom mpl.rcParams

	Note:
		Default behavior is to use the ``mpl`` defaults. ``PlotStyle.set_style()`` changes the defaults
		stored in this class for the current session

	Examples:

		The main instance of this class used by all plotting functions is accessible with
		``gs.PlotStyle``. The only things that the user should require from this instance is to
		either save the configured plot style to a given file:

		>>> gs.PlotStyle.save("plottingstyles.cfg")

		Load a saved configuration from a saved output:

		>>> gs.PlotStyle.load("plottingstyles.cfg")

		or to write system defaults, or load system defaults:

		>>> gs.PlotStyle.set_systemdefault()
		>>> gs.PlotStyle.get_systemdefault()

		one can optionally restore the defaults for the current session with:

		>>> gs.PlotStyle.restore_defaults()

		To set ``ccgpaper`` as the system default, the following is required:

		>>> gs.PlotStyle.set_style("ccgpaper")  # this is now the session default
		>>> gs.gsPlotStyle.set_systemdefault()

		``gsPlotStyle`` params are autoloaded if ``Parameters["config.autoload.gsplotstyle"] == True``

		>>> gs.Parameters["config.autoload.gsplotstyle"] = True
		>>> gs.Parameters.set_systemdefault()

	"""

	import matplotlib as mpl

	def __init__(self, plot_style, custom_dict=None):
		if plot_style in (None, False):
			plot_style = "nochanges"
		self.style_string = plot_style

		self._default_style = 'mpldefault'
		self._default_dict = plot_styles[self._default_style]

		self._default_descriptions = None
		self._section_name = 'PlotStyle'

		self.update(self._generate(plot_style, custom_dict))
	
	def set_style(self, plot_style, custom_dict=None):
		'''
		Sets the plot style based on the available styles
		'''
		dict.clear(self)
		if plot_style in (None, False):
			plot_style = "nochanges"
		self.update(self._generate(plot_style, custom_dict))

	def update(self, *args,**kwargs):
		super().update(*args,**kwargs)
		self.update_mplrcParams()

	def _check_style_parameter(self, key):
		
		if not key in self.mpl.rcParams.keys():
			raise KeyError
		# Ignore a list of deprecated matplotlib rcParams attributes
		if key in ['text.latex.unicode', 'axes.prop_cycle', 'examples.directory']:
			return False
		else:
			return True

	def __setitem__(self, key, val):
		'''
		Overriding the __setitem__ method to avoid setting parameter key words that are not required.
		'''
		try:
			if self._check_style_parameter(key):
				dict.__setitem__(self, key, val)
		except KeyError:
			raise KeyError('%s is not a valid matplotlib plotting parameter.' % (key,))

	def _generate(self, style=None, custom=None):
		"""
		The purpose of this method is to generate dictionaries using the ``style``, and
		allow ``custom`` keys to modify what is specified in the custom styles
		"""
		import copy

		if not isinstance(style, str):
			return {}
		if style in plot_styles.keys():
			newstyle = copy.copy(plot_styles[style])
			if isinstance(custom, dict):
				newstyle.update(custom)

		# allow 'pt' followed by a int or float number for the size of the font for all font related plot style prameters
		elif style.startswith('pt'):
			newstyle = copy.copy(ccgpaper)
			idx = style.rfind('t')
			newstyle.update(ccgpaper.fromkeys(fontparms, float(style[idx + 1:])))
			if isinstance(custom, dict):
				newstyle.update(custom)
		else:
			tempstyles = list(plot_styles.keys())
			tempstyles.append('pt##')
			raise ValueError("The style '%s' is not available in pygeostat. Please use one of the"
							 " following:\n%s" % (style, tempstyles))
		# make sure the keys are permissible mpl keys
		for key in newstyle:
			if key not in mpl.rcParams:
				raise KeyError("ERROR: {} not a valid rcParams key!".format(key))
		return newstyle
		

	def update_mplrcParams(self, style=None, custom=None):
		"""
		Update the current mpl rcParams by drawing from the set of pygeostat plotting styles,
		with optional custom parameters.

		Parameters:
			style (str): One of the pygeostat default plotting styles
			custom (dict): dictionary of mpl rcParams and values to set
		"""
		from matplotlib import cbook
		import warnings

		if Parameters['config.ignore_mpl_warnings']:
			warnings.filterwarnings("ignore", category=cbook.MatplotlibDeprecationWarning)
		else:
			warnings.filterwarnings('default', category=cbook.MatplotlibDeprecationWarning)

		# Required in order to restore mlp rcParams
		self._old_rcParams = self.mpl.rcParams.copy()
		if style in (None, False):
			if custom in (None, False):
				self.mpl.rcParams.update(self)
			else:
				tempkwargs = self.copy()
				tempkwargs.update(custom)
				self.mpl.rcParams.update(tempkwargs)
		else:
			tempkwargs = self._generate(style, custom)
			self.mpl.rcParams.update(tempkwargs)

	
	def restore_defaults(self):
		"""
		Restore the matplotlib defaults for this session
		"""
		self.update(self._default_dict)

	def restore_mplrcParams(self):
		""" Restores the `old` rcParams saved prior to the style update """
		mpl.rcParams.update(self._old_rcParams)


	def get_systemdefault(self, confirm_autoload=False, verbose = True):
		"""
		Get the PygeostatParams stored in the system defaults in ``%USER%/.PygeostatParams``. ``confirm_autoload``
		is passed only on the very first call on module init.
		"""
		if os.path.isfile(_params_cfg_file):
			newdict = _parse_pardict(_params_cfg_file, section=self._section_name)
			
			if newdict: # If dictionary is not empty
				if confirm_autoload:
					if Parameters.get("config.autoload.plot_style", False):
						self.update(newdict)
						if Parameters.get("config.verbose", False):
							print("Loading default Pygeostat Parameters from {}".format(_params_cfg_file))
				else:
					self.update(newdict)
					if Parameters.get("config.verbose", False):
						print("Loading default Pygeostat Parameters from {}".format(_params_cfg_file))

#--------------------Initilizations-------------------
#-----------------------------------------------------
# 1) Parameters comes first
def _parameters_initialize():
	ret = Parameters([(key, default) for key, (default, _) in
					six.iteritems(default_parameters)])
	ret.get_systemdefault(confirm_autoload=True)
	return ret

Parameters = _parameters_initialize()

# 2) Setting the plot styles 
def _plot_style_initilizer():
	ret = PlotStyle('ccgpaper')
	ret.get_systemdefault(confirm_autoload=True)
	return ret

PlotStyle = _plot_style_initilizer()

