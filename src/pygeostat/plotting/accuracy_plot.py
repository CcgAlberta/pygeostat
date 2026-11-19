# !/usr/bin/env python
#  -*- coding: utf-8 -*-
"""accuracy_plot plots the fraction of the true values that fall within a specified probability interval"""

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from . set_style import set_plot_style

@set_plot_style
def accuracy_plot(truth=None, reals=None, mik_thresholds=None, acctype='sim', probability_increment=0.05,
		   figsize=None, title=None, xlabel=None, ylabel=None, stat_blk='standard',
		   stat_xy=(.95, .05), stat_fontsize=None, ms=5, color='k', label=None,
		   grid=None, axis_xy=None, ax=None, plot_style=None, custom_style=None, output_file=None,
		   **kwargs):
	"""
	Accuracy plot based on probability intervals quantified using an estimation technique e.g kriging.
	The distribution at each location is given either by realizations or by MIK probabilities.


	Three statistics block sets are available: ``'minimal'``, the default ``'standard'`` and
	``'all'``. The statistics block can be customized to a user defined list and order. Available
	statistics are as follows:

	>>> ['ndat', 'nint', 'avgvar', 'mse', 'acc', 'pre', 'goo']

	Please review the documentation of the :func:`gs.PlotStyle.set_style()
	<pygeostat.pygeostat_parameters.PlotStyle.set_style>` and :func:`gs.export_image()
	<pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
	their use in this function can be understood.

	Keyword Arguments:
	   
		truth: Tidy (long-form) 1D data where a single column containing the true values.
			A pandas dataframe/series or numpy array can be passed
		reals: Tidy (long-form) 2D data where a single column contains values from a single
			realizations and each row contains the simulated values from a single truth location. “reals” has the same 
			number of data points as “truth”, that is , each row of “reals” corresponds to the true value in the same row of “truth”.
			A pandas dataframe or numpy matrix can be passed
		mik_thresholds (np.ndarray): 1D array of the z-vals ``mik_thresholds`` corresponding to the
			probabilities defined in reals for each location
		acctype (str): Currently ``sim`` and ``mik`` are valid. if ``mik_thresholds`` is passed the
			type is assumed to be ``mik``
		probability_increment (float): Probability increment used during accuracy_plot calculation
		figsize (tuple): Figure size (width, height)
		title (str): Title for the plot
		xlabel (str): X-axis label
		ylabel (str): Y-axis label
		stat_blk (str, list or bool): A statistics block set name, a list of statistics, or
			``False`` to plot no statistics
		stat_xy (float tuple): X, Y coordinates of the annotated statistics in figure space. The
			coordinates specify the top right corner of the text block
		stat_fontsize (float): the fontsize for the statistics block. If None, based on
			gs.Parameters['plotting.stat_fontsize']. If less than 1, it is the fraction of the
			matplotlib.rcParams['font.size']. If greater than 1, it the absolute font size.
		ms (float): Size of scatter plot markers
		color: Colour of the markers. The lines from each marker to the 1:1 line are red for
			black markers, otherwise the same colour as the markers
		label (str): Legend label, also the first line of the statistics block. Useful when
			several models are compared on the same ``ax``
		grid(bool): plots the major grid lines if True. Based on gs.Parameters['plotting.grid']
			if None.
		axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
			visible) if axis_xy is True. Based on gs.Parameters['plotting.axis_xy'] if None.
		ax (mpl.axis): Matplotlib axis to plot the figure
		plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
			to turn it off
		custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
		output_file (str): Output figure file name and location
		**kwargs: Optional permissible keyword arguments to pass to
			:func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`

	Returns:
		ax (ax): Matplotlib Axes object with the accuracy plot


	Note: 
		The simulated values in each row of “reals” are used to calculate an empirical CDF using the standard method of 
		the midpoint of a histogram bin and the CDF is constructed using all points (as done in GSLIB). When "mik_thresholds" is None, 
		“reals” contains the simulated values from which a cdf should be computed. Otherwise the “reals” contains the distribution values, 
		F(mik_thresholds), for each location (each row of the data file).


	Examples:


	A simple call using truth and realization data: in this example, the first column of data_file is "truth" 
	and the rest of columns are "reals".

	.. plot::

		import pygeostat as gs
		data_file = gs.ExampleData('accuracy_plot')
		reals = data_file[list(data_file.columns[1:])].values
		truth = data_file[list(data_file.columns)[0]].values
		gs.accuracy_plot(truth=truth, reals=reals)

	|

	Compare two models on one plot. Here the second model uses only the first 10 realizations.

	.. plot::

		import pygeostat as gs
		data_file = gs.ExampleData('accuracy_plot')
		reals = data_file[list(data_file.columns[1:])].values
		truth = data_file[list(data_file.columns)[0]].values
		ax = gs.accuracy_plot(truth=truth, reals=reals, label='100 realizations',
		                      stat_blk='minimal')
		gs.accuracy_plot(truth=truth, reals=reals[:, :10], label='10 realizations', color='C1',
		                 stat_blk='minimal', stat_xy=(0.95, 0.2), ax=ax)
		ax.legend(loc='upper left')

	"""
	from matplotlib.colors import same_color
	from ..statistics.utils import accsim, accmik
	from .utils import setup_plot, get_statblk, format_plot, _set_stat_fontsize
	from .export_image import export_image

	# Sanity checks based on the provided input arguments
	simcheck = [(val is not None) for val in [truth, reals]]
	mikcheck = [(val is not None) for val in [mik_thresholds]]
	if any(simcheck) and not all(simcheck):
		raise ValueError("The arguments `truth` and `reals` are required")
	if any(mikcheck) and not all(simcheck):
		raise ValueError("Must have `truth`, `reals` and `mik_thresholds` for MIK !")
	if all(simcheck) and acctype.lower() != 'sim':
		acctype = 'sim'
	if all(mikcheck) and acctype == 'sim':
		acctype = 'mik'
	# Handle accuracy_plot-sim functionality
	if acctype.lower() == 'sim':
		propavg, statlist = accsim(truth, reals, probability_increment)
	elif acctype.lower() == 'mik':
		propavg, statlist = accmik(truth, mik_thresholds, reals, probability_increment)
	else:
		raise ValueError ('Wrong acctype was provided.')
	x = propavg['ProbInt'].values
	y = propavg['FracIn'].values
	ndat = len(truth)

	fig, ax, cax = setup_plot(ax, cbar=False, figsize=figsize)
	ax.set_aspect('equal')
	# Plot the figure
	line_color = 'red' if same_color(color, 'k') else color
	ax.vlines(x, x, y, color=line_color, lw=1, zorder=1)
	ax.scatter(x, y, color=color, s=ms, label=label, zorder=2)
	ax.plot((1, 0), (1, 0), zorder=1, lw=0.5, color='k')
	# Configure plot
	if xlabel is None:
		xlabel = 'Probability Interval'
	if ylabel is None:
		ylabel = 'Fraction in Interval'
	format_plot(ax, xlabel=xlabel, ylabel=ylabel, title=title,
			  grid=grid, axis_xy=axis_xy, xlim=(0, 1.0), ylim=(0, 1.0))
	# Plot stats
	if stat_blk:
		statlist = {'ndat': '$n = {}$'.format(ndat),
					'nint': '$n_p = {}$'.format(len(x)),
					'avgvar': '$U = {:0.3g}$'.format(statlist['avgvar']),
					'mse': '$MSE = {:0.3g}$'.format(statlist['mse']),
					'acc': '$Accuracy = {:0.3f}$'.format(statlist['acc']),
					'pre': '$Precision = {:0.3f}$'.format(statlist['pre']),
					'goo': '$Goodness = {:0.3f}$'.format(statlist['goo'])}
		statsets = {'minimal': ['nint'],
					'standard': ['ndat', 'nint', 'avgvar', 'mse'],
					"all": ['ndat', 'nint', 'avgvar', 'mse', 'acc', 'pre', 'goo']}
		# Get the stat block text
		txtstats, stat_xy, ha, va = get_statblk(stat_blk, statsets, statlist, stat_xy)
		if label is not None:
			txtstats = label + '\n' + txtstats
		# Write the stat block to the plot
		stat_fontsize = _set_stat_fontsize(stat_fontsize)
		ax.text(stat_xy[0], stat_xy[1], txtstats, va=va, ha=ha, transform=ax.transAxes,
				fontsize=stat_fontsize, color=color)
	# Export figure
	if output_file or ('pdfpages' in kwargs):
		export_image(output_file, **kwargs)

	return ax
