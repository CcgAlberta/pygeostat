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
def accuracy_plot(truth=None, reals=None, mik_thresholds=None, acctype='sim', probability_increment=None,
		   figsize=None, title=None, xlabel=None, ylabel=None, stat_blk='standard',
		   stat_xy=(.95, .05), stat_fontsize=None, ms=40, color=None, label=None,
		   grid=None, axis_xy=None, ax=None, plot_style=None, custom_style=None, output_file=None,
		   reals_gaussian=None, data=None, weights=None, global_variance=1.0, **kwargs):
	"""
	Accuracy plot based on probability intervals quantified using an estimation technique e.g kriging.
	The distribution at each location is given either by realizations or by MIK probabilities.


	The statistics block shows Clayton Deutsch's accuracy A, precision P and goodness G (in
	percent), labelled A26, P26 and G26 to tell them from the original measures of 1996, with ``p``
	the probability intervals and ``t`` the fraction of the true values in each:

	* ``A = (1 - 2 * mean(max(0, p - t))) * 100``: 0 to 100. Points at or above the 1:1 line add
	  nothing; only points below it lower A
	* ``P = (1 - mean(local variance) / global variance) * 100``, with the variances in Gaussian
	  units. It is negative if the mean local variance is larger than the global variance (it is not
	  clipped). The true values are not used. Give the realizations in Gaussian units
	  (``reals_gaussian``), or the data that the normal score transform was made from (``data``)
	* ``G = 2 / (1/A + 1/P)``, the harmonic average of A and P. It lies between them

	P and G are shown as n/a without ``reals_gaussian`` or ``data``, and for MIK (not available
	yet). They replace the earlier accuracy, precision and goodness of pygeostat.

	The definitions are those of Deutsch, C. V. (2026), Refresh on the Accuracy Plot and Measures of
	Performance of Probabilistic Prediction, Paper 40, CCG Annual Report 28.

	Three statistics block sets are available: ``'minimal'``, the default ``'standard'`` and
	``'all'``. The statistics block can be customized to a user defined list and order. Available
	statistics are as follows:

	>>> ['ndat', 'nint', 'acc', 'pre', 'goo', 'avgvar', 'mse']

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
		probability_increment (float): Probability increment used during accuracy_plot calculation. If
			None, 0.1 (9 intervals, 10% to 90%) for fewer than 300 true values and 0.05 (19 intervals,
			5% to 95%) otherwise
		figsize (tuple): Figure size (width, height)
		title (str): Title for the plot
		xlabel (str): X-axis label (default "Predicted Probability")
		ylabel (str): Y-axis label (default "True Fraction")
		stat_blk (str, list or bool): A statistics block set name, a list of statistics, or
			``False`` to plot no statistics
		stat_xy (float tuple): X, Y coordinates of the annotated statistics in figure space. The
			coordinates specify the top right corner of the text block
		stat_fontsize (float): the fontsize for the statistics block. If None, based on
			gs.Parameters['plotting.stat_fontsize']. If less than 1, it is the fraction of the
			matplotlib.rcParams['font.size']. If greater than 1, it the absolute font size.
		ms (float): Size of scatter plot markers
		color: Colour of the markers. If None (the default), a marker and its line to the 1:1 line
			are green where the fraction in the interval is at or above the 1:1 line and red where
			it is below. With a colour, all markers have it, and the lines are the same colour
			(red for black markers)
		label (str): Legend label, also the first line of the statistics block. Useful when
			several models are compared on the same ``ax``
		grid(bool): plots light grey major grid lines behind the plot if True. On if None.
		axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
			visible) if axis_xy is True. Based on gs.Parameters['plotting.axis_xy'] if None.
		ax (mpl.axis): Matplotlib axis to plot the figure
		plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
			to turn it off
		custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
		output_file (str): Output figure file name and location
		reals_gaussian: The realizations in Gaussian units (before the back transform), the same shape
			as ``reals``, for P and G
		data: The data values that the normal score transform of the simulation was made from. If
			``reals_gaussian`` is None, the same transform is applied to ``reals`` for P and G
		weights: Declustering weights of ``data``, if the transform used them
		global_variance (float): Variance of the transformed data, 1 for a normal score transform
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
	import numpy as np
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
		propavg, statlist = accsim(truth, reals, probability_increment, reals_gaussian=reals_gaussian,
								   data=data, weights=weights, global_variance=global_variance)
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
	if color is None:
		# Green at or above the 1:1 line, red below it
		point_colors = line_colors = np.where(y >= x, 'tab:green', 'tab:red')
		text_color = 'k'
	else:
		point_colors, text_color = color, color
		line_colors = 'red' if same_color(color, 'k') else color
	ax.vlines(x, x, y, color=line_colors, lw=1, zorder=1)
	ax.scatter(x, y, facecolors=point_colors, edgecolors='k', linewidths=0.5, s=ms, label=label,
			   zorder=2)
	ax.plot((1, 0), (1, 0), zorder=1, lw=0.5, color='k')
	# Configure plot
	if xlabel is None:
		xlabel = 'Predicted Probability'
	if ylabel is None:
		ylabel = 'True Fraction'
	if grid is None:
		grid = True
	format_plot(ax, xlabel=xlabel, ylabel=ylabel, title=title,
			  grid=grid, axis_xy=axis_xy, xlim=(0, 1.0), ylim=(0, 1.0))
	if grid:
		ax.grid(True, color='0.7', lw=0.6)
		ax.set_axisbelow(True)
	# Plot stats
	if stat_blk:
		statlist = {'ndat': '$n = {}$'.format(ndat),
					'nint': '$n_p = {}$'.format(len(x)),
					'avgvar': '$U = {:0.3g}$'.format(statlist['avgvar']),
					'mse': '$MSE = {:0.3g}$'.format(statlist['mse']),
					'acc': r'$A_{{26}} = {:0.1f}\%$'.format(statlist['acc']),
					'pre': r'$P_{{26}} = {:0.1f}\%$'.format(statlist['pre']) if np.isfinite(statlist['pre'])
					else '$P_{26} = n/a$',
					'goo': r'$G_{{26}} = {:0.1f}\%$'.format(statlist['goo']) if np.isfinite(statlist['goo'])
					else '$G_{26} = n/a$'}
		statsets = {'minimal': ['nint'],
					'standard': ['ndat', 'nint', 'acc', 'pre', 'goo', 'avgvar', 'mse'],
					"all": ['ndat', 'nint', 'acc', 'pre', 'goo', 'avgvar', 'mse']}
		# Get the stat block text
		txtstats, stat_xy, ha, va = get_statblk(stat_blk, statsets, statlist, stat_xy)
		if label is not None:
			txtstats = label + '\n' + txtstats
		# Write the stat block to the plot
		stat_fontsize = _set_stat_fontsize(stat_fontsize)
		ax.text(stat_xy[0], stat_xy[1], txtstats, va=va, ha=ha, transform=ax.transAxes,
				fontsize=stat_fontsize, color=text_color)
	# Export figure
	if output_file or ('pdfpages' in kwargs):
		export_image(output_file, **kwargs)

	return ax
