#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot multi dimensional scaling (MDS) coordinates.
"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from . set_style import set_plot_style

@set_plot_style
def mds_plot(data, variables=None, enumerate_vars=True, ax=None, figsize=(3, 3), title=None, aspect = None,
			cmap='viridis', s=15, varlabels=None, plot_style=None, custom_style=None, output_file=None, out_kws=None,
			**kwargs):
	"""
	Python implementation of the MDS coordinates ploted by the CCG program corrmat when set to
	'ordination' mode.

	MDS coordinates are calculate from the correlation matrix.

	.. seealso::

		1. Deutsch, M. V, & Deutsch, C. V. (2013). A Program to Calculate and Display Correlation
		   Matrices with Relevant Analysis. Edmonton AB. Retrieved from http://www.ccgalberta.com

	Parameters:
		data: Tidy (long-form) 2-D data where each column is a variable and each row is an
			observation. A pandas dataframe or numpy array may be passed.
		variables (list): Variables from the pd.DataFrame passed with ``data`` to calculate
			coordinates for
		enumerate_vars (bool): Indicate if the scatter points should be annotated with an
			enumerator rather then the variable name. A reference legend is placed to the right
			of the plot if set to ``True``
		ax (mpl.axis): Matplotlib axis to plot the figure on
		figsize (tuple): Figure size (width, height) in inches
		title (str): Title for the plot
		cmap (str): A matplotlib colormap object or a registered matplotlib or pygeostat colormap
			name.
		s (float): Size of location map markers
		plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
			to turn it off
		custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
		output_file (str): Output figure file name and location
		out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
			:func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
		**kwargs: Optional dictionary of permissible keyword arguments to pass to the scatter plot

	Examples:
		A simple call:

		>>> gs.mds_plot(data.data)

		.. image:: ./figures/mds_plot1.png

	"""

	from .cmaps import avail_cmaps
	from .utils import get_contcbarargs, setup_plot, smart_annotate, get_cmap
	from .export_image import export_image
	from .. multivariate.utils import mds
	import numpy as np

   
	# Default dictionary
	if not out_kws:
		out_kws = dict()
	# Function set-up and defaults
	if variables is None:
		try:
			variables = list(data.columns)
		except:
			variables = list(np.arange(data.data.shape[1]))
	if enumerate_vars is True:
		annots = list(range(1, len(variables) + 1))
	else:
		annots = variables
	if cmap in avail_cmaps:
		cmap = get_cmap(cmap)
	# Calculate the MDS coordinates
	coords = mds(data, variables)
	xcoord = coords['X1']
	ycoord = coords['X2']
	zcoord = coords['X3']
	# Make the figure object
	fig, ax, cax = setup_plot(ax, cax=None, cbar=True, figsize=figsize, **kwargs)
	ax.set_aspect(aspect)
	# Plot the scatter plot using the MDS coordinates
	vlim, ticklocs, ticklabels = get_contcbarargs(zcoord, 2, None)
	plot = ax.scatter(xcoord, ycoord, c=zcoord, cmap=cmap, vmin=vlim[0], vmax=vlim[1], s=s,
					  alpha=0.8)
	# Plot scatter plot outlines to better show overlapping
	ax.scatter(xcoord, ycoord, c='none', s=s)
	# Adjust the annotations of each point
	smart_annotate(ax, xcoord, ycoord, annots)
	# Set-up the rest of the plots
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	if title is not None:
		ax.set_title(title)
	ax.tick_params(axis='both', pad=2)
	ax.set_xlabel('X1')
	ax.set_ylabel('X2')
	# Configure the color bar
	cbar = fig.colorbar(plot, cax=cax, ticks=ticklocs)
	cbar.ax.set_yticklabels(ticklabels, ha='left')
	cbar.ax.tick_params(axis='y', pad=2)
	cbar.set_label('X3', ha='center', va='top', labelpad=2)
	# Write the legend if the variables are enumerated
	if enumerate_vars:
		renderer = fig.canvas.get_renderer()
		fig.draw(renderer)
		pad = len('%s' % len(variables))
		colheight = cbar.ax.yaxis.get_tightbbox(renderer).height
		txt = cbar.ax.yaxis.label
		txtx = ax.transAxes.inverted().transform((txt.get_window_extent(renderer).xmax,
												  txt.get_window_extent(renderer).xmax))[0] + 0.02
		legend = ''
		nlines = 0
		if varlabels is None:
			varlabels = variables
		for i, var in enumerate(varlabels):
			legend = legend + '{number:{pad}d}: {var}   \n'.format(pad=pad, number=(i + 1),
																   var=var)
			if nlines == 0:
				txt = ax.text(txtx, 1, legend + '\nTest', ha='left', va='top',
							  transform=ax.transAxes)
				height = txt.get_window_extent(renderer).height
				txt.remove()
				if height > colheight:
					nlines = (i + 1)
			if (nlines > 0) and (((i % nlines) == (nlines - 1)) or
								 (i == len(variables) - 1)):
				txt = ax.text(txtx, 1, legend, ha='left', va='top', transform=ax.transAxes)
				txtx = ax.transAxes.inverted().transform((txt.get_window_extent(renderer).xmax,
													  txt.get_window_extent(renderer).xmax))[0]
				legend = ''
			elif var == varlabels[-1]:
				txt = ax.text(txtx, 1, legend, ha='left', va='top', transform=ax.transAxes)
	# Export figure
	if output_file or ('pdfpages' in out_kws):
		export_image(output_file, **out_kws)

	return ax
