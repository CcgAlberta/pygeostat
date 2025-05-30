#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""plotting a 2D gridded dataset or a slice of a 3D gridded dataset with the ability to include another layer to show point data"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from . set_style import set_plot_style
import numpy as np
import matplotlib as mpl

@set_plot_style
def slice_plot(data, griddef=None, var=None, catdata=None, pointdata=None, pointvar=None,
			 pointtol=None, pointkws=None, pointcmap=None, orient='xy', slice_number=0, ax=None,
			 figsize=None, vlim=None, clim=None, title=None, xlabel=None, ylabel=None, unit=None,
			 rotateticks=None, cbar=True, cbar_label=None, catdict=None, cbar_label_pad=None,
			 cax=None, sigfigs=3, cmap=None, interp='none', aspect=None, grid=None, axis_xy=None,
			 rasterize=False, plot_style=None, custom_style=None, output_file=None, out_kws=None,
			 return_cbar=False, return_plot=False, slice_thickness=None, cbar_nticks=5,
			 plotformat_dict=None, **kwargs):
	"""
	slice_plot displays a 2D gridded dataset or a slice of a 3D gridded dataset. To plot only
	scattered data, please see :func:`gs.location_plot()<pygeostat.plotting.location_plot>`

	The only required parameters are ``data`` and ``griddef``. All other parameters are optional or
	calculated automatically. Axis tick labels are automatically checked for overlap and if needed,
	are rotated. The figure instance is always returned. To allow for easier modification, the
	colorbar object and the data used to generate the plot can also be returned. Examples of their
	use are provided bellow.

	The values use to bound the data (i.e., vmin and vmax) are automatically calculated by default.
	These values are determined based on the number of significant figures and the sliced data;
	depending on data and the precision specified, scientific notation may be used for the colorbar
	tick lables. When point data shares the same colormap as the gridded data, the points displayed
	are integrated into the above calculation.

	Categorical data can be used, however ``catdata`` will need to be set to ``True`` for proper
	implementation. Categorical colour palettes are available within pygeostat. See the
	documentation for :func:`gs.get_palette() <pygeostat.plotting.utils.get_palette>` for more
	information.

	Please review the documentation of the :func:`gs.set_style()
	<pygeostat.PlotStyle.set_style>` and :func:`gs.export_image()
	<pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
	their use in this function can be understood.

	Parameters:
		data: A numpy ndarray, pandas DataFrame or pygeostat DataFile, where each column is a
			variable and each row is an observation
		griddef (GridDef): A pygeostat GridDef class, which must be provided if a
			DataFile is not passed as data with a valid internal GridDef
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`
		var (str,int): The name of the column within data to plot. If an int is provided, then it
			corresponds with the column number in data. If None, the first column of data
			is used.
		catdata (bool): Indicate if the data is categorical or not. Will be automatically set if
			less than Parameters['plotting.assumecat'] unique values are found within ``data``
		pointdata (DataFile): A pygeostat DataFile class created using
			:class:`gs.DataFile <pygeostat.data.data.DataFile>` containg point data to overlay the
			gridded data with. Must have the necessary coordinate column headers stored within the
			class.
		pointvar (str): Column header of variable to plot within ``pointdata`` or a permissible
			matplotlib colour
		pointtol (float): Slice tolerance to plot point data (i.e. plot +/- ``pointtol`` from the
			center of the slice). Any negative value plots all data. Default is to plot all data.
		pointkws (dict): Optional dictionary of permissible keyword arguments to pass to
			matplotlib's scatter function. Default values is ``{'marker':'o', 's':15}``
		orient (str): Orientation to slice data. ``'xy'``, ``'xz'``, ``'yz'`` are the only accepted
			values
		slice_number (int): Grid cell location along the axis not plotted to take the slice of data to
			plot
		ax (mpl.axis): Matplotlib axis to plot the figure
		figsize (tuple): Figure size (width, height)
		vlim (float tuple): Data minimum and maximum values
		clim (int tuple) or (list): Categorical data minimum and maximum values, Forces categorical
			colorbar to plot the full range of categorical values - even if none show in the plot.
			Can be either a tuple of integer values OR a list of integers.
		title (str): Title for the plot. If left to it's default value of ``None`` or is set to
			``True``, a logical default title will be generated for 3-D data. Set to ``False`` if
			no title is desired.
		xlabel (str): X-axis label
		ylabel (str): Y-axis label
		unit (str): Unit to place inside the axis label parentheses
		rotateticks (bool tuple): Indicate if the axis tick labels should be rotated (x, y)
		cbar (bool): Indicate if a colorbar should be plotted or not
		cbar_label (str): Colorbar title
		catdict (dict): Dictionary to map enumerated catagories to names (e.g., {100: 'Unit1'}).
			Taken from Parameters['data.catdict'] if catdata=True and its keys align.
		sigfigs (int): Number of sigfigs to consider for the colorbar
		cmap (str): Matplotlib or pygeostat colormap or palette
		interp (str): plt.imshow interpolation option; ``'spline36'`` (continuous) and
			``'hermite'`` (categorical) are good starting points if smoothing is desired.
			``'none'`` is the default setting
		aspect (str): Set a permissible aspect ratio of the image to pass to matplotlib. If None,
			it will be 'equal' if each axis is within 1/5 of the length of the other. Otherwise,
			it will be 'auto'.
		grid (bool): Plots the major grid lines if True. Based on Parameters['plotting.grid']
			if None.
		axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
			visible) if axis_xy is True. Based on Parameters['plotting.axis_xy'] if None.
		rasterize (bool): Indicate if the gridded image should be rasterized during export. The
			output resolution is depended on the DPI setting used during export.
		plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
			to turn it off
		custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
		output_file (str): Output figure file name and location
		out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
			:func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
		return_cbar (bool): Indicate if the colorbar axis should be returned. A tuple is returned
			with the first item being the axis object and the second the cbar object.
		return_plot (bool): Indicate if the plot from imshow should be returned. It can be used
			to create the colorbars required for subplotting with the ImageGrid()
		slice_thickness (int): Number of slices around the selected slice_number to average the attributes
		**kwargs: Optional permissible keyword arguments to pass to matplotlib's imshow function

	Returns:
		Default:
		ax (ax): Matplotlib axis instance which contains the gridded figure
		return_cbar:
		ax (ax), cbar (cbar): Optional, default False. Matplotlib colorbar object
		retrun_plot:
		ax (ax), plot(?): Optional, default False. Matplotlib colorbar object
		return_cbar & return_plot
		ax, plot, cbar: default False

	**Examples:**

		A simple call by passing the data file

		.. plot::

			import pygeostat as gs
			griddef = gs.GridDef([40,0.5,1,40,0.5,1,40,0.5,1])
			data_file = gs.ExampleData('3d_grid', griddef)
			gs.slice_plot(data_file, orient='xy', cmap='viridis')

		|

		Plot with slice thickness averaging/upscaling:

		.. plot::

			import pygeostat as gs
			griddef = gs.GridDef([40,0.5,1,40,0.5,1,40,0.5,1])
			data_file = gs.ExampleData('3d_grid', griddef)
			gs.slice_plot(data_file, orient='xy', cmap='viridis', slice_thickness=20)


	.. codeauthor:: pygeostat development team 2016-04-11
	"""

	import matplotlib.pyplot as plt
	from . utils import (_spatial_labels, _spatial_pointdata, _spatial_griddata, _spatial_slice,
						 _format_tick_labels, format_plot, _spatial_aspect, setup_plot, format_subplot_axis)

	from .. pygeostat_parameters import Parameters
	from .export_image import export_image
	from ..datautils.utils import slice_grid
	from ..data.data import DataFile


	# Handle dictionary defaults
	if out_kws is None:
		out_kws = dict()
	# Handle var
	if var is None:
		if isinstance(data, DataFile):
			if data.variables is not None:
				if isinstance(data.variables, str):
					var = data.variables
	# Handle pointvar
	if pointdata is not None:
		if pointvar is None:
			if isinstance(var, str):
				if var in list(pointdata.columns):
					pointvar = var
			if pointvar is None:
				pointvar = Parameters['plotting.location_plot.c']
	# Handle title
	if title is None:
		if hasattr(data, 'columns'):
			if var in data.columns:
				title = var
	# Parse the data, var and griddef input to determine the data and griddef
	data, griddef = _spatial_griddata(data, var, griddef)
	# some threshold of uniqueness means its likely categorical?
	# also the number of max categories in the palettes
	nunique = len(np.unique(data)[np.isfinite(np.unique(data))])
	if catdata is None and nunique < Parameters['plotting.assumecat']:
		catdata = True
	elif catdata is None and isinstance(catdict, dict):
		catdata = True
	elif catdata is None:  # back to original default argument
		catdata = False
	# Slice the data
	if orient in ['xy', 'xz', 'yz']:
		view = slice_grid(data, griddef, orient, slice_number, slice_thickness=slice_thickness)
	else:
		raise ValueError("Not a valid orientation! {}".format(orient))
	# Validate then slice the point data if required
	try:
		if pointtol < 0:
			pointtol = None
	except:
		pass
	if pointdata is not None:
		# Convert the pointdata into the data frame, x, y and z
		pointdata, pointx, pointy, pointz, _ = _spatial_pointdata(pointdata, orient)
		# Slice the data
		pointx, pointy, pointvari = _spatial_slice(pointdata, pointvar, pointx, pointy,
												   pointz, griddef, orient, slice_number, pointtol)
		if pointx is not None:
			if pointvar in pointdata.columns:
				if pointcmap is None:
					pointcmap = True
		else:
			pointcmap = None
	else:
		pointvari = None
		pointcmap = None
	# Handel grid extent
	a = orient[0]
	b = orient[1]
	xmin = getattr(griddef, a + 'limits')[0]
	xmax = getattr(griddef, a + 'limits')[1]
	ymin = getattr(griddef, b + 'limits')[0]
	ymax = getattr(griddef, b + 'limits')[1]
	if aspect is None:
		aspect = _spatial_aspect([xmin, xmax], [ymin, ymax])
	# handle parms for colorbars and colormaps
	output = _color_handling_gridded(data, cmap=cmap, catdata=catdata,
									 pointcmap=pointcmap, pointvar=pointvari,
									 view=view, sigfigs=sigfigs, vlim=vlim, clim=clim,
									 nticks=cbar_nticks, catdict=catdict)
	ticklabels, vlim, ticklocs, ncat, pointvari, cmap = output
	# Make sure the chosen colormap is set to the point data if required
	if pointcmap is True or pointcmap is None:
		pointcmap = cmap
	# Handel categorical data
	if catdata:
		tview = view.copy()
		for i in range(ncat):
			tview[view == ticklabels[i]] = i
		if catdict is None:
			catdict = Parameters['data.catdict']
		if isinstance(catdict, dict):
			ticklabels = [catdict[cat] for cat in ticklabels]
	else:
		tview = view.copy()
	# Set-up plot if no axis is supplied using the ImageGrid method if required or the regular way
	fig, ax, cax = setup_plot(ax, cbar=cbar, cax=cax, figsize=figsize)
	# Handle interpolation in aspect has changed
	if aspect != 'equal' and interp == 'none':
		if catdata:
			interp = 'hermite'
		else:
			interp = 'spline36'
	# Plot gridded data
	plot = ax.imshow(tview, origin='lower', extent=[xmin, xmax, ymin, ymax], aspect=aspect,
					 interpolation=interp, vmin=vlim[0], vmax=vlim[1], cmap=cmap, **kwargs)
	# Rasterize the gridded plot
	if rasterize:
		plot.set_rasterized(True)
	# Plot scattered data if required
	if pointvari is not None:
		# Handel default dictionary
		if pointkws is None:
			pointkws = {'marker': 'o', 's': Parameters['plotting.location_plot.s']}
		if 'marker' not in pointkws:
			pointkws['marker'] = 'o'
		if 's' not in pointkws:
			pointkws['s'] = Parameters['plotting.location_plot.s']
		# if 'edgecolors' not in pointkws:
		#     pointkws['edgecolors'] = 'k'
		ax.scatter(pointx, pointy, c=pointvari, cmap=pointcmap,
				   vmin=vlim[0], vmax=vlim[1], **pointkws)
		ax.set_xlim((xmin, xmax))
		ax.set_ylim((ymin, ymax))
	# Plot labels
	ax = _spatial_labels(ax, orient, griddef, slice_number=slice_number, title=title, xlabel=xlabel,
						 ylabel=ylabel, unit=unit)
	# Note on Tick Labels:
	#   If a group of subplots are put together which share a x-axis, the rotation may not work. By
	#   getting the tick labels generated by matplotlib as a set of label objects, they can be
	#   looped through and have their settings individually fixed. This appears to be the only way
	#   to have the shared axis labels formated properly. The labels are also adjusted closer to
	#   the axis for esthetics. --Warren E. Black

	# The plots tick labels will not be properly accessible until the figure is "drawn", once the
	# command below is run, ax.get_ticklabel() will actually work properly.
	plt.draw()
	# Check to see if the ytick labels need to be rotated if the rotate argument was not passed.
	# This doesn't work if the figure is being placed into a subplot (i.e., is not a standalone
	# figure). Don't know why...
	_format_tick_labels(ax, rotateticks=rotateticks)
	# Note on Colorbar:
	#   It is my opinion that colorbars are more esthetically pleasing if they extend the full
	#   height of the figure it represents. The method that shifts the plot to the left distorts
	#   its aspect. A divider is used to accomplish the desired effect without distorting the
	#   figure. It is also important that the tick labels indicate that - in the case of a
	#   continuous variable - the min and max include their respective greater/less than or equal
	#    to symbols if they are not the absolute max/min values. --Warren E. Black
	# Plot colorbar, the tick locations and labels were generated above using gs.get_vlimdata
	if cbar:
		# Plot the colorbar
		cbar = fig.colorbar(plot, cax=cax, ticks=ticklocs)
		# Configure the color bar
		if cbar_label_pad is None:
			cbar_label_pad = 2
			for lbl in ticklabels:
				if isinstance(lbl, str):
					if 'geq' in lbl or 'leq' in lbl:
						cbar_label_pad = -6 + sigfigs
		cbar.ax.set_yticklabels(ticklabels, ha='left')
		cbar.ax.tick_params(axis='y', pad=2)
		if cbar_label is not None:
			cbar.set_label(cbar_label, ha='center', va='top', labelpad=cbar_label_pad)

	# Format plot based on plotformat_dict if it's available
	if plotformat_dict and cbar:
		fig, ax, cbar = format_subplot_axis(fig, ax, cbar=cbar,
											   plotformat_dict=plotformat_dict)
	elif plotformat_dict:
		fig, ax = format_subplot_axis(fig, ax, cbar=None,
										 plotformat_dict=plotformat_dict)
	if axis_xy is None:
		axis_xy = Parameters['plotting.axis_xy_spatial']
	format_plot(ax, grid=grid, axis_xy=axis_xy)

	# Export figure
	if output_file or ('pdfpages' in out_kws and out_kws['pdfpages']):
		export_image(output_file, **out_kws)
	# Return whats needed
	if return_cbar:
		if return_plot:
			return ax, plot, cbar
		else:
			return ax, cbar
	elif return_plot:
		return ax, plot
	else:
		return ax


def _color_handling_gridded(data, cmap=None, catdata=None, pointcmap=None,
							pointvar=None, view=None, sigfigs=3, vlim=None,
							clim=None, ticklabels=None, ticklocs=None, nticks=5,
							catdict=None):
	'''A small utility function called from many of the plotting functions. This sets the colorbar
	and colormap parameters and loads the specified colormap.

	Parameters:
		data: Tidy (long-form) dataframe where each column is a variable and each row is an
			observation. Pandas dataframe or numpy array
		cmap (str): Matplotlib or pygeostat colormap or palette
		catdata (bool): Indicate if the data is categorical or not. Will be automatically set if
			less than Parameters['plotting.assumecat'] unique values are found within ``data``
		pointcmap (bool): True to color the points using cmap
		pointvar: a pandas series of the point data used for the colormap. If categorical data is
			passed this will get converted to sequential numbers used in the colormap object.
		view: A slice of the data using ``gs.slice_grid()``
		sigfigs (int): Number of sigfigs to consider for the colorbar
		vlim (float tuple): Data minimum and maximum values
		clim (int tuple) or (list): Categorical data minimum and maximum values, Forces categorical
			colorbar to plot the full range of categorical values - even if none show in the plot.
			Can be either a tuple of integer values OR a list of integers.
		ticklabels: Filler. Otherwise you seem to have to set this variable before calling the
			function otherwise it gives you an error about assigning a value before the variable is
			defined
		ticklocs: Filler. Otherwise you seem to have to set this variable before calling the
			function otherwise it gives you an error about assigning a value before the variable is
			defined
		nticks: the number of ticks on the colorbar

	Returns:
		ticklabels: Matplotlib ticklabels
		vlim: If catdata is passed then this is reset to match the colormap. None otherwise
		ticklocs: Matplotlib ticklocations
		ncat: if categorical data is found then this is the number of categories. None otherwise
		pointvar: If categorical data is found then the point view is mapped to the colormap values
		cmap: a matplotlib colormap object
		'''
	from . utils import _get_cmap, get_contcbarargs

	if catdata:
		if clim:
			if isinstance(clim, tuple):
				ticklabels = np.arange(clim[0], (clim[1] + 1))
			elif isinstance(clim, list):
				ticklabels = [int(x) for x in clim]
			else:
				raise TypeError('Invalid clim type. Must be tuple or list')
		else:
			ticklabels = np.unique(data[~np.isnan(data.astype(float))])
			ticklabels = [int(x) for x in ticklabels]
		ncat = len(ticklabels)
		if isinstance(catdict, dict) and ncat != len(catdict):
			# modify the cmap to store the colors as if cmap is generated from all cats
			cmap = _get_cmap(cmap, catdata, len(catdict))
			cmap = mpl.colors.ListedColormap(
				[clr for cat, clr in zip(catdict, cmap.colors) if cat in ticklabels]
			)
		else:
			cmap = _get_cmap(cmap, catdata, ncat)
		vlim = (0, ncat)
		ticklocs = np.arange(ncat) + 0.5
		if pointcmap is True:
			dump = pointvar.copy()
			for i in range(ncat):
				dump[pointvar == ticklabels[i]] = i
			pointvar = dump
	else:
		if (pointvar is not None) and (pointcmap is not None):
			try:
				temp_pointdat = pointvar.values()
			except:
				temp_pointdat = pointvar
			temp = np.append(temp_pointdat, view)
			vlim, ticklocs, ticklabels = get_contcbarargs(temp, sigfigs, vlim,
															 nticks=nticks)
		else:
			vlim, ticklocs, ticklabels = get_contcbarargs(view, sigfigs, vlim,
															 nticks=nticks)
		ncat = None
		cmap = _get_cmap(cmap, catdata, ncat)
	return ticklabels, vlim, ticklocs, ncat, pointvar, cmap
