#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""Location map plotting routine using matplotlib"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters

@set_plot_style
def location_plot(data, x=None, y=None, z=None, var=None, dhid = None, catdata=None, allcats=True, cbar=True,
		   cbar_label=None, catdict=None, cmap=None, cax=None, vlim=None, title=None,  plot_collar = True, collar_marker='x',
		   collar_offset = 0, lw = None, xlabel=None, ylabel=None, unit=None, griddef=None, slice_number=None, orient='xy',
		   slicetol=None, xlim=None, ylim=None, ax=None, figsize=None, s=None, marker='o',
		   rotateticks=None, sigfigs=3, grid=None, axis_xy=None, aspect=None, plot_style=None,
		   custom_style=None, output_file=None, out_kws=None, return_cbar=False, return_plot=False,
		   **kwargs):
	"""
	location_plot displays scattered data on a 2-D XY plot. To plot gridded data with or without
	scattered data, please see :func:`gs.slice_plot()<pygeostat.plotting.slice_plot>`.

	The only required parameter is ``data`` if it is a
	:class:`gs.DataFile <pygeostat.data.data.DataFile>` that contains the necessary coordinate
	column headers, data, and if required, a pointer to a valid
	:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>` class. All other parameters are
	optional. If ``data`` is a :class:`gs.DataFile <pygeostat.data.data.DataFile>` class and does
	not contain all the required parameters or if it is a long-form table, the following
	parameters will need to be pass are needed: ``x``, ``y``, ``z``, and ``griddef``. The three
	coordinate parameters may not be needed depending on what ``orient`` is set to and of course
	if the dataset is 2-D or 3-D. The parameter ``griddef`` is required if ``slicetol`` or
	`` slice_number`` is used. If parameter ``slice_number`` and ``slicetol`` is not set then the default
	slice tolerance is half the cell width. If a negative ``slicetol`` is passed or slice_number is set
	to None then all data will be plotted. ``slicetol`` is based on coordinate units.

	The values used to bound the data (i.e., vmin and vmax) are automatically calculated by default.
	These values are determined based on the number of significant figures and the sliced data;
	depending on data and the precision specified, scientific notation may be used for the colorbar
	tick lables. When point data shares the same colormap as the gridded data, the points displayed
	are integrated into the above calculation.

	Please review the documentation of the :func:`gs.set_style()
	<pygeostat.plotting.set_style.set_style>` and :func:`gs.export_image()
	<pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
	their use in this function can be understood.

	Parameters:
		data (pd.DataFrame or gs.DataFile): data containing coordinates and (optionally) var
		x (str): Column header of x-coordinate. Required if the conditions discussed above are not
			met
		y (str): Column header of y-coordinate. Required if the conditions discussed above are not
			met
		z (str): Column header of z-coordinate. Required if the conditions discussed above are not
			met
		var (str): Column header of the variable to use to colormap the points. Can also be a list
			of or single permissible matplotlib colour(s). If None and data is a DataFile,
			based on DataFile.variables if len(DataFile.variables) == 1. Otherwise, based on
			Parameters['plotting.location_plot.c']
		dhid (str): Column header of drill hole ID.
		catdata (bool): Force categorical data
		catdict (dict): Dictionary containing the enumerated IDs alphabetic equivalent, which is
			drawn from Parameters['data.catdict'] if None
		allcats (bool): ensures that if categorical data is being plotted and plotted on slices,
			that the categories will be the same color between slices if not all categories are
			present on each slice
		cbar (bool): Indicate if a colorbar should be plotted or not
		cbar_label (str): Colorbar title
		cmap (str): A matplotlib colormap object or a registered matplotlib or pygeostat colormap
			name.
		cax(Matplotlib.ImageGrid.cbar_axes): color axis, if a previously created one should be used
		vlim (float tuple): Data minimum and maximum values
		title (str): Title for the plot. If left to it's default value of ``None`` or is set to
			``True``, a logical default title will be generated for 3-D data. Set to ``False`` if
			no title is desired.
		plot_collar (bool): Option to plot the collar if the orient is xz or yz and the dhid is provided/inferred
		collar_marker (str): One of the permissible matplotlib markers, like 'o', or '+'... and others.
		lw (float): Line width value if the orient is xz or yz and the dhid is provided/inferred. Because lines or plotted instead of points. 
		xlabel (str): X-axis label
		ylabel (str): Y-axis label
		unit (str): Unit to place inside the axis label parentheses
		griddef (GridDef): A pygeostat GridDef class created using
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Required if using the
			argument ``slicetol``
		orient (str): Orientation to slice data. ``'xy'``, ``'xz'``, ``'yz'`` are t he only accepted
			values
		slice_number (int): Grid cell location along the axis not plotted to take the slice of data to
			plot. None will plot all data
		slicetol (float): Slice tolerance to plot point data (i.e. plot +/- ``slicetol`` from the
			center of the slice). Any negative value plots all data. Requires ``slice_number``. If a
			``slice_number`` is passed and no ``slicetol`` is set, then the default will half the cell
			width based on the griddef.
		xlim (float tuple): X-axis limits. If None, based on data.griddef.extents(). If
			data.griddef is None, based on the limits of the data.
		ylim (float tuple): Y-axis limits. If None, based on data.griddef.extents(). If
			data.griddef is None, based on the limits of the data.
		ax (mpl.axis): Matplotlib axis to plot the figure
		figsize (tuple): Figure size (width, height)
		s (float): Size of location map markers
		marker (str): One of the permissible matplotlib markers, like 'o', or '+'... and others.
		rotateticks (bool tuple): Indicate if the axis tick labels should be rotated (x, y)
		sigfigs (int): Number of sigfigs to consider for the colorbar
		grid (bool): Plots the major grid lines if True. Based on Parameters['plotting.grid']
			if None.
		axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
			visible) if axis_xy is True. Based on Parameters['plotting.axis_xy'] if None.
		aspect (str): Set a permissible aspect ratio of the image to pass to matplotlib. If None,
			it will be 'equal' if each axis is within 1/5 of the length of the other. Otherwise,
			it will be 'auto'.
		plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
			to turn it off
		custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
		output_file (str): Output figure file name and location
		out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
			:func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
		return_cbar (bool): Indicate if the colorbar axis should be returned. A tuple is returned
			with the first item being the axis object and the second the cbar object.
		return_plot (bool): Indicate if the plot from scatter should be returned. It can be used
							to create the colorbars required for subplotting with the ImageGrid()
		**kwargs: Optional permissible keyword arguments to pass to matplotlib's scatter function

	Returns:
		ax (ax): Matplotlib axis instance which contains the gridded figure

	Returns:
		cbar (cbar): Optional, default False. Matplotlib colorbar object

	**Examples:**

	A simple call:

	.. plot::

		import pygeostat as gs
		data_file = gs.ExampleData('point3d_ind_mv')
		gs.location_plot(data_file)

	|

	A simple call using a variable to color the data:

	.. plot::

		import pygeostat as gs
		data_file = gs.ExampleData('point3d_ind_mv')
		gs.location_plot(data_file, var = 'Phi')

	|

	Plotting along xz/yz orientation and use line plots based on drill hole id

	.. plot::

		import pygeostat as gs
		# load data
		data_file = gs.ExampleData('point3d_ind_mv')
		gs.location_plot(data_file, var='Phi', orient='yz', aspect =5, plot_collar = True)
		
	|

	Location plot for a categorical variable

	.. plot::

		import pygeostat as gs
		# load data
		data_file = gs.ExampleData('point3d_ind_mv')
		gs.location_plot(data_file, var='Lithofacies', orient='yz', aspect =5, plot_collar = True)

	|
	
	"""
	from . utils import (_spatial_labels, _spatial_pointdata, _spatial_slice, _format_tick_labels, setup_plot,
						 _spatial_orient2fig, format_plot, _get_cmap, _spatial_aspect, get_contcbarargs)

	from .. pygeostat_parameters import Parameters
	from . export_image import export_image
	from ..data.data import DataFile
	from .cmaps import avail_palettes
	from matplotlib.collections import LineCollection
	
	# Handle dictionary defaults
	if not out_kws:
		out_kws = dict()

	# Infer dhid column if not provided
	drill_plot_format = False
	if dhid is None:
		if isinstance(data, DataFile):
			dhid = data.dh
	else:
		if not dhid in data.columns:
			raise ValueError('The provided drill hole id column, {} does not exist'.format(dhid))

	# Handle var
	if var is None:
		if isinstance(data, DataFile):

			if data.variables is not None:
				if isinstance(data.variables, str):
					var = data.variables
			elif isinstance(data.cat, str):
				var = data.cat
		if var is None:
			var = Parameters['plotting.location_plot.c']
	
	# Handle title
	if title is None:
		if var in data.columns:
			title = var
	# Determine the x, y, z and griddef with error checking
	data, x, y, z, griddef = _spatial_pointdata(data, orient, x, y, z, griddef)

	# Determine the data to plot based on grid def and slice info
	if dhid is None:
		pointx, pointy, pointvar = _spatial_slice(data, var, x, y, z, griddef, orient, slice_number, slicetol)
	else:
		pointx, pointy, pointvar, point_dh = _spatial_slice(data, var, x, y, z, griddef, orient, slice_number, slicetol, dhid)

	if pointx is None:
		if griddef is None:
			raise ValueError('At least two coordinates are required for the location_plot')
		else:
			print('Note: There is no data point detected within the provided slice number and/or slice tolerance')

	# Setup the color
	if pointvar is None:
		try:
			pointvar = mpl.colors.ColorConverter().to_rgb(var)
			if int(mpl.__version__.split('.')[0]) >= 3:
				pointvar = [pointvar]
			cbar = False
		except Exception:
			if griddef is None:
				raise ValueError(
					'var is not a variable in the passed data and not a valid Matplotlib color {}'.format(var))
			else:
				pointvar = [(0.0, 0.0, 0.0)] # Use a default color since there is no pint within the grid
				cbar = False

	# Get default colormaps or palettes. If there are Parameters['plotting.assumecat'] or less
	# unique data, assume categorical
	if isinstance(var, str):
		if var in data.columns:
			if allcats:
				alldata = data[var]
			else:
				alldata = pointvar
			ncat = len(np.unique(alldata[np.isfinite(alldata)]))
		else:
			ncat = None
	else:
		ncat = None
		cmap = None
	if catdata is not True:
		if cmap is None and catdata is None and ncat is not None and \
				ncat < Parameters['plotting.assumecat']:
			catdata = True
		# Get palette from pygeostat
		elif cmap in avail_palettes:
			catdata = True
		elif ncat is not None and ncat <= Parameters['plotting.assumecat'] and catdata is None:
			catdata = True
		else:
			catdata = False
	if cmap is None:
		cmap = _get_cmap(cmap, catdata, ncat)
	# Set-up categorical parameter if required
	if catdata:
		ticklabels = np.unique(alldata[np.isfinite(alldata)]).astype(int)
		if catdict is None:
			catdict = Parameters['data.catdict']
		if isinstance(catdict, dict):
			if len(ticklabels) != len(catdict):
				# modify the cmap to store the colors as if cmap is generated from all cats
				cmap = _get_cmap(cmap, True, len(catdict))
				cmap = mpl.colors.ListedColormap(
					[clr for cat, clr in zip(catdict, cmap.colors) if cat in ticklabels]
				)
		vlim = (0, ncat)
		ticklocs = np.arange(ncat) + 0.5
		if cmap:
			if isinstance(cmap, str):
				cmap = _get_cmap(cmap, catdata, ncat)
			dump = pointvar.copy()
			for i in range(ncat):
				dump[pointvar == ticklabels[i]] = i
			pointvar = dump
		if isinstance(catdict, dict):
			ticklabels = [catdict[cat] for cat in ticklabels]
	# Set-up some parameters
	if not catdata and cmap is not False and cbar:
		vlim, ticklocs, ticklabels = get_contcbarargs(pointvar, sigfigs, vlim)
	if vlim is None:
		vlim = (None, None)
	# Set-up plot if no axis is supplied using the ImageGrid method if required or the regular way
	fig, ax, cax = setup_plot(ax, cax=cax, cbar=cbar, figsize=figsize)
	if s is None:
		s = Parameters['plotting.location_plot.s']

	# used for plot extents and get the two coordinates being plotted
	figx, figy, _ = _spatial_orient2fig(orient, x, y, z)

	if orient.lower() == 'xy' or not dhid:
		if pointvar is not None and len(pointvar) > 0:
			plot = ax.scatter(pointx, pointy, c=pointvar, cmap=cmap, vmin=vlim[0],
					 vmax=vlim[1], s=s, marker=marker, **kwargs)
		else:
			plot = ax.scatter(pointx, pointy)
	else: # Plot lines if the dhid was inferred/provided and the orientation is xz or yz
		drill_plot_format = True
		data_temp = pd.DataFrame({figx:pointx, figy: pointy, var:pointvar, dhid: point_dh})
		grouped = data_temp.groupby(dhid, sort=False)

		if lw is None:
			lw = Parameters['plotting.location_plot.lw']
		for _, group in grouped:
			if plot_collar:
				x_points = group[figx][group[figx].index.min()]
				y_points = group[figy][group[figy].index.min()] + collar_offset
				if var in data.columns:
					collar_var = group[var][group[var].index.min()]
					norm = plt.cm.colors.Normalize(*vlim)
					_cmap = plt.colormaps.get_cmap(cmap)
					collar_var = [_cmap(norm(collar_var))]
					plot = ax.scatter(x_points, y_points, c=collar_var, cmap=cmap, vmin=vlim[0], vmax=vlim[1], s=s,
							marker=collar_marker, **kwargs)
				else:
					plot = ax.scatter(x_points, y_points, c=pointvar, vmin=vlim[0], vmax=vlim[1], s=s,
							marker=collar_marker, **kwargs)

			points = np.array([group[figx], group[figy]]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)

			if var in data.columns:
				lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(vlim[0], vlim[1]))
				lc.set_array(group[var])
			else:
				lc = LineCollection(segments, color=var, norm=plt.Normalize(vlim[0], vlim[1]))
			lc.set_lw(lw)

			ax.add_collection(lc)

	# Plot labels
	ax = _spatial_labels(ax, orient, griddef, slice_number, title, xlabel, ylabel, unit, sigfigs)
	
	
	if xlim is None:
		if griddef:
			if orient in ['xy', 'xz']:
				xlim = griddef.extents()[0]
			else:
				xlim = griddef.extents()[1]
		else:
			pad = (data[figx].max() - data[figx].min()) * 0.025
			xlim = (data[figx].min() - pad, data[figx].max() + pad)
	if ylim is None:
		if griddef:
			if orient == 'xy':
				ylim = griddef.extents()[1]
			else:
				ylim = griddef.extents()[2]
		else:
			pad = (data[figy].max() - data[figy].min()) * 0.025
			ylim = (data[figy].min() - pad, data[figy].max() + pad)
	if aspect is None:
		aspect = _spatial_aspect(xlim, ylim)

	if aspect:
		ax.set_aspect(aspect)  # Plot the figure
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
	ax = _format_tick_labels(ax, rotateticks)

	# Plot colorbar, the tick locations and labels were generated above using get_vlimdata
	if cbar:
		# Set-up axis to place colorbar into
		# cax = grid.cbar_axes[0]
		# Plot the colorbar
		if drill_plot_format:
			cbar = lc.axes.figure.colorbar(lc, cax=cax, ticks=ticklocs)
		else:
			cbar = plot.figure.colorbar(plot, cax=cax, ticks=ticklocs)
		# Configure the color bar
		cbar.ax.set_yticklabels(ticklabels, ha='left')
		cbar.ax.tick_params(axis='y', pad=2)
		if cbar_label is not None:
			cbar.set_label(cbar_label, ha='center', va='top', labelpad=2)
	if axis_xy is None:
		axis_xy = Parameters['plotting.axis_xy_spatial']
	format_plot(ax, grid=grid, axis_xy=axis_xy, xlim=xlim, ylim=ylim)
	# Export figure
	if output_file or ('pdfpages' in out_kws):
		export_image(output_file, **out_kws)
	# Return whats needed
	if return_cbar:
		return ax, cbar
	if return_plot:
		return plot
	return ax
