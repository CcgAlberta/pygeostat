#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Basic variogram plotting function using matplotlib reminiscent of varplt from gslib"""

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt


import matplotlib as mpl
import matplotlib.pyplot as plt
from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters
import numpy as np



@set_plot_style
def variogram_plot(data, index=None, sill=1, experimental=True, label=None, ax=None, figsize=None,
		   xlim=None, ylim=None, title=None, xlabel=None, unit=None, ylabel=None, color=None,
		   marker=None, ms=None, ls=None, lw=None, minpairs=40, pairnumbers=False,
		   grid=None, axis_xy=None, plot_style=None, custom_style=None, output_file=None, out_kws=None, **kwargs):
	"""
	This function uses matplotlib to create a variogram plot. Input dataframe structure is
	important as the required data is found within columns that have recognizable headers.

	The only parameter needed is ``data`` and must be a pandas dataframe. All other arguments are
	optional or automatically determined. The settings for experimental and modeled variogram
	plotting is controlled by the ``experimental`` parameter.

	Please review the documentation of the :func:`gs.set_style()
	<pygeostat.plotting.set_style.set_style>` and :func:`gs.export_image()
	<pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
	their use in this function can be understood.

	Parameters:
		data (pd.DataFrame/gs.DataFile): Dataframe/DataFile containing the variogram value, variogram distance, and
			variogram index (if required) data as columns. The dataframe must contain the correct
			column IDs. The column header containing the variogram distance can be: 'h', 'Lag
			Distance', or 'Distance.' The column header containing the variogram values can be:
			'vario', 'Variogram Value', or 'Variogram'
		index (int): Point to which variogram you would like to plot if there are multiple
			variogram within your dataframe. The dataframe must contain the correct column ID. The
			column header containing the variogram index values can be: 'Variogram Index' or
			'Index'
		sill (float): Value to plot a horizontal line representing the variograms sill
		experimental (bool): Indicator if the variogram is experimental ``True`` or modeled
			``False``
		label (str or bool): String to pass to Matplotlib's auto legend function. A default value
			will be generated; however, to prevent this, set label to ``False``
		ax (mpl.axis): Matplotlib axis to plot the figure
		figsize (tuple): Figure size (width, height)
		xlim (float tuple): Minimum and maximum limits of data along the x axis
		ylim (float tuple): Minimum and maximum limits of data along the y axis
		title (str): Title for the plot
		xlabel (str): X-axis label
		unit (str): Distance units used for lag distance. Only used if the keyword parameter
			``xlabel`` is left to its default value of ``None``.
		yalabl (str): Y-axis label
		color (str): Any Matplotlib color
		marker (str): A valid Matplotlib marker style
		ms (float): Marker size in points
		ls (float): A valid Matplotlib line style
		lw (float): Line width in points
		minpairs (int or bool): Any experimental variogram values that were calculated using fewer
			pairs then what is specified by the argument ``minpairs``, is highlighted red. To turn
			this functionality off, set ``minpairs`` to ``False``.
		grid (bool): Plots the major grid lines if True. Based on Parameters['plotting.grid']
			if None.
		axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
			visible) if axis_xy is True. Based on Parameters['plotting.axis_xy'] if None.
		plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
			to turn it off
		custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
		output_file (str): Output figure file name and location
		out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
			:func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
		**kwargs: Optional permissible keyword arguments to pass to matplotlib's plot function

	Returns:
		ax (ax): matplotlib Axes object with the variogram

	**Examples:**

		

		A simple call for experimental variograms, plotting only one direction:

		.. plot::
			
			import pygeostat as gs

			#Load the data from output from varcal/varmodel file
			varcalcdat = gs.ExampleData('experimental_variogram')
			gs.variogram_plot(varcalcdat, index=1)

		|

		A simple call for modeled variograms, plotting only one direction:

		.. plot::

			import pygeostat as gs
			
			varmodeldat = gs.ExampleData('variogram_model')
			gs.variogram_plot(varmodeldat, index=1, experimental=False)


		|

		Plot both experimental and modeled variograms for one direction:

		.. note::
			Some odd behavior may occur if the sill is repeatedly plotted. In the case when
			variograms are being plotted iteratively on the same figure, set the parameter ``sill``
			to ``False`` on all but the last figure.

		.. plot::

			import pygeostat as gs

			varcalcdat = gs.ExampleData('experimental_variogram')
			varmodeldat = gs.ExampleData('variogram_model')
			ax = gs.variogram_plot(varcalcdat.data, index=1, sill=False)
			gs.variogram_plot(varmodeldat.data, index=1, experimental=False, ax=ax)

		|

		Plot both directions experimental and modeled variograms with a legend, grab 2 colors
		from :func:`gs.get_palette() <pygeostat.plotting.utils.get_palette>` to use for the plots,
		and prevent points calculated using a low amount of pairs from being highlighted for one of
		the plots:

		.. plot::

			import pygeostat as gs

			varcalcdat = gs.ExampleData('experimental_variogram')
			varmodeldat = gs.ExampleData('variogram_model')

			colors = gs.get_palette('cat_dark', 2, cmap=False)
			ax = gs.variogram_plot(varcalcdat.data, index=1, color=colors[0], minpairs=False, label=False)
			gs.variogram_plot(varmodeldat.data, index=1, experimental=False, ax=ax, color=colors[0], label='Minor')
			gs.variogram_plot(varcalcdat.data, index=2, ax=ax, color=colors[1], label=False)
			gs.variogram_plot(varmodeldat.data, index=2, experimental=False, ax=ax, color=colors[1], label='Major')

			plt.legend(loc=4)

	"""
	from .export_image import export_image
	from ..data.data import DataFile
	from ..utility.logging import printerr
	from . utils import format_plot, setup_plot, _spatial_aspect

	
	if isinstance (data, DataFile):
		data = data.data

	# Handle dictionary defaults
	if out_kws is None:
		out_kws = dict()
	# Infer the data columns assuming pygeostat or gslib's varcalc function or pygeostat's model
	# function was used
	if 'h' in data.columns:
		h_name = 'h'
	elif 'Lag Distance' in data.columns:
		h_name = 'Lag Distance'
	elif 'Distance' in data.columns:
		h_name = 'Distance'
	else:
		print("Error: A distance column could not be found. Ensure its header is one of the",
			  "following: 'h',\n       'Lag Distance', or 'Distance.'")
		return
	if 'vario' in data.columns:
		vario_name = 'vario'
	elif 'Variogram Value' in data.columns:
		vario_name = 'Variogram Value'
	elif 'Variogram' in data.columns:
		vario_name = 'Variogram'
	else:
		print("Error: A variogram column could not be found. Ensure its header is one of the",
			  "\n       following: 'vario', 'Variogram Value', or 'Variogram.'")
		return
	if index:
		if 'Variogram Index' in data.columns:
			index_name = 'Variogram Index'
		elif 'Index' in data.columns:
			index_name = 'Index'
		else:
			print("Error: A variogram index column could not be found. Ensure its header is one of"
				  " the\n       following: 'Variogram Index' or 'Index.'")
			return
	if minpairs and experimental:
		if experimental or pairnumbers:
			pairs_name = False
			if 'numpairs' in data.columns:
				pairs_name = 'numpairs'
			if 'Number of Pairs' in data.columns:
				pairs_name = 'Number of Pairs'
			if not pairs_name:
				printerr("The argument `minpairs` is being used; however, the column 'numpairs'"
						 " or 'Number of Pairs' was not found in `data`. `minpairs` has been"
						 " set to a value of `False`.", errtype='warning')
				minpairs = False
		else:
			pairs_name = False
			minpairs = False
			pairnumbers = False
	else:
		pairs_name = False
		minpairs = False
		pairnumbers = False
	# Set-up some default plotting parameters
	if color is None:
		color = Parameters['plotting.variogram_plot.color']
	if experimental:
		if ms is None:
			ms = Parameters['plotting.variogram_plot.ms']
		if marker is None:
			marker = 'o'
		if lw is None:
			lw = 0
		if ls is None and lw > 0:
			ls = '-'
		else:
			ls = 'None'
	else:
		if lw is None:
			lw = 1
		if ls is None:
			ls = '-'
		if ms is None:
			ms = 0
	if label is None:
		if index:
			label = 'Var %s' % index
		else:
			label = 'Var'
	elif not label:
		label = "_nolegend_"
	# Get the right data
	data = data[data[vario_name] != np.nan]
	# This shouldn't be necessary once everything is converted to an nan standard, but for now..
	data = data[data[vario_name] != Parameters['data.null']]
	if experimental:
		data = data[data[h_name] > 0]
	pairs_label = None
	varpairs = None
	if index is None:
		vardata = data[vario_name]
		vardist = data[h_name]
		if minpairs or pairnumbers and experimental:
			varpairs = data[data[pairs_name] <= minpairs]
			pairs_label = data[pairs_name]
	else:
		vardata = data[vario_name][data[index_name] == index]
		vardist = data[h_name][data[index_name] == index]
		if minpairs or pairnumbers and experimental:
			varpairs = data[(data[pairs_name] <= minpairs) & (data[index_name] == index)]
			pairs_label = data[pairs_name]
		if len(vardata) == 0 or len(vardist) == 0:
			print("Error: The index value provided wasn't found")
			return
	# Set-up plot if no axis is supplied
	fig, ax, cax = setup_plot(ax, figsize=figsize, aspect =False)
	# Highlight the pairs with fewer pairs then minpairs if required
	if minpairs:
		for i in range(len(varpairs)):
			try:
				ax.plot(varpairs.iloc[i][h_name], varpairs.iloc[i][vario_name], marker=marker,
						color='#fb8072', ms=round((ms * 1.4 + 0.5), 0))
			except:
				ax.plot(varpairs.iloc[i][h_name], varpairs.iloc[i][vario_name], marker=marker,
						ms=round((ms * 1.4 + 0.5), 0))
	# Plot the variogram and identifier
	ax.plot(vardist, vardata, color=color, lw=lw, marker=marker, ls=ls, ms=ms, label=label,
			**kwargs)
	if experimental and pairnumbers and pairs_label is not None:
		if ylim is None:
			ymin, ymax = ax.get_ybound()
		else:
			ymin, ymax = ylim
		for x, y, lbl in zip(vardist, vardata, pairs_label):
			if ymin <= y <= ymax:
				ax.text(x, y - 0.025, int(lbl), rotation=-55)
	# Plot figure and axis labels
	if unit is None:
		unit = Parameters['plotting.unit']
	if unit is None or unit == '':
		unit = ''
	else:
		unit = ' ({})'.format(unit)
	if xlabel is None:
		xlabel = Parameters['plotting.lagname'] + unit
	if ylabel is None:
		ax.set_ylabel(r'$\gamma$   ', fontsize=mpl.rcParams['font.size'] *
					  Parameters['plotting.gammasize'], rotation=0)
	elif ylabel:
		ax.set_ylabel(ylabel)
	# Draw Sill
	if sill:
		ax.axhline(sill, color='black')
	# Set plot limits
	if ylim is None:
		if sill == 1:
			ylim = (0, 1.2)
		elif abs(sill) < 1:
			# If the sill is bellow 1, it's likely a cross variogram
			if sill < 0:
				# Positive default which also draws a line at 0
				ylim = (-1, 0.2)
				ax.axhline(0, color='black', linestyle='--', lw=0.5, dashes=(2, 2))
			if sill >= 0:
				# Positive defaults
				ylim = (0, 1.2)
	if xlim is None:
		xmin, xmax = ax.get_xbound()
		xlim = [0, xmax]

	ax = format_plot(ax, xlabel=xlabel, title=title, grid=grid, axis_xy=axis_xy, xlim=xlim, ylim=ylim)

	# Export figure
	if output_file or ('pdfpages' in out_kws):
		export_image(output_file, **out_kws)
	# Return axis
	return ax
