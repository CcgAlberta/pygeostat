#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""drill_plot for a drill hole data set"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from . set_style import set_plot_style

@set_plot_style
def drill_plot(z, var, cat=None, categorical_dictionary=None, lw=2, line_color='k', barwidth=None, namelist=None,
			   legend_fontsize=10, title=None, ylabel=None, unit=None, grid=None, axis_xy=None,
			   reverse_y=False, xlim=None, ylim=None,  figsize=(2,10), ax=None,plot_style=None,
			   custom_style=None, output_file=None, out_kws=None, **kwargs):
	'''
	 A well log plot for both continuous and categorical variables. This plot handles one well log
	 plot at a time and the user can choose to generate subplots and pass the axes to this function
	 if multiple well log plots are required.

	Parameters:
		z (Elevation/Depth or distance along the well): Tidy (long-form) 1D data where a single
			column of the variable exists with each row is an observation. A pandas dataframe/series
			or numpy array can be passed.
		var (Variable being plotted): Tidy (long-form) 1D data where a single column of the variable
			exists with each row is an observation. A pandas dataframe/series or numpy array can be
			passed.
		lw (float): line width for log plot of a continuous variable
		line_color (string): line color for the continuous variable
		barwidth(float): width of categorical bars
		categorical_dictionary (dictionary): a dictionary of colors and names for categorical codes. E.g.
		{1: {'name': 'Sand', 'color': 'gold'},
		 2: {'name': 'SHIS','color': 'orange'},
		 3: {'name': 'MHIS','color': 'green'},
		 4: {'name': 'MUD','color': 'gray'}}
		legend_fontsize(float): fontsize for the legend plot rleated to the categorical codes. set
			this parameter to 0 if you do not want to have a legend
		title (str): title for the variable
		ylabel (str): Y-axis label, based on ``Parameters['plotting.zname']`` if None.
		unit (str): Unit to place inside the y-label parentheses, based on
			``Parameters['plotting.unit']`` if None.
		grid (bool): Plots the major grid lines if True. Based on ``Parameters['plotting.grid']`` if
			None.
		axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
			visible) if axis_xy is True. Based on ``Parameters['plotting.axis_xy']`` if None.
		reverse_y(bool): if true, the yaxis direction is set to reverse(applies to the cases that
			depth is plotted and not elevation)
		aspect (str): Set a permissible aspect ratio of the image to pass to matplotlib.
		xlim (float tuple): X-axis limits
		ylim (float tuple): Y-axis limits
		figsize (tuple): Figure size (width, height)
		ax (mpl.axis): Existing matplotlib axis to plot the figure onto
		out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
			:func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`

	Returns:
		ax (ax): Matplotlib axis instance which contains the gridded figure

	Examples:

		A simple call to plot a continuous variable

	.. plot:: 

		import pygeostat as gs
		dat = gs.ExampleData('point3d_ind_mv')
		data = dat.data[dat.data['HoleID'] == 3] 
		gs.drill_plot(data['Elevation'], data['Sw'], grid = True)

	|

	Plot a categorical variable

	.. plot:: 

		import pygeostat as gs
		dat = gs.ExampleData('point3d_ind_mv')
		data = dat.data[dat.data['HoleID'] == 3] 
		gs.drill_plot(data['Elevation'], data['Lithofacies'])

	|

	Plot a categorical variable and provide a categorical dictionary
	
	.. plot:: 

		import pygeostat as gs
		dat = gs.ExampleData('point3d_ind_mv')
		data = dat.data[dat.data['HoleID'] == 3]
		cat_dict = {1: {'name': 'Sand', 'color': 'gold'},
			3: {'name': 'SHIS','color': 'orange'},
			4: {'name': 'MHIS','color': 'green'},
			5: {'name': 'MUD','color': 'gray'}}
		gs.drill_plot(data['Elevation'], data['Lithofacies'], categorical_dictionary=cat_dict)

	'''

	from matplotlib.ticker import FormatStrFormatter
	import matplotlib.patches as mpatches
	import matplotlib.pyplot as plt
	from .export_image import export_image
	import numpy as np
	import pandas as pd
	from .. pygeostat_parameters import Parameters
	from .utils import format_plot, get_label, setup_plot, catcmapfromcontinuous, _spatial_aspect

	# Sanity checks
	if not isinstance(z, (np.ndarray, np.generic, pd.Series)):
		raise ValueError("The parameter `z` must be a numpy array or a pandas series")
	if not isinstance(var, (np.ndarray, np.generic, pd.Series)):
		raise ValueError("The parameter `var` must be a numpy array or a pandas series")
	if title is None:
		title = get_label(var)
	if unit is None:
		unit = Parameters['plotting.unit']
	if unit is None or unit == '':
		unit = ''
	else:
		unit = ' ({})'.format(unit)
	if ylabel:
		ax.set_ylabel(ylabel)
	elif ylabel is None:
		ylabel = Parameters['plotting.zname'] + unit
	if not out_kws:
		out_kws = dict()

	fig, ax, cax = setup_plot(ax, figsize=figsize, aspect = False)

	# Create a new data frame
	df = pd.DataFrame()
	df['Z'] = np.array(z)
	df['Var'] = np.array(var)
	df.sort_values('Z', axis=0, ascending=True, inplace=True)
	df.reset_index(inplace=True, drop=True)


	if barwidth is None:
		barwidth = (np.max(df['Z'].values)  - np.min(df['Z'].values))/10
	
	if cat is None:
		if (len(set(df['Var'].values)) < (len(df) / 5)) :
			cat = True

	# Start of the plot
	if not cat:

		ax.plot(df['Var'], df['Z'], lw=lw, c=line_color, **kwargs)
		xlim = [np.nanmin(df['Var']), np.nanmax(df['Var'])]
		ylim = [np.nanmin(df['Z']), np.nanmax(df['Z'])]
		ax.set_xticks(np.linspace(np.nanmin(df['Var']),
								  np.nanmax(df['Var']), 2))
		ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
		ax.xaxis.set_label_position('top')
		ax.xaxis.set_ticks_position('top')
	else:
		categories = sorted(set(df['Var'].values))
		ncat = len(categories)
		
		# make sure grid option is not triggered
		grid = False

		if categorical_dictionary is None:
			colorlist = catcmapfromcontinuous("Spectral", ncat).colors
			namelist = []
			for i in range(ncat):
				string = 'Category Code %i' % (i + 1)
				namelist.append(string)
		
			categorical_dictionary = {}
			for i, item in enumerate(categories):
				categorical_dictionary.update({item: {'name': namelist[i], 'color': colorlist[i]} })
		else:
			if len(categorical_dictionary)<ncat:
				raise ValueError('The provided categorical_dictionary does not match the number of categorical codes')
			for item in categorical_dictionary.keys():
				if item not in categories:
					raise ValueError('The provided categorical code, {} does not exist in the provided data'.format(item))

		
		# First bar from bottom based on sorting the z column
		category = df['Var'][0]
		color = categorical_dictionary[category]['color']
		offset = (df['Z'][1] - df['Z'][0])/2
		z_b = df['Z'][0]
		z_t = z_b + offset
		ax.fill_between([0, barwidth], [z_b, z_b], [z_t,z_t], color = color)

		for i in range(1,len(df)-1):
			upper_offset = (df['Z'][i+1] - df['Z'][i])/2
			lower_offset = (df['Z'][i-1] - df['Z'][i])/2
			z_b = df['Z'][i] + lower_offset
			z_t = df['Z'][i] + upper_offset
			category = df['Var'][i]
			color = categorical_dictionary[category]['color']
			ax.fill_between([0, barwidth], [z_b, z_b], [z_t,z_t], color = color)

		# Last bar from bottom based on sorting the z column
		category = df['Var'][i+1]
		color = categorical_dictionary[category]['color']
		z_b = z_t
		z_t = df['Z'][i+1]
		ax.fill_between([0, barwidth], [z_b, z_b], [z_t,z_t], color = color)

		Patch_list = []
		for item in categorical_dictionary.values():
			Patch_list.append(mpatches.Patch(color=item['color'], label=item['name']))

		if legend_fontsize > 0:
			ax.legend(handles=Patch_list, loc='upper left',
					  bbox_to_anchor=(1, 1), fontsize=legend_fontsize)

		ylim = [np.nanmin(df['Z']), np.nanmax(df['Z'])]
		xlim = [0, barwidth]
		ax.xaxis.set_label_position('top')
		ax.get_xaxis().set_visible(False)

	format_plot(ax, xlabel='', ylabel = ylabel, title = title, grid = grid, axis_xy = axis_xy, xlim = xlim, ylim = ylim)
	
	if reverse_y:
		ax.invert_yaxis()

	if output_file or ('pdfpages' in out_kws):
		export_image(output_file, **out_kws)

	return ax
