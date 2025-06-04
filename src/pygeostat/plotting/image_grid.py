#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""image_grid.py: Provides a means to easily generate grids of plots"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from . set_style import set_plot_style


@set_plot_style
def image_grid(ncol, nrow, symmetric=False, nvar=None, gridfunc=None, upperfunc=None,
			  lowerfunc=None, diagfunc=None, tight=False, axes_pad=0.2, cbar=False, vlim=None,
			  cbar_label=None, figsize=None, aspect=None, xlim=None, ylim=None, axislabels=None,
			  xlabel=None, ylabel=None, suptitle=None, ntickbins=2, rotateticks=False,
			  plot_style=None, custom_style=None, output_file=None, labelmode='all',
			  unequal_aspects=False, direction='row', **kwargs):
	"""
	Create either a symmetric or non-symmetric plot matrix. This function interprets a symmetric
	plot matrix as one that plots multivariate data in different types of bivariate plots in the
	lower and/or upper triangles and univariate plots along the diagonal. Conversely, a non-
	symmetric plot matrix assumes a single plotting function is used to populate all of the
	subplots.

	To provide a very flexible plotting skeleton, this function does not actually instruct any
	plotting functions. Instead, the user is required to define python generators that contain a
	loop that will produce the desired subplots. Use of generators allows the user to
	customize the subplots as desired, without having to use this function as a middle man.

	If plotting a symmetric plot matrix, the keyword argument ``nvar`` is required and one of the
	following is required: ``upperfunc``, ``owerfunc``, or ``diagfunc``. If plotting a
	non-symetric plot matrix, the following keyword arguments are required: ``nvar``, ``ncol``,and
	``gridfunc``.

	.. note::

		The subplots regardless of their location, plot left to right, top to bottom.

		The plots aspect will need to be considered. If your plots appear flat, your aspect is
		wrong. As an example, variograms are classically plotted at a 4:3 ratio, meaning if
		the y-axis limits are left to their default of 0 to 1.2 and say your x-axis is being
		plotted to 1000, you would require a aspect of :math:`{1000 / 1.2 / (4/3)}` or 62.5. If
		this function detects a y-axis limit of 1.2, it will calculate an aspect automatically
		unless ``aspect`` is set manually.

	Colorbars can be plotted for each individual subplot by setting keyword argument ``cbar`` to
	``'each'``; however, with this setting, this function does not handle any colorbar plotting.
	Therefore, the color bar axes are passed to the iterator and the plotting function within the
	iterator deals with the colorbar. If a single colorbar is desired, ``cbar`` is set to
	``'single'``. When using this setting, the keyword argument ``vlim`` must be passed. It is also
	important that all of the subplots have their colormaps limited to this same range. If no color
	bar(s) are desired, the keyword argument is set to ``'none'``.

	Please review the documentation of the
	:func:`gs.set_style() <pygeostat.plotting.set_style>` and
	:func:`gs.export_image() <pygeostat.plotting.export_image>` functions for details on their
	parameters so that their use in this function can be understood.

	Keyword Arguments:
		symmetric (bool):
		nvar (int): Only used for symmetric grids. Number of to variables to plot
		ncol (int): Only used for non-symmetric grids. Number of to columns to plot
		nrow (int): Only used for non-symmetric grids. Number of to rows to plot
		gridfunc (generator): Python generator that contains a loop that can be used to plot the
			desire subplots for the whole grid
		upperfunc (generator): Python generator that contains a loop that can be used to plot the
			desire subplots in the upper triangle of the grid
		lowerfunc (generator): Python generator that contains a loop that can be used to plot the
			desire subplots in the lower triangle of the grid
		diagfunc (generator): Python generator that contains a loop that can be used to plot the
			desire subplots along the diagonal of the grid
		tight (bool): Indicate if the whitespace between the subplots should be removed.
			Reminiscent of R's scatter
		axes_pad (float or tuple): Padding in iches to place between the plots. Can pass a tuple to
			indicate different padding in the horizontal and veritcal directions (width_pad,
			height_pad)
		cbar (str): Indicate what colorbar mode to use. The available options are ``['none',
			'single', 'each']``. See above for instructions
		vlim (tuple): If ``cbar`` is set to ``'single'``, this value is required and instructs the
			limits of the colorbar.
		cbar_label (str): Colorbar title
		figsize (tuple): Figure size (width, height)
		aspect (str): Set a permissible aspect ratio of the image to pass to matplotlib. The
			function will try and detect what aspect is best as described above.
		xlim (float tuple): X-axis limits applied to all axes in the grid
		ylim (float tuple): Y-axis limits applied to all axes in the grid
		axislabels (list): Only used for symmetric grids. Labels for each row and column
		xlabel (str): Super x-axis label
		yalabl (str): Super y-axis label
		suptitle (str): Super title
		ntickbins (int or tuple): int: applied to both x and y, or tuple, applied to x and y
			respectively
		rotateticks (bool or float tuple): Bool or float values to control axis label rotations.
			See above for more info.
		plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
			:class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
			to turn it off
		custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
		output_file (str): Output figure file name and location
		labelmode (str): Labeling input to the image_grid function. Default is ``'all'``. ``'L'``
			labels the left column and bottom row only. There may be other valid parameters. The
			last plot in each column will get x-axis labels in ``'L'`` mode as unused subplots in
			the last row are removed. If ``tight`` is set to ``True``, ``'L'`` is always used.
		unequal_aspects (bool): Indicate `True` if the limits of the plots will differe between
			subplots. This will then use plt.subplots() rather than image_grid() for plotting
			functions
		**kwargs: Optional permissible keyword arguments to pass to
			:func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`

	Examples:

		Below are a few images illustrating the results of this function.

		Accuracy plot matrix using the results from multivariate simulation:

		.. image:: ./figures/advancedplotting_imggrid_accplt.png

	|
	"""
	import types
	import numpy as np
	from .export_image import export_image
	from .utils import get_contcbarargs, get_supaxislocs
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	from mpl_toolkits.axes_grid1 import ImageGrid, LocatableAxes, SubplotDivider

	# Sanity checks
	if symmetric:
		if nvar is None:
			raise ValueError("When plotting a symmetric grid, the keyword argument `nvar` must be"
							 " passed")
		if all([upperfunc is None, lowerfunc is None, diagfunc is None]):
			raise ValueError("When plotting a symmetric grid, one of the keyword arguments"
							 "`upperfunc`, `lowerfunc`, or `diagfunc` must be passed")
	else:
		if any([ncol is None, nrow is None]):
			raise ValueError("When plotting a non-symmetric grid, the keyword arguments `ncol` and"
							 " `nrow` must be pased")
		if gridfunc is None:
			raise ValueError("When plotting a non-symmetric grid, the keyword argument `gridfunc`"
							 " must be passed")
	# ------------------------------------------------
	# Setup Styles, Parameters, and Image Grid
	# ------------------------------------------------
	# Sort out some figure parameters
	if symmetric:
		if ncol is None:
			ncol = nvar
		if nrow is None:
			nrow = nvar
	if tight:
		labelmode = 'L'
		axes_pad = 0
		fontsize = mpl.rcParams['font.size']
		if fontsize > 8:
			scalar = (0.42 * (fontsize / 16)) + 1
		elif fontsize < 8:
			scalar = (0.83 * (fontsize / 5))
		else:
			scalar = 1
		cbar_pad = 0.2 + 0.18 * scalar
		cbar_location = 'bottom'
	else:
		cbar_pad = None
		cbar_location = 'right'
	if cbar in ['single', 'each']:
		if cbar == 'each':
			cbar_pad = 0.05
			cbar_size = 0.075
		else:
			cbar_size = 0.075
	else:
		cbar_size = None
	# ------------------------------------------------
	# Setup the axes either using plt.subplots() or image_grid based on `unequal_aspects`
	# ------------------------------------------------
	if unequal_aspects:
		# checks for unhandled combinations
		if tight:
			raise ValueError('`tight=True` and `unequal_aspects=True` cannot be used together '
							 'with `gs.image_grid()`')
		if cbar:
			raise ValueError('Cannot specify `cbar` with `unequal_aspects=True`, use the `cbar` '
							 'plotting options in the respective iterator functions instead')
		# setup this function to use plt.subplots
		fig, grid = plt.subplots(nrow, ncol, figsize=figsize, facecolor='none')
		ax_all = grid.flatten().tolist()
		if isinstance(axes_pad, tuple):
			plt.subplots_adjust(wspace=axes_pad[0], hspace=axes_pad[1])
		elif isinstance(axes_pad, (int, float)):
			plt.subplots_adjust(wspace=axes_pad, hspace=axes_pad)
		else:
			plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
	else:
		# Set-up Figure object and image grid
		fig = plt.figure(figsize=figsize)
		grid = ImageGrid(fig, 111, nrows_ncols=(nrow, ncol), axes_pad=axes_pad, cbar_mode=cbar,
						 cbar_size=cbar_size, cbar_pad=cbar_pad, cbar_location=cbar_location,
						 share_all=False, label_mode=labelmode, direction=direction)
		ax_all = grid.axes_all
	# Fix some figure aesthetics prior to the plotting function so the user is able to move away
	# from these defaults if desired
	for ax in ax_all:
		if isinstance(ntickbins, tuple):
			ax.locator_params(axis='x', nbins=ntickbins[0])
			ax.locator_params(axis='y', nbins=ntickbins[1])
		else:
			ax.locator_params(nbins=ntickbins)
	# Get the axes for the lower, upper, and diagonal
	if symmetric:
		indexes = np.arange(nvar * nvar).reshape(nvar, nvar)
		lower = np.tril_indices(nvar, k=-1)
		upper = np.triu_indices(nvar, k=1)
		diag = np.diag_indices(nvar)
		ax_lower = [ax_all[iax] for iax in indexes[lower].tolist()]
		ax_upper = [ax_all[iax] for iax in indexes[upper].tolist()]
		ax_diag = [ax_all[iax] for iax in indexes[diag].tolist()]
	cax_lower = cax_upper = cax_diag = []
	if cbar == 'each':
		if symmetric:
			args = tuple([grid.cbar_axes])
			cax_lower = [grid.cbar_axes[iax] for iax in indexes[lower].tolist()]
			cax_upper = [grid.cbar_axes[iax] for iax in indexes[upper].tolist()]
			cax_diag = [grid.cbar_axes[iax] for iax in indexes[diag].tolist()]
		else:
			args = tuple([grid.cbar_axes])
	else:
		args = tuple()
	# ------------------------------------------------
	# Plot data using operator functions
	# ------------------------------------------------
	if symmetric:
		if isinstance(lowerfunc, types.FunctionType):
			for ax in lowerfunc(ax_lower, *args):
				continue
		else:
			for ax in ax_lower:
				fig.delaxes(ax)
			for cax in cax_lower:
				fig.delaxes(cax)
		if isinstance(diagfunc, types.FunctionType):
			for ax in diagfunc(ax_diag, *args):
				continue
		else:
			for ax in ax_diag:
				fig.delaxes(ax)
			for cax in cax_diag:
				fig.delaxes(cax)
		if isinstance(upperfunc, types.FunctionType):
			for ax in upperfunc(ax_upper, *args):
				continue
		else:
			for ax in ax_upper:
				fig.delaxes(ax)
			for cax in cax_upper:
				fig.delaxes(cax)
	else:
		nplot = 0
		for ax in gridfunc(ax_all, *args):
			nplot += 1
	# ------------------------------------------------
	# Fix some figure aesthetics after plotting
	# ------------------------------------------------

	if aspect is None:
		ax = ax_all[0]
		pltxlim = ax.get_xlim()[1]
		pltylim = ax.get_ylim()[1]
		# Sort out a 4:3 ratio for variograms
		if pltylim == 1.2:
			aspect = pltxlim / 1.2 / (4 / 3)
		# or the 4:3 ratio for CDF's
		elif pltylim == 1.0:
			aspect = pltxlim / 1.0 / (4 / 3)
	for ax in ax_all:
		ax.tick_params(axis='both', pad=2)
		if xlim:
			ax.set_xlim(xlim)
		if ylim:
			ax.set_ylim(ylim)
		if aspect:
			ax.set_aspect(aspect)
		if tight:
			for spine in ax.spines:
				ax.spines[spine].set_visible(True)
		# Set-up the tick labels
		if symmetric and (labelmode.lower() == 'l' or tight):
			ax.set_xlabel('')
			ax.set_ylabel('')
			# y-axis labels
			if isinstance(lowerfunc, types.FunctionType):
				if tight:
					# Handle left labels
					if isinstance(diagfunc, types.FunctionType):
						for iax in indexes[:, 0][::2]:
							ax = ax_all[iax]
							ax.axis['left'].toggle(ticklabels=True, label=True)
						for iax in indexes[:, 0][1::2]:
							ax = ax_all[iax]
							ax.axis['left'].toggle(ticklabels=False, label=False)
					else:
						for iax in indexes[1:, 0][::2]:
							ax = ax_all[iax]
							ax.axis['left'].toggle(ticklabels=True, label=True)
						for iax in indexes[1:, 0][1::2]:
							ax = ax_all[iax]
							ax.axis['left'].toggle(ticklabels=False, label=False)
					# Handle right labels
					if isinstance(diagfunc, types.FunctionType):
						for ax in ax_diag[1::2]:
							ax.axis['right'].toggle(ticklabels=True, label=True)
					elif isinstance(upperfunc, types.FunctionType):
						for iax in indexes[1:, nvar - 1][1::2]:
							ax = ax_all[iax]
							ax.axis['right'].toggle(ticklabels=True, label=True)
					else:
						for iax in indexes[diag][:-1][::2] - 1:
							ax = ax_all[iax]
							ax.axis['right'].toggle(ticklabels=True, label=True)
				else:
					for iax in indexes[:, 0]:
						ax = ax_all[iax]
						ax.axis['left'].toggle(ticklabels=True, label=True)
			elif isinstance(upperfunc, types.FunctionType):
				if tight:
					# Handle right labels
					for iax in indexes[:, -1][::2]:
						ax = ax_all[iax]
						ax.axis['right'].toggle(ticklabels=True, label=True)
					# Handle left labels
					if isinstance(diagfunc, types.FunctionType):
						for ax in ax_diag[1::2]:
							ax.axis['left'].toggle(ticklabels=True, label=True)
					else:
						for iax in indexes[diag][:-1][1::2] + 1:
							ax = ax_all[iax]
							ax.axis['left'].toggle(ticklabels=True, label=True)
				else:
					if isinstance(diagfunc, types.FunctionType):
						for ax in ax_diag:
							ax.axis['left'].toggle(ticklabels=True, label=True)
					else:
						for iax in indexes[diag][:-1] + 1:
							ax = ax_all[iax]
							ax.axis['left'].toggle(ticklabels=True, label=True)
			else:
				for ax in ax_diag:
					ax.axis['left'].toggle(ticklabels=True, label=True)
			# x-axis labels
			if isinstance(lowerfunc, types.FunctionType):
				if tight:
					# Handle bottom labels
					for iax in indexes[-1, :][::2]:
						ax = ax_all[iax]
						ax.axis['bottom'].toggle(ticklabels=True, label=True)
					for iax in indexes[-1, :][1::2]:
						ax = ax_all[iax]
						ax.axis['bottom'].toggle(ticklabels=False, label=False)
					# Handle top labels
					if isinstance(diagfunc, types.FunctionType):
						for ax in ax_diag[1::2]:
							ax.axis['top'].toggle(ticklabels=True, label=True)
					elif isinstance(upperfunc, types.FunctionType):
						for iax in indexes[1:, nvar - 1][1::2]:
							ax = ax_all[iax]
							ax.axis['top'].toggle(ticklabels=True, label=True)
					else:
						for iax in indexes[diag][:-1][::2] - 1:
							ax = ax_all[iax]
							ax.axis['top'].toggle(ticklabels=True, label=True)
				else:
					for iax in indexes[-1, :]:
						ax = ax_all[iax]
						ax.axis['bottom'].toggle(ticklabels=True, label=True)
			elif isinstance(upperfunc, types.FunctionType):
				if tight:
					# Handle top labels
					if isinstance(diagfunc, types.FunctionType):
						for iax in indexes[0, :][::2]:
							ax = ax_all[iax]
							ax.axis['top'].toggle(ticklabels=True, label=True)
					else:
						for iax in indexes[0, 1:][::2]:
							ax = ax_all[iax]
							ax.axis['top'].toggle(ticklabels=True, label=True)
					# Handle bottom labels
					if isinstance(diagfunc, types.FunctionType):
						for iax in indexes[diag][:-1][::2] + 1:
							ax = ax_all[iax]
							ax.axis['bottom'].toggle(ticklabels=True, label=True)
					else:
						for iax in indexes[diag][:-1][1::2] + 1:
							ax = ax_all[iax]
							ax.axis['bottom'].toggle(ticklabels=True, label=True)
				else:
					if isinstance(diagfunc, types.FunctionType):
						for ax in ax_diag:
							ax.axis['bottom'].toggle(ticklabels=True, label=True)
					else:
						for iax in indexes[diag][:-1] + 1:
							ax = ax_all[iax]
							ax.axis['bottom'].toggle(ticklabels=True, label=True)
			else:
				for ax in ax_diag:
					ax.axis['bottom'].toggle(ticklabels=True, label=True)
		elif not symmetric and tight:
			if ((ncol * nrow) - nplot) % 2 == 0:
				icol = 0
			else:
				icol = 1
			for ax in grid.axes_column[ncol - 1][1::2]:
				ax.yaxis.set_visible(True)
				ax.yaxis.set_ticks_position('right')
			for ax in grid.axes_column[0][1::2]:
				ax.axis['left'].toggle(ticklabels=False, label=False)
			for ax in grid.axes_row[0][icol::2]:
				ax.xaxis.set_ticks_position('top')
				ax.xaxis.set_visible(True)
			for ax in grid.axes_row[nrow - 1][icol::2]:
				ax.axis['bottom'].toggle(ticklabels=False, label=False)
	plt.draw()
	# ------------------------------------------------
	# Remove unused plots and move the x-axis tick labels are required for non-symmetric plots
	# ------------------------------------------------
	if not symmetric:
		if nplot < (nrow * ncol):
			for i in range(nplot, (nrow * ncol)):
				# Delete the unused plot
				row, col = divmod(i, ncol)
				delaxes = False   # True == default..but perhaps they shouldnt be deleted?
				if delaxes:
					if hasattr(ax_all[row * ncol + col], 'cax'):
						fig.delaxes(ax_all[row * ncol + col].cax)
					fig.delaxes(ax_all[row * ncol + col])
				else:
					if hasattr(ax_all[row * ncol + col], 'cax'):
						fig.axes[row * ncol + col].cax.set_axis_off()
					fig.axes[row * ncol + col].set_axis_off()
				# Plot the x-axis on the plots above the one being deleted
				row = row - 1
				ax = ax_all[row * ncol + col]
				if not unequal_aspects:
					ax.axis["bottom"].toggle(ticklabels=True, label=True)
				if tight:
					for ax in grid.axes_row[nrow - 2][icol::2]:
						ax.axis['bottom'].toggle(ticklabels=False)
	# remove the labels if the unequal_aspects keyword and labelmode='l'
	if labelmode.lower() == 'l' and unequal_aspects:
		for iplot, ax in enumerate(ax_all):
			irow = np.floor(iplot / ncol)
			icol = iplot - irow * ncol
			if icol == 0 and irow < (nrow - 1):
				ax.set_xticklabels([])
				ax.set_xlabel('')
			elif irow == (nrow - 1) and icol > 0:
				ax.set_yticklabels([])
				ax.set_ylabel('')
			elif icol > 0 and irow < (nrow - 1):
				ax.set_xticklabels([])
				ax.set_yticklabels([])
				ax.set_xlabel('')
				ax.set_ylabel('')
	# ------------------------------------------------
	# Add super axis labels if not using the symmetric variable matrix or variable labels if so
	# ------------------------------------------------
	if symmetric and axislabels:
		# Add the axis labels so they are inline...for some reason the bbox of the ylabels are not
		#  providing the right coordinates so a different method is needed. This will not work if
		#  the bar is put on the left side.
		if tight:
			# Add xaxis labels
			ax = grid.axes_llc
			ax.set_xlabel(axislabels[0])
			ax.set_ylabel(axislabels[nvar - 1])
			renderer = fig.canvas.get_renderer()
			plt.draw()
			label = ax.xaxis.get_label()
			bbox = label.get_window_extent(renderer)
			ymax = bbox.ymax
			ymax = ax.transAxes.inverted().transform((0, ymax))[1]
			for i, ax in enumerate(grid.axes_row[nvar - 1]):
				ax.text(0.5, ymax, axislabels[i], ha='center', va='top', transform=ax.transAxes)
			# Add yaxis labels
			xmin, xmid, ymin, ymid, ymax = get_supaxislocs(fig, nvar, nvar, figsize, 0.01)
			for i, ax in enumerate(grid.axes_column[0]):
				bbox = ax.get_position()
				ax_ymax = bbox.ymax
				ax_ymin = bbox.ymin
				ax_ymid = ((ax_ymax - ax_ymin) / 2) + ax_ymin
				fig.text(xmin, ax_ymid, axislabels[i], rotation=90, ha='right', va='center')
		else:
			for i, ax in enumerate(grid.axes_column[0]):
				ax.set_ylabel(axislabels[i])
			for i, ax in enumerate(grid.axes_row[nvar - 1]):
				ax.set_xlabel(axislabels[i])
	elif (ylabel is not None) or (xlabel is not None) or (suptitle is not None):
		# Add superaxis titles
		for ax in ax_all:
			if xlabel is not None:
				ax.set_xlabel('')
			if ylabel is not None:
				ax.set_ylabel('')
		xmin, xmid, ymin, ymid, ymax = get_supaxislocs(fig, nrow, ncol, figsize, 0.01)
		if ylabel is not None:
			fig.text(xmin, ymid, ylabel, ha='right', va='center', rotation=90)
		if xlabel is not None:
			fig.text(xmid, ymin, xlabel, va='top', ha='center')
		if suptitle is not None:
			fig.text(xmid, ymax, suptitle, va='bottom', ha='center')

	# ------------------------------------------------
	# Mess with the colorbar
	# ------------------------------------------------

	if cbar == 'single':
		# Check to see if the limits of each plot are the same. May only work for pixelplt
		try:
			vlimcheck = False
			for ax in ax_all:
				if not vlimcheck:
					vmin = ax.get_images()[0].norm.vmin
					vmax = ax.get_images()[0].norm.vmax
					vlimcheck = ((vmin, vmax) != tuple(vlim)) or vlimcheck
		except:
			from ..scriptnotifier.utils import printerr as printerr
			printerr(("Please ensure that the data value limits of all subplots is set to the"
					  " same limits passed to image_grid"), errtype='warning')
		if vlimcheck:
			raise ValueError("The single colorbar is being set to a limit of `%s`; however, a plot"
							 " was found to have a limit of `%s`" % (vlim, (vmin, vmax)))
		cax = grid.cbar_axes[0]
		# Grab the image from the first subplot to pass to the colorbar function
		plot = ax_all[0].get_images()[0]
		vlim, ticklocs, ticklabels = get_contcbarargs(np.array([vlim[0], vlim[1]]),
														 sigfigs=3, vlim=vlim)
		# Plot the colorbar
		cbar = fig.colorbar(plot, cax=cax, ticks=ticklocs)
		# Configure the color bar
		cbar.ax.set_yticklabels(ticklabels, ha='left')
		cbar.ax.tick_params(axis='y', pad=2)
		if cbar_label is not None:
			cbar.set_label(cbar_label, ha='center', va='top', labelpad=2)
	# Export figure
	if output_file or ('pdfpages' in kwargs):
		export_image(output_file, **kwargs)

	return fig
