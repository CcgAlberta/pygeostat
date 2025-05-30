#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Pit plotting routine with 2D slices using matplotlib"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters


@set_plot_style
def pit_plot(arr, griddef, ax=None, orient='xz', slice_number=0, lineweight=1, color='k', iso=0.5,
            linestyle='solid', figsize=None, xlim=None, ylim=None, title=None, xlabel=None,
            ylabel=None, unit=None, rotateticks=None, randomize=False, grid=None, axis_xy=None,
            cust_style=None, label=None, output_file=None, out_kws=None):
    """This funcion will take an array of indicators (from lg3d) and an orientation, and plot the
    pit shell lines for a given cross section view.

    Note:
        This function can only deal with 1 realization in the file
        so if you have multiple realizations you need to either pass a slice to this function or copy 1
        realization to a separate file.

    Parameters:
        arr (Array): Array (DataFrame Column) passed to the program with indicator
            values (i.e. 0 & 1)
        ax (mpl.axis): Matplotlib axis to plot the figure
        orient (str): Orientation to slice data. 'xy', 'xz', 'yz' are the only accepted values
        slice_number (int): Location of slice to plot
        lineweight (float): Any Matplotlib line weight
        color (str): Any Matplotlib color
        iso (float): Inside or Outside of Pit limit (i.e. if greater than 0.5 inside of pit)
        linestyle (str): Any Matplotlib linestyle
        randomize (bool): True or False... obviously
        figsize (tuple): Figure size (width, height)
        title (str): Title for the plot
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        unit (str): Unit to place inside the axis label parentheses
        rotateticks (bool tuple): Indicate if the axis tick labels should be rotated (x, y)
        grid (bool): Plots the major grid lines if True. Based on Parameters['plotting.grid']
            if None.
        axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
            visible) if axis_xy is True. Based on Parameters['plotting.axis_xy'] if None.
        label (str): Legend label to be added to Matplotlib axis artist
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`


    Returns:
        fig (fig): Matplotlib figure instance

    Examples:
        A simple call:

        >>> gs.pit_plot(data.data, data.griddef, title='Pit Outline Using LG3D output')

        .. image:: ./figures/PlottingGallery/pitplot.png

        In order to plot multiple pits (say from a file with multiple realizations) you have can
        plot to the same matplotlib axis. For multiple realizations using a loop is the easiest
        as shown below.

        >>> sim = gs.DataFile(flname='SGS.lg3d', griddef=grid_5m)

        Loop through the SGSIM LG3D output file

        First plot the first realization and grab the matplotlib axis

        >>> import matplotlib.pyplt as plt
        ... rmin = 0
        ... rmax = pit.griddef.count()
        ... fig = gs.pit_plot(sim.data[rmin:rmax], sim.griddef, title='Pit Outline Using LG3D output
        ...                 with Multiple Realizations')
        ... ax = fig.gca()

        Then loop through the rest of the realizations (Say 50) and plot them on current axis

        >>> for i in range (1, 50):
        ...     rmin = i*pit.griddef.count()
        ...     rmax = rmin + pit.griddef.count()
        ...     gs.pit_plot(sim.data[rmin:rmax], sim.griddef, ax=ax)

        Save the figure

        >>> gs.export_image('pitplot_mr.png', format='png')

        .. image:: ./figures/PlottingGallery/pitplot_mr.png
    """
    from .export_image import export_image
    from . utils import format_plot, _spatial_labels, _tickoverlap
    
    # Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()
    # Set-up plot if no axis is supplied
    if ax is None:
        fig, ax = plt.subplots(1, figsize=figsize)
    # Handle panda input
    if hasattr(arr, 'values'):
        arr = arr.values

    if arr.shape != (griddef.nz, griddef.ny, griddef.nx):
        arr = arr.reshape((griddef.nz, griddef.ny, griddef.nx))

    binarr = np.zeros_like(arr, dtype=bool)
    binarr[arr > iso] = True

    if orient == 'xy':
        # view = binarr[slice_number, :, :]
        print('Orientation of ', orient, ' not yet supported - to be added!')
        return
    elif orient == 'xz':
        view = binarr[:, slice_number, :]
        # Get limits for orient
    elif orient == 'yz':
        view = binarr[:, :, slice_number]
        # Get limits for orient
    else:
        print('Orientation of ', orient, ' not supported')
        return

    a = orient[0]
    b = orient[1]

    xmin = getattr(griddef, a + 'limits')[0]
    xmax = getattr(griddef, a + 'limits')[1]
    ymin = getattr(griddef, b + 'limits')[0]
    ymax = getattr(griddef, b + 'limits')[1]

    x = []
    y = []

    xmn = getattr(griddef, a + 'mn')
    nx = getattr(griddef, 'n' + a)
    xsiz = getattr(griddef, a + 'siz')

    for i in range(1, nx - 1):
        for v, k in zip(view[:, i], range(griddef.nz)):
            if v:
                x = x + [xmn + i * xsiz]
                y = y + [griddef.zmn + k * griddef.zsiz - griddef.zsiz / 2]
                break
            if k == griddef.nz - 1:
                if view[k, i + 1] or view[k, i - 1]:
                    x = x + [xmn + i * xsiz]
                    y = y + [griddef.zmn + k * griddef.zsiz + griddef.zsiz / 2]

    if randomize:
        r = (np.random.random() - 0.5) * xsiz
        for i in range(len(x)):
            x[i] = x[i] + r

    if len(x) > 2:
        verts = list(zip(x, y))
        codes = [Path.MOVETO] + [Path.LINETO] * (len(x) - 1)

        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', ec=color, lw=lineweight, ls=linestyle)

        ax.add_patch(patch)
        fig = plt.gcf()
    # Add Legend artist to axis if label is passed
    if label:
        line1, = ax.plot([-9999999, -9999999, -9999999], label=label, linestyle=linestyle,
                         linewidth=lineweight, color=color)
        line1.set_label(label)
    # Set axis limits
    if xlim is None:
        xlim = (xmin, xmax)
    if ylim is None:
        ylim = (ymin, ymax)
    ax = _spatial_labels(ax, orient, griddef, slice_number, title, xlabel, ylabel, unit)
    if axis_xy is None:
        axis_xy = Parameters['plotting.axis_xy_spatial']
    format_plot(ax, xlim=xlim, ylim=ylim, grid=grid, axis_xy=axis_xy)
    # Note on Tick Labels:
    #   If a group of subplots are put together which share a x-axis, the rotation may not work. By
    #   getting the tick labels generated by matplotlib as a set of label objects, they can be
    #   looped through and have their settings individually fixed. This appears to be the only way
    #   to have the shared axis labels formated properly. The labels are also adjusted closer to the
    #   axis for esthetics. --Warren E. Black

    # The plots tick labels will not be properly accessible until the figure is "drawn", once the
    # command below is run, ax.get_ticklabel() will actually work properly.
    plt.draw()
    # Check to see if the ytick labels need to be rotated if the rotate argument was not passed.
    # This doesn't work if the figure is being placed into a subplot (i.e., is not a standalone
    # figure). Don't know why...
    if rotateticks is None:
        rotateticks = Parameters['plotting.rotateticks']
    if rotateticks is None:
        rotateticks = _tickoverlap(ax)
    # Configure y-axis tick labels
    ylabels = ax.get_yticklabels()
    for ylabel in ylabels:
        ylabel.set_ha('right')
        ylabel.set_va('center')
        # Rotate if required
        if rotateticks[1]:
            ylabel.set_rotation(90)
    # Configure x-axis tick labels
    xlabels = ax.get_xticklabels()
    for xlabel in xlabels:
        xlabel.set_ha('center')
        xlabel.set_va('top')
        # Rotate if required
        if rotateticks[0]:
            xlabel.set_rotation(45)
    # Fix tick label padding
    ax.tick_params(axis='both', pad=2)
    if rotateticks[1]:
        ax.tick_params(axis='y', pad=0)
    # Export figure
    if output_file or ('pdfpages' in out_kws):
        export_image(output_file, **out_kws)
    if ax is None:
        return fig
    else:
        return ax
