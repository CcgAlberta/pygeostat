#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''contour_plot.py: Contains a basic contour plotting routine using matplotlib'''
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np

from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters



@set_plot_style
def contour_plot(data, griddef=None, var=None, orient='xy', slice_number=0, ax=None, output_file=None, c='k',
               figsize=None, xlabel=None, ylabel=None, title=None, unit=None, leg_label=None,
               aspect=None, clabel=False, lw=1.0, plot_style=None, custom_style=None, axis_xy=None, grid=None, return_ax=True, return_csi=False):
    """
    Contains a basic contour plotting routine using matplotlib

    Parameters:
        data: A numpy ndarray, pandas DataFrame or pygeostat DataFile, where each column is a
            variable and each row is an observation
        griddef (GridDef): A pygeostat GridDef class, which must be provided if a DataFile is
            not passed as data with a valid internal
            GridDef :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`
        var (str,int): The name of the column within data to plot. If an int is provided, then it
            corresponds with the column number in data. If None, the first column of data
            is used.
        orient (str): Orientation to slice data. ``'xy'``, ``'xz'``, ``'yz'`` are the only accepted
            values
        slice_number (int): Grid cell location along the axis not plotted to take the slice of data to
            plot
        ax (mpl.axis): Matplotlib axis to plot the figure
        output_file (str): Output figure file name and location
        show (bool): ``True`` will use plt.show() at end. Typically don't need this.
        c (str): Matplotlib color
        figsize (tuple): Figure size (width, height)
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        title (str): title for the plot
        unit (str): Distance unit, taken from Parameters if ``None``
        leg_label (str): Adds a single label to the legend for the contour lines
        aspect (str): Set a permissible aspect ratio of the image to pass to matplotlib.
        clabel (bool): Whether or not to label the contours wth their values
        lw (float): the weight of the contour lines
        plot_style (str): Optional pygeostat plotting style
        custom_style (dict): Custom dictionary for plotting styles
        grid (bool): Plots the major grid lines if True. Based on Parameters['plotting.grid']
            if None.
        axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
            visible) if axis_xy is True. Based on Parameters['plotting.axis_xy'] if None.
        return_ax (bool): specify if the plotting axis should be returned
        return_csi (bool): specify if the contour instance should be returned


    Returns:
        csi (ax): Matplotlib ax.contour instance

    **Examples:**

    A basic contour plotting example:

    .. plot::

        import pygeostat as gs
        grid_str = '''120 5.00  10.00 -nx, xmn, xsiz
                    110 1205.00  10.00 -ny, ymn, ysiz
                    1 0.5 1.0  -nz, zmn, zsiz'''
        
        griddef = gs.GridDef(grid_str=grid_str)
        data_fl = gs.ExampleData("grid2d_surf", griddef=griddef)
        gs.contour_plot(data_fl, var="Thickness", clabel=True)

    Contour plot on top of slice plot:

    .. plot::

        import pygeostat as gs
        grid_str = '''120 5.00  10.00 -nx, xmn, xsiz
                    110 1205.00  10.00 -ny, ymn, ysiz
                    1 0.5 1.0  -nz, zmn, zsiz'''
        
        griddef = gs.GridDef(grid_str=grid_str)

        data_fl = gs.ExampleData("grid2d_surf", griddef=griddef)
        ax = gs.slice_plot(data_fl, var="Thickness")
        _ = gs.contour_plot(data_fl, var="Thickness", ax = ax, clabel=True)

    """
    from . utils import format_plot, _spatial_griddata, _spatial_labels, _spatial_aspect
    from .export_image import export_image
    from ..datautils.utils import slice_grid

    # Parse the data, var and griddef input to determine the data and griddef
    data, griddef = _spatial_griddata(data, var, griddef)
    # Slice the data
    if orient in ['xy', 'xz', 'yz']:
        view = slice_grid(data, griddef, orient, slice_number)
    else:
        raise Exception("Error: no orientation set! {}".format(orient))

    a = orient[0]
    b = orient[1]
    xmin = getattr(griddef, a + 'limits')[0]
    xmax = getattr(griddef, a + 'limits')[1]
    ymin = getattr(griddef, b + 'limits')[0]
    ymax = getattr(griddef, b + 'limits')[1]

    if aspect is None:
        aspect = _spatial_aspect([xmin, xmax], [ymin, ymax])

    # Discretize for contour map
    X = np.linspace(xmin, xmax, len(view[0, :]))
    Y = np.linspace(ymin, ymax, len(view[:, 0]))

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=figsize)

    # Create the 'artists' that draw the contours
    csi = ax.contour(X, Y, np.reshape(view, (len(view[:, 0]), len(view[0, :]))), linewidths=lw,
                     colors=c)
    if leg_label is not None:
        ax.plot(np.nan, lw=lw, color=c, label=leg_label)
    if clabel:
        ax.clabel(csi, inline=1, fontsize=9, fmt='%.0f')
    if aspect:
        ax.set_aspect(aspect)
    _spatial_labels(ax, orient, griddef, slice_number, title, xlabel, ylabel, unit)
    if axis_xy is None:
        axis_xy = Parameters['plotting.axis_xy_spatial']
    format_plot(ax, axis_xy=axis_xy, grid=grid)

    # Done Plotting, save figure if required and show it
    if output_file is not None:
        export_image(output_file)

    # return some things
    if return_ax and return_csi:
        return ax, csi
    elif return_ax:
        return ax
    elif return_csi:
        return csi
