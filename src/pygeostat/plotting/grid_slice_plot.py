#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Function that emulates the GSLIB program plotem when used with slice_plot"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters


@set_plot_style
def grid_slice_plot(data, griddef = None, super_xlabel=True, super_ylabel=True, super_title=None, ncol=None, nrow=None,
                    start_slice=None, end_slice=None, figsize=None, n_slice=None, slice_title=True, unit=None,
                    plot_style=None, custom_style=None, output_file=None, out_kws=None, cbar_label=None,
                    axpad=0.15, cbar_cats=None, ntickbins=None, axfuncs=None, label_mode='L', **kwargs):
    """
    grid_slice_plot can be used to automatically generate set of slices
    through a 3D gridded realization. Given some target number of rows, columns, orientation and
    slice ranges through the realization, this function will automatically generate the slice_plot
    slices and arrange them according to the specified dimensions. It is possible to pass keyword
    arguments for slice_plot to this function in order to specify the format of the slice_plots. So,
    adding data locations, different colormaps, and other slice_plot formatting is permitted for
    all subplots by passing those slice_plot arguments to this function. See
    :func:`gs.slice_plot()<pygeostat.plotting.slice_plot>` for the full list of permissable kwargs.

    Updates April 2016 - use a ImageGrid subplots to get the job done

    Parameters:
        data (array, dataframe) : array of data, passed directly to slice_plot()
        griddef (pygeostat griddef) : pygeostat grid definitions, passed directly to slice_plot()
        super_xlabel (str) :  super x axis label
        super_ylabel (str) :  super y axis label
        super_title (str) : super title for the subplots
        ncol (int) : the number of columns considered for the subplots (may change)
        nrow (int) : the number of rows considered for the subplots (may change)
        start_slice (int) : the starting slice to be plotted
        end_slice (int) : the end slice to be plotted
        figsize (tuple) : size of the figure to be created
        n_slice (int) : the number of desired slices
        slice_title (bool) : either plot the orientation and slice no on top of each slice, or dont!
        unit (str): Unit to place inside the axis label parentheses
        plot_style (str) : Use a predefined set of matplotlib plotting parameters as specified by
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
            to turn it off
        custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.exportimg() <pygeostat.plotting.exportimg.exportimg>`
        cbar_label (str): colorbar title
        axpad (float): figure points padding applied to the axes, vertical padding is further
            modified to account for the slice titles, if required
        ntickbins (int or tuple): The number of tick bins for both axes, or the (x, y) respectively
        axfuncs (function or list of functions): External function(s) that takes `ax` `slice_number` and
            `orient` as keyword arguments, does not return anything
        label_mode (str): default `'L'`, or `'all'`, passed to the ImageGrid() constructor
        **kwargs  : **NOTE** the arguments here are either valid slice_plot (including all keyword)
            dictionary arguments, or valid imshow and valid imshow keyword arguments. If errors are
            thrown from invalid arguments it is likely that something that shouldnt have been
            passed to imshow was passed. Check and double check those \**kwargs!

    Returns:
        fig (plt.figure) : figure handle


    **Examples:**

    A simple call: generates a set of slices through the model

    .. plot::

        import pygeostat as gs
        # load some data
        data = gs.ExampleData('3d_grid',griddef = gs.GridDef([40,1,2,40,1,2,40,0.5,1]))
        # plot the grid slices
        _ = gs.grid_slice_plot(data)

    |

    Possible to specify the orientation and the number of slices:

    .. plot::

        import pygeostat as gs
        # load some data
        data = gs.ExampleData('3d_grid',griddef = gs.GridDef([40,1,2,40,1,2,40,0.5,1]))
        # plot the grid slices
        _ = gs.grid_slice_plot(data, orient='xz', n_slice=5)
        

    |

    Can specify the number of rows or columns required for the slices:

    .. plot::

        import pygeostat as gs
        # load some data
        data = gs.ExampleData('3d_grid',griddef = gs.GridDef([40,1,2,40,1,2,40,0.5,1]))
        # plot the grid slices
        _ = gs.grid_slice_plot(data, orient='xz', n_slice=6, ncol=2, nrow=3)

    |

    Also able to specify slice_plot kwargs using this function, so we can apply consistent
        custom formatting to all of the subplots:

    .. plot::
    
        import pygeostat as gs
        # load some data
        data = gs.ExampleData('3d_grid',griddef = gs.GridDef([40,1,2,40,1,2,40,0.5,1]))
        # plot the grid slices
        _ = gs.grid_slice_plot(data, nrow=2, ncol=5, start_slice=10, end_slice=25, n_slice=10, cmap='hot', vlim=(-3,3))

    """
    from . slice_plot import slice_plot
    from .utils import (addxticks, get_contcbarargs, get_supaxislocs)
    from ..datautils.utils import slice_grid
    from ..utility.logging import printerr
    from .export_image import export_image

    from ..data.data import DataFile

    if axfuncs is not None and not isinstance(axfuncs, list):
        axfuncs = [axfuncs]
    if out_kws is None:
        out_kws = dict()

    if griddef is None:
        if isinstance(data, DataFile):
            griddef = data.griddef
        else:
            raise ValueError ('griddef must be provided')
    
    if isinstance(data, DataFile):
        data = data.data

    if 'orient' in kwargs:
        ori = kwargs.get('orient')
        del kwargs['orient']
    else:
        ori = 'xy'
    if slice_title and 'title' in kwargs:
        del kwargs['title']
    if 'slice_number' in kwargs:
        del kwargs['slice_number']
    if 'ax' in kwargs:
        del kwargs['ax']

    if 'catdata' in kwargs:
        catdata = kwargs.get('catdata')
    else:
        if len(np.unique(data)) <= 12:
            catdata = True
        else:
            catdata = False

    if start_slice is None or start_slice < 0:
        start_slice = 0

    if hasattr(data, 'values'):
        data = np.array(data)

    di = ori[1]
    if end_slice is None or end_slice > getattr(griddef, 'n' + di):
        if ori == 'xy':
            end_slice = griddef.nz
        elif ori == 'xz':
            end_slice = griddef.ny
        elif ori == 'yz':
            end_slice = griddef.nx

    if end_slice <= 0:
        print('Bad end_slice parameter, defaulting to 1')
        end_slice = 1

    if n_slice is None:
        n_slice = 4

    if nrow is None:
        addrows = True
        nrow = 2
    else:
        addrows = False
    if ncol is None:
        addcols = True
        ncol = 2
    else:
        addcols = False

    # Try to automatically expand the rows,columns if we dont have enough for the number of slices
    alternate = 1
    if not addrows and not addcols and nrow * ncol < n_slice:

        printerr('Please ensure there are enough rows and columns explicitly defined for '
                 'the target number of plots! (n_slice = %i, nrow*ncol = %i)' %
                 (n_slice, nrow * ncol))
        raise Exception()
    while ncol * nrow < n_slice:
        if alternate == 1:
            if addcols:
                ncol += 1
            alternate = 0
        else:
            if addrows:
                nrow += 1
            alternate = 1

    # Deal with some of the kwarg inputs (add any other defaults that are needed?)
    if 'sigfigs' in kwargs:
        sigfigs = kwargs.get('sigfigs')
    else:
        sigfigs = 3

    # set the vlim to the min max of the entire dataset if none is provided
    if 'vlim' in kwargs:
        vlim = kwargs.get('vlim')
        del kwargs['vlim']
    else:
        vlim = (np.nanmin(data), np.nanmax(data))

    # Turn off slice_plot axis labels if we want super axis labels:
    if super_xlabel:
        kwargs['xlabel'] = False
    if super_ylabel:
        kwargs['ylabel'] = False

    # Create the subplots on the new ImageGrid
    if figsize is None:
        figsize = mpl.rcParams.get("figure.figsize")  # but also needed for later
    fig = plt.figure(figsize=figsize)  # this should just pull from rcParams?
    if slice_title:
        if isinstance(mpl.rcParams['axes.titlesize'], str):
            scalar = 1.1
        else:  # assume this operation will work..
            scalar = mpl.rcParams['axes.titlesize'] / 8
        vaxpad = (axpad * scalar) + 0.1
    else:
        vaxpad = axpad
    imgrid = ImageGrid(fig, 111, nrows_ncols=(nrow, ncol), axes_pad=(axpad, vaxpad),
                       label_mode=label_mode, direction='row', cbar_mode='single', cbar_size=0.15)

    axind = 0
    # ensure a suitable number of slices
    steps = max(int(np.floor((end_slice - start_slice) / n_slice)), 1)
    slrange = np.arange(start_slice, end_slice, steps)[0:n_slice]

    num_slices = 0
    num_errors = 0
    for sl in slrange:
        num_slices += 1
        if slice_title:
            ptitle = ori + ' slice at: ' + str(round(griddef.get_slice_coordinate(ori, sl), 1))
        else:
            ptitle = ""
        if axind < nrow * ncol:
            try:
                ax, plot = slice_plot(data, griddef, ax=imgrid[axind], slice_number=sl, orient=ori,
                                      vlim=vlim, title=ptitle, cbar=False,return_plot=True, plot_style=False,
                                       output_file=None, out_kws=None, **kwargs)
                if ntickbins is not None:
                    if isinstance(ntickbins, tuple):
                        ax.locator_params(axis='x', nbins=ntickbins[0])
                        ax.locator_params(axis='y', nbins=ntickbins[1])
                    else:
                        ax.locator_params(nbins=ntickbins)
                if axfuncs:
                    for func in axfuncs:
                        func(slice_number=sl, ax=imgrid[axind], orient=ori)
                axind += 1
            except ValueError as error:  # hopefully this error was raised by an empty slice..
                num_errors += 1
                err = error
                pass

    if num_errors == num_slices:
        raise ValueError(err)

    # remove unused axes from the figure
    [x.set_axis_off() for x in imgrid[axind:]]

    # Add xticklabels to the last subplot in a column when it is not in the last row
    addxticks(fig, nrow, ncol, axind, griddef.extents(orient=ori[0]))

    # the vlim is either given or set to the max range for the entire array of data..
    # could get more clever about this by preanalyzing the data slices and determining the
    # correct range from all the slices...
    if catdata:  # this section is straight from the slice_plot section for categorical variables
        ticklabels = np.unique(data[~np.isnan(data)])
        if isinstance(cbar_cats, dict):
            ticklabels = [cbar_cats[cat] for cat in ticklabels]
        else:
            ticklabels = [int(x) for x in ticklabels]
        ncat = len(ticklabels)
        vlim = (0, ncat)
        ticklocs = np.arange(ncat) + 0.5
    else:
        vlims, ticklocs, ticklabels = get_contcbarargs(slice_grid(data, griddef, ori, sl),
                                                          sigfigs, vlim)
    # get the min and max colorbar arguments
    cbar = fig.colorbar(plot, imgrid.cbar_axes[0], ticks=ticklocs)
    cbar.ax.set_yticklabels(ticklabels, ha='left')
    cbar.ax.tick_params(axis='y', pad=2)
    if cbar_label is not None:
        cbar.set_label(cbar_label, ha='center', va='top', labelpad=2)

    # Get the figure coordates for the superaxis titles and the figure title
    if super_xlabel or super_ylabel or super_title:
        xmin, xmid, ymin, ymid, ymax = get_supaxislocs(fig, nrow, ncol, figsize, 0.01)

    if unit is None:
        unit = Parameters['plotting.unit']
    if unit is None or unit == '':
        unit = ''
    else:
        unit = ' ({})'.format(unit)
    if super_xlabel:
        if not isinstance(super_xlabel, str):
            if ori == 'xy':
                super_xlabel = Parameters['plotting.xname'] + unit
            elif ori == 'xz':
                super_xlabel = Parameters['plotting.xname'] + unit
            elif ori == 'yz':
                super_xlabel = Parameters['plotting.yname'] + unit
        plt.draw()  
        fig.text(xmid, ymin, super_xlabel, va='top', ha='center')
    if super_ylabel:
        if not isinstance(super_ylabel, str):
            if ori == 'xy':
                super_ylabel = Parameters['plotting.yname'] + unit
            elif ori == 'xz':
                super_ylabel = Parameters['plotting.zname'] + unit
            elif ori == 'yz':
                super_ylabel = Parameters['plotting.zname'] + unit
        plt.draw()
        fig.text(xmin, ymid, super_ylabel, ha='right', rotation=90, va='center')
    if super_title:
        plt.draw()
        fig.text(xmid, ymax, super_title, va='bottom', ha='center')

    # Export figure
    if output_file or ('pdfpages' in out_kws):
        export_image(output_file=output_file, **out_kws)

    return fig
