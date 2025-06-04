#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Generate a correlation matrix plot"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


from scipy import cluster
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters


@set_plot_style
def correlation_matrix_plot(correlation_data, figsize=None, ax=None, cax=None, title=None, xticklabels=None,
                            ticklabels=None, yticklabels=None, rotateticks=None, cbar=None, annotation=None, lower_matrix=False,
                            lw=0.5, hierarchy=None, dendrogram=False, vlim=(-1, 1), cbar_label=None, cmap=None,
                            plot_style=None, custom_style=None, output_file=None, out_kws=None, sigfigs=3, **kwargs):
    """
    This function uses matplotlib to create a correlation matrix heatmap illustrating the
    correlation coefficient between each pair of variables.

    The only parameter needed is the correlation matrix. All of the other arguments are optional.
    Figure size will likely have to be manually adjusted. If the label parameters are left to their
    default value of ``None`` and the input matrix is contained in a pandas dataframe, the
    index/column information will be used to label the columns and rows. If a numpy array is
    passed, axis tick labels will need to be provided. Axis tick labels are automatically checked
    for overlap and if needed, are rotated. If rotation is necessary, consider condensing the
    variables names or plotting a larger figure as the result is odd. If ``cbar`` is left to its
    default value of ``None``, a colorbar will only be plotted if the ``lower_matrix`` is set to True. It
    can also be turned on or off manually. If ``annotation`` is left to its default value of ``None``,
    annotations will only be placed if a full matrix is being plotted. It can also be turned on or
    off manually.

    The parameter ``ticklabels`` is odd in that it can take a few forms, all of which are a tuple
    with the first value controlling the x-axis and second value controlling the y-axis (x, y). If
    left to its default of ``None``, another pygeostat function will check to see if the labels
    overlap, if so it will rotate the axis labels by a default angle of (45, -45) if required. If
    a value of ``True`` is pased for either axis, the respective default values previously stated
    is used. If either value is a float, that value is used to rotate the axis labels.

    The correlation matrix can be ordered based on hierarchical clustering. The following is a list
    of permissible arguments: ``'single', 'complete', 'average', 'weighted', 'centroid', 'median',
    'ward'``. The denrogram if plotted will have a height equal to 15% the height of the
    correlation matrix. This is currently hard coded.

    .. seealso::
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    Please review the documentation of the :func:`gs.set_style()
    <pygeostat.plotting.set_style.set_style>` and :func:`gs.export_image()
    <pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
    their use in this function can be understood.

    Parameters:
        correlation_data: Pandas dataframe or numpy matrix containing the required loadings or correlation
            matrix
        figsize (tuple): Figure size (width, height)
        ax (mpl.axis): Matplotlib axis to plot the figure
        title (str): Title for the plot
        ticklabels (list): Tick labels for both axes
        xticklabels (list): Tick labels along the x-axis (overwritten if ticklabels is passed)
        yticklabels (list): Tick labels along the y-axis (overwritten if ticklabels is passed)
        rotateticks (bool or float tuple): Bool or float values to control axis label rotations.
            See above for more info.
        cbar (bool): Indicate if a colorbar should be plotted or not
        annotation (bool): Indicate if the cells should be annotationated or not
        lower_matrix (bool): Indicate if only the lower matrix should be plotted
        lw (float): Line width of lines in correlation matrix
        hierarchy (str): Indicate the type of hieriarial clustering to use to reorder the
            correlation matrix. Please see above for more details
        dendrogram (bool): Indicate if a dendrogram should be plotted. The argument ``hierarchy``
            must be set to ``true`` for this argument to have any effect
        vlim (tuple): vlim for the data on the correlation_matrix_plot, default = (-1, 1)
        cbar_label (str): string for the colorbar label
        cmap (str): valid Matplotlib colormap
        plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
            to turn it off
        custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
        sigfigs (int): significant digits for labeling of colorbar and cells
        **kwargs: Optional permissible keyword arguments to pass to matplotlib's pcolormesh
            function

    Returns:
        ax (ax): matplotlib Axes object with the correlation matrix plot

    **Examples:**
        Calculate the correlation matrix variables in a pandas dataframe

        .. plot::
        
            import pygeostat as gs
            data_file = gs.ExampleData("point3d_ind_mv")
            data = data_file[data_file.variables]
            data_cor = data.corr()
            gs.correlation_matrix_plot(data_cor, cmap = 'bwr')

        |

        Again for illustration, convert the correlation dataframe into a numpy matrix. By using a
        numpy matrix, the axis labels will need to me manually entered. Reduce the figure size as
        well:

        .. plot::
        
            import pygeostat as gs
            data_file = gs.ExampleData("point3d_ind_mv")
            data = data_file[data_file.variables]
            data_cor = data.corr()
            gs.correlation_matrix_plot(data_cor.values, cmap = 'bwr')

        |

        Plotting a lower correlation matrix while having
        annotations:

        .. plot::

            import pygeostat as gs

            data_file = gs.ExampleData("point3d_ind_mv")
            data = data_file[data_file.variables]
            data_cor = data.corr()
            gs.correlation_matrix_plot(data_cor, lower_matrix=True, annotation=True)
                
    """
  
    from . export_image import export_image
    from . utils import titleoverlap, _tickoverlap, get_contcbarargs
    from mpl_toolkits.axes_grid1 import ImageGrid
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Sanity checks
    if lower_matrix and dendrogram:
        raise NotImplementedError("Dendrogram plotting while using the ``lower_matrix`` functionality is"
                                  " not currently supported")
    # Convert the numpy array to a pd.DataFrame (needed for hierarchy)
    if isinstance(correlation_data, np.ndarray):
        correlation_data = pd.DataFrame(data=correlation_data)

    if not isinstance (correlation_data, pd.DataFrame):
        raise ValueError('correlation_data must be convertable to pandas dataframe')

    # Set-up some parameters
    nx = correlation_data.shape[1]
    ny = correlation_data.shape[0]
    #  Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()
    # Plot Colourbar if required
    if cbar is None:
        if lower_matrix or annotation is False:
            cbar = True
        else:
            cbar = False
    #  Determine hierarchy if needed
    if hierarchy in ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']:
        # Determine the clustering
        linkage = cluster.hierarchy.linkage(correlation_data, method=hierarchy)
        dendo = cluster.hierarchy.dendrogram(linkage, no_plot=True)
        hierarchy = True
        # Reorder the correlation matrix and the tick labels
        ind = dendo['leaves']
        xlab = list(correlation_data.columns)
        xorder = [xlab[i] for i in ind]
        correlation_data = correlation_data[xorder]
        ylab = list(correlation_data.index)
        yorder = [ylab[i] for i in ind]
        correlation_data = correlation_data.loc[yorder]
        if ticklabels:
            ticklabels = [ticklabels[i] for i in ind]
        if xticklabels:
            xticklabels = [xticklabels[i] for i in ind]
        if yticklabels:
            yticklabels = [yticklabels[i] for i in ind]
    # Copy the data and convert it to a numpy matrix if a pandas dataframe was passed
    if isinstance(correlation_data, pd.DataFrame):
        plot_data = correlation_data.values
    else:
        plot_data = np.asarray(correlation_data)
    # Set-up plot else coerce the passed axes as required
    if ax is None:
        if cbar:
            if hierarchy:
                cbar_mode = 'each'
            else:
                cbar_mode = 'single'
        else:
            cbar_mode = None
        if hierarchy and dendrogram:
            nrows_ncols = (2, 1)
        else:
            nrows_ncols = (1, 1)
        # Setup up a new plot
        fig = plt.figure(figsize=figsize)
        imggrid = ImageGrid(fig, 111, nrows_ncols, axes_pad=(0.07, 0), cbar_mode=cbar_mode, cbar_size=0.075)
        ax = imggrid[0]
        if hierarchy and dendrogram:
            ax_dendo = imggrid[1]
            fig.delaxes(imggrid.cbar_axes[1])
        if cbar:
            cax = imggrid.cbar_axes[0]
    elif (cbar or hierarchy):
        fig = plt.gcf()
        divider = make_axes_locatable(ax)
        if cbar:
            cax = divider.append_axes("right", size=0.075, pad=0.07, aspect='auto')
        if hierarchy and dendrogram:
            ax_dendo = divider.append_axes("bottom", size=0.4, pad=0.0, aspect='auto')
    elif hierarchy:
        raise ValueError("`ax` cannotation be divided meaning the dendrogram cannotation be plotted")
    elif cbar and cax is None:
        raise ValueError("A colorbar axes `cax` must be passed as the passed `ax` cannotation be"
                         " divided.")
    ax.set_aspect('equal')
    # Set the axis ticklabels if possible
    if ticklabels:
        xlabels = ticklabels
        ylabels = ticklabels
    else:
        # Handle xticklabels
        if xticklabels is None and hasattr(correlation_data, 'columns'):
            xticklabels = list(correlation_data.columns)
        elif isinstance(xticklabels, bool) and xticklabels:
            xticklabels = list(correlation_data.columns)
        if xticklabels is not None:
            xlabels = xticklabels
        else:
            xlabels = None
        # Handle yticklabels
        if yticklabels is None and hasattr(correlation_data, 'index'):
            yticklabels = list(correlation_data.index)
        elif isinstance(yticklabels, bool) and yticklabels:
            yticklabels = list(correlation_data.index)
        if yticklabels is not None:
            ylabels = yticklabels
        else:
            ylabels = None
    # Set-up figure estetics
    if lower_matrix:
        gridclr = 'white'
        ax.set(xlim=(0, nx - 1), ylim=(0, ny - 1))
        xticklocs = np.arange(nx - 1)
        yticklocs = np.arange(ny - 1)
        # Trim the labels
        if isinstance(xlabels, list):
            xlabels = xlabels[1:]
        if isinstance(ylabels, list):
            ylabels = ylabels[1:]
    else:
        gridclr = 'black'
        ax.xaxis.tick_top()
        ax.set(xlim=(0, nx), ylim=(0, ny))
        xticklocs = np.arange(nx)
        yticklocs = np.arange(ny)
    # Set-up x-axis labels and grid locations
    ax.set_xticks(xticklocs + 0.5)
    ax.set_xticks(xticklocs, minor=True)
    if xlabels is not None:
        if lower_matrix:
            va = 'top'
        else:
            va = 'bottom'
        ax.set_xticklabels(xlabels, va=va, ha='center', rotation='horizontal')
        ax.tick_params(axis='x', pad=2)
    else:
        ax.get_xaxis().set_ticks([])
    # Set-up y-axis labels and grid locations
    ax.invert_yaxis()
    ax.set_yticks(yticklocs + 0.5)
    ax.set_yticks(yticklocs, minor=True)
    if ylabels is not None:
        ax.set_yticklabels(ylabels, va='center', ha='right', rotation='vertical')
        ax.tick_params(axis='y', pad=1)
    else:
        ax.get_yaxis().set_ticks([])
    # Set-up the figure spines
    for spine in ax.spines:
        if lower_matrix:
            ax.spines[spine].set_color('white')
    # Mask the data if lower_matrix is being used
    if lower_matrix:
        mask = np.zeros_like(correlation_data, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True
        plot_data = np.ma.masked_where(mask, plot_data)
        plot_data = plot_data[1:, :(ny - 1)]
    # Plot the grid
    ax.grid(True, which='minor', color=gridclr, zorder=3, lw=lw)
    # Tick rotation
    plt.draw()
    # Check if the axis tick labels overlap, if so rotate them
    if rotateticks is None:
        rotateticks = Parameters['plotting.rotateticks']
    if rotateticks is None:
        # The plots tick labels will not be properly accessible until the figure is "drawn", once
        # the command below is run, ax.get_ticklabel() will actually work properly.
        rotateticks = _tickoverlap(ax)
    # Rotate if required
    if rotateticks[0] is not False or rotateticks[0] is not None:
        if rotateticks[0] is True:
            rotateticks[0] = 45
        xlabels = ax.get_xticklabels()
        for xlabel in xlabels:
            if lower_matrix:
                xlabel.set_ha('center')
                xlabel.set_va('top')
            else:
                xlabel.set_ha('center')
                xlabel.set_va('bottom')
            xlabel.set_rotation(rotateticks[0])
    if rotateticks[1] is not False or rotateticks[1] is not None:
        if rotateticks[1] is True:
            rotateticks[1] = -45
        ylabels = ax.get_yticklabels()
        for ylabel in ylabels:
            ylabel.set_ha('right')
            ylabel.set_va('center')
            ylabel.set_rotation(rotateticks[1])
        ax.tick_params(axis='y', pad=2)
    # Plot the title if required
    if title:
        titletxt = ax.set_title(title)
        # Due to placing the xticklaebls on the top, if a title is plotted it is going to overlap.
        # The following code checks for overlap and bumps the title up step by step until there
        # isn't any
        if not lower_matrix:
            # Check to see if there is overlap with the title and the xticklabels The plots tick
            # labels will not be properly accessible until the figure is "drawn", once the command
            # below is run, ax.titleoverlap() will actually work properly.
            plt.draw()
            _titleoverlap = titleoverlap(ax, titletxt)
            # If there is, clear the title and start moving a suptitle up until there isn't overlap
            if _titleoverlap:
                titletxt.set_text('')
                shifttitle = True
                y = 1.01
            else:
                shifttitle = False
            while _titleoverlap:
                titletxt = ax.set_title(title, y=y)
                plt.draw()
                _titleoverlap = titleoverlap(ax, titletxt)
                y = y + 0.01
            # Now that a spot without overlap has been found, add the pad of 0.015 (0.01 alread
            # added to y) so that it is slightly farther away from the axis labels
            if shifttitle:
                titletxt = ax.set_title(title, y=(y + 0.005))
    # Plot the figure
    if cmap is None:
        cmap = Parameters['plotting.cmap']
    plot = ax.pcolormesh(plot_data, cmap=cmap, norm=plt.Normalize(vmin=vlim[0], vmax=vlim[1]),
                         zorder=0)
    # annotationate if required
    if not isinstance(annotation, bool):
        if lower_matrix:
            annotation = False
        else:
            annotation = True
    if annotation:
        # Set-up a colormap to use to color the annotationations. This was manually tuned so if you can
        # come up with something better do it
        clrvals = [0.85, 0.85, 0.85, 0.85, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        if cmap is None:
            cmap = []
            for clr in clrvals:
                cmap.append(str(clr))
        if isinstance(cmap, list):
            cmap = mpl.colors.ListedColormap(cmap)
            clrnorm = mpl.colors.Normalize(vmin=vlim[0], vmax=vlim[1])
            clrmap = mpl.cm.ScalarMappable(norm=clrnorm, cmap=cmap)
        for y in range(plot_data.shape[0]):
            for x in range(plot_data.shape[1]):
                if isinstance(plot_data[y, x], float):
                    # try:
                    #     color = clrmap.to_rgba(plot_data[y, x])
                    # except (KeyError, ValueError):
                    color = 'black'
                    ax.text(x + 0.5, y + 0.5, ('{:.%ig}' % sigfigs).format(plot_data[y, x]),
                            ha='center', va='center', color=color)
    if cbar:
        # handle parms for colorbars and colormaps
        vlim, ticklocs, ticklabels = get_contcbarargs(np.linspace(vlim[0], vlim[1], 5),
                                                         sigfigs, vlim)
        # Plot the colorbar
        cbar = plt.colorbar(plot, cax=cax, ticks=ticklocs)
        # Configure the color bar
        cbar.ax.set_yticklabels(ticklabels, ha='left')
        cbar.ax.tick_params(axis='y', pad=2)
        if cbar_label is not None:
            cbar.set_label(cbar_label, ha='center', va='top', labelpad=2)
    # Plot the dendrogram if required
    if hierarchy and dendrogram:
        # Get the line coordinates and scale
        ylines = np.array(dendo['dcoord'])
        ydendo = ny * 0.15
        ylines = ylines * (ydendo / ylines.max())
        xlines = np.array(dendo['icoord']) / 10
        lines = []
        for (xline, yline) in zip(xlines, ylines):
            lines.append(list(zip(xline, yline)))
        # Plot the dendrogram
        coll = mpl.collections.LineCollection(lines, color='k', lw=lw)
        ax_dendo.add_collection(coll)
        # Fix the subplots aesthetics
        ax_dendo.set_ylim(ydendo, 0)
        for spine in ax_dendo.spines:
            ax_dendo.spines[spine].set_visible(False)
        ax_dendo.yaxis.set_visible(False)
        ax_dendo.xaxis.set_visible(False)
        ax_dendo.patch.set_visible(False)
    # Export figure
    if output_file or ('pdfpages' in out_kws):
        export_image(output_file, **out_kws)

    return ax
