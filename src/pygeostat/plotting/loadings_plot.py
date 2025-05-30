#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Generate a loadings plot depicting the loadings or correlation between the original variables and
their transformed counterparts."""

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
def loadings_plot(loadmat, figsize=None, ax=None, title=None, xticklabels=None, yticklabels=None,
                rotateticks=None, plot_style=None, custom_style=None, output_file=None, **kwargs):
    """
    This function uses matplotlib to create a loadings plot with variably sized colour mapped
    boxes illustrating the contribution of each of the input variables to the transformed variables.

    The only parameter needed ``loadmat`` containing the loadings or correlation matrix. All of the
    other arguments are optional. Figure size will likely have to be manually adjusted. If
    ``xticklabels`` and/or ``yticklabels`` are left to their default value of ``None`` and the input
    matrix is contained in a pandas dataframe, the index/column information will be used to label
    the columns and rows. If a numpy array is passed, axis tick labels will need to be provided.
    Axis tick labels are automatically checked for overlap and if needed, are rotated. If rotation
    is necessary, consider condensing the variable names or plotting a larger figure if the result
    appears odd.

    Please review the documentation of the :func:`gs.set_style()
    <pygeostat.plotting.set_style.set_style>` and :func:`gs.export_image()
    <pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
    their use in this function can be understood.

    Parameters:
        loadmat: Pandas dataframe or numpy matrix containing the required loadings or correlation
            matrix
        figsize (tuple): Figure size (width, height)
        ax (mpl.axis): Matplotlib axis to plot the figure
        title (str): Title for the plot.
        xticklabels (list): Tick labels along the x-axis
        yticklabels (list): Tick labels along the y-axis
        rotateticks (bool tuple): Indicate if the axis tick labels should be rotated (x, y)
        plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
            to turn it off
        custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
        output_file (str): Output figure file name and location
        **kwargs: Optional permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`

    Returns:
        ax (ax): matplotlib Axes object with the loadings plot


    Examples:

	Grab the correlation between the PCA variables and their corresponding input variables as
	a pandas dataframe:

	.. plot::

		import pygeostat as gs
		data_file = gs.ExampleData('3d_correlation')
		loadmat = data_file.data.corr().iloc[3:6,6:9] 
		gs.loadings_plot(loadmat.values, figsize=(5,5), xticklabels=['PC1', 'PC2', 'PC3'], yticklabels=['InputVariable1', 'InputVariable2', 'InputVariable3']) 

    """
    from .utils import _tickoverlap
    from .export_image import export_image

    # Set-up plot if no axis is supplied
    if ax is None:
        _, ax = plt.subplots(1, figsize=figsize)
    ax.set_aspect('equal')
    # Set title if required
    if title:
        ax.set_title(title)
    # Set the labels if possible
    if xticklabels and not isinstance(xticklabels, bool):
        xlabels = xticklabels
    elif hasattr(loadmat, 'columns'):
        xlabels = list(loadmat.columns)
    else:
        xlabels = None
    if yticklabels and not isinstance(yticklabels, bool):
        ylabels = yticklabels
    elif hasattr(loadmat, 'index'):
        ylabels = list(loadmat.index)
    else:
        ylabels = None
    # Set-up some parameters
    nx = loadmat.shape[1]
    ny = loadmat.shape[0]
    # Set-up colour mapping
    clrnorm = mpl.colors.Normalize(vmin=-1, vmax=1)
    clrmap = mpl.cm.ScalarMappable(norm=clrnorm, cmap='bwr')
    # Set-up axis parameters
    ax.set(xlim=(0, nx), ylim=(0, ny))
    xticklocs = np.arange(nx)
    yticklocs = np.arange(ny)
    # Set-up x axis labels and grid
    ax.xaxis.tick_top()
    ax.set_xticks(xticklocs + 0.5)
    ax.set_xticks(xticklocs, minor=True)
    if xlabels is not None:
        ax.set_xticklabels(xlabels, va='bottom', ha='center', rotation='horizontal')
        ax.tick_params(axis='x', pad=2)
    else:
        ax.get_xaxis().set_ticks([])
    # Set-up y axis labels and grid
    ax.invert_yaxis()
    ax.set_yticks(yticklocs + 0.5)
    ax.set_yticks(yticklocs, minor=True)
    if ylabels is not None:
        ax.set_yticklabels(ylabels, va='center', ha='right', rotation='vertical')
        ax.tick_params(axis='y', pad=2)
    else:
        ax.get_yaxis().set_ticks([])
    # The plots tick labels will not be properly accessible until the figure is "drawn", once the
    # command below is run, ax.get_ticklabel() will actually work properly.
    plt.draw()
    # Check if the axis tick labels overlap, if so rotate them
    if rotateticks is None:
        rotateticks = Parameters['plotting.rotateticks']
    if rotateticks is None:
        rotateticks = _tickoverlap(ax)
    # Rotate if required
    if rotateticks[0]:
        xlabels = ax.get_xticklabels()
        for xlabel in xlabels:
            xlabel.set_ha('center')
            xlabel.set_va('bottom')
            xlabel.set_rotation(45)
    if rotateticks[1]:
        ylabels = ax.get_yticklabels()
        for ylabel in ylabels:
            ylabel.set_ha('right')
            ylabel.set_va('center')
            ylabel.set_rotation(-45)
    # Plot grid
    plt.grid(False)
    plt.grid(True, which='minor')
    # Convert pandas dataframe to a numpy matrix
    if isinstance(loadmat, pd.DataFrame):
        loadmat = loadmat.values
    # Plot the loading blocks
    for y in range(ny):
        for x in range(nx):
            corr = loadmat[y, x]
            length = abs(corr) * 0.8
            boxstart = (0.5 - (length / 2))
            rectangle = mpl.patches.Rectangle((boxstart + x, (y + 0.15)), length, 0.35)
            rectangle.set(edgecolor='grey', facecolor=clrmap.to_rgba(corr))
            ax.annotate(round(corr, 2), xy=(x + .5, y + 0.8), ha='center')
            ax.add_patch(rectangle)
    # Export figure
    if output_file is not None:
        export_image(output_file, **kwargs)

    return ax
