#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""validation_plot plots a scatter plot of cross validation results with relevant statistics"""

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from .utils import get_label, get_statblk
from .export_image import export_image
from ..statistics.kde import kde_scipy
from . set_style import set_plot_style


@set_plot_style
def validation_plot(x, y, figsize=None, vlim=None, xlabel=None, ylabel=None, title=None, stat_blk='all',
                    stat_xy=(0.95, 0.05), stat_ha=None, stat_fontsize=None, mc='k', ms=None, plot_style=None,
                    lw=None, grid=None, axis_xy=None, custom_style=None, output_file=None, ax=None, dens=False, 
                    rasterized=False, **kwargs):
    """
    This function uses numpy to calculate the regression model and matplotlib to plot the
    scatter plot, regression line, and 45 degree line. Statistics are calculated using numpy.

    The only parameters needed are the ``x`` and ``y``. All of the other arguments are optional. If
    the label parameters are left to their default value of ``None``, the column information will
    be used to label the axes.

    Two statistics block sets are available: ``'minimal'`` and the default ``'all'``. The
    statistics block can be customized to a user defined list and order. Available statistics are
    as follows:

    >>> ['ndat', 'ymean', 'ystdev', 'xmean', 'xstdev', 'cov', 'rho', 'mse', 'sor']

    Please review the documentation of the :func:`gs.set_style()
    <pygeostat.plotting.set_style.set_style>` and :func:`gs.export_image()
    <pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
    their use in this function can be understood.

    Parameters:
        x: Tidy (long-form) 1D data where a single column of the variable to plot along the x-axis
            exists with each row is an observation. A pandas dataframe/series or numpy array can be
            passed.
        y: Tidy (long-form) 1D data where a single column of the variable to plot along the y-axis
            exists with each row is an observation. A pandas dataframe/series or numpy array can be
            passed.
        figsize (tuple): Figure size (width, height)
        vlim (float tuple): A single tuple for the minimum and maximum limits of data along both
            axes. Will not be a symmetrical plot if they are not the same value
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        title (str): Title for the plot
        stat_blk (str or list): Indicate what preset statistics block to write or a specific list
        stat_xy (str or float tuple): X, Y coordinates of the annotated statistics in figure
            space.
        stat_ha (str): Horizontal alignment parameter for the annotated statistics. Can be
            ``'right'``, ``'left'``, or ``'center'``. The value ``None`` can also be used to allow
            the parameter ``stat_xy`` to determine the alignment automatically.
        stat_fontsize (float): the fontsize for the statistics block. If None, based on
            gsParams['plotting.stat_fontsize']. If less than 1, it is the fraction of the
            matplotlib.rcParams['font.size']. If greater than 1, it the absolute font size.
        mc (str): Any permissible matplotlib color value for the scatter plot markers
        ms (float): Size of scatter plot markers
        grid(bool): plots the major grid lines if True. Based on gsParams['plotting.grid']
            if None.
        axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
            visible) if axis_xy is True. Based on gsParams['plotting.axis_xy'] if None.
        plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
            to turn it off
        custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
        output_file (str): Output figure file name and location
        **kwargs: Optional permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`

    Returns:
        ax (ax): Matplotlib Axes object with the cross validation plot

    **Examples:**
    
        A simple call:

        .. plot::

            import pygeostat as gs
            data = gs.ExampleData('3d_estimate')
            gs.validation_plot(data.data['Estimate'], data.data['True'], stat_blk='minimal')

        |

        Fixing the value limits, moving the statistics block, and exporting the figure.

        .. plot::

            import pygeostat as gs
            import numpy as np
            mean = [0, 0]
            cov = [[1, 0.8], [0.8, 1]]  # diagonal covariance
            x, y = np.random.multivariate_normal(mean, cov, 5000).T
            gs.validation_plot(x,y,vlim=(-3.5, 3.5) ,grid=True, stat_xy=(1, 0.68))

    |

    .. codeauthor:: pygeostat development team 2015-08-05
    """
    import types
    from .utils import format_plot, _set_stat_fontsize
    # Convert list to pd.series
    if isinstance(x, list):
        x = pd.DataFrame(x)
        x = x[0]
    if isinstance(y, list):
        y = pd.DataFrame(y)
        y = y[0]
    if isinstance(x, pd.DataFrame):
        x = x[x.columns[0]]
    if isinstance(y, pd.DataFrame):
        y = y[y.columns[0]]
    # Set-up figure
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    # Configure plot
    ax.axis('equal')
    if xlabel is None:
        xlabel = get_label(x)
        if xlabel == 0:
            xlabel = None
        if xlabel is None:
            xlabel = 'Estimate'
    if ylabel is None:
        ylabel = get_label(y)
        if ylabel == 0:
            ylabel = None
        if ylabel is None:
            ylabel = 'True'
    ax = format_plot(ax, xlabel=xlabel, ylabel=ylabel, title=title, grid=grid, axis_xy=axis_xy, xlim=vlim, ylim=vlim)
    #ax.tick_params(axis='both', pad=2)
    if isinstance(dens, types.FunctionType):  # fucntion takes (x, y, #comps) as positional args and returns `colors`
        mc = dens(x, y, 8)
    else:
        if dens:
            # plot the points as density colors
            values = np.vstack([x, y])
            mc = kde_scipy(values, values)
    if lw is None:
        lw = 0
    # Plot the scatter plot
    ax.scatter(x, y, c=mc, s=ms, lw=lw, rasterized=rasterized)
    # Get the axis limits determined by plt.scatter()
    if ax.get_ylim()[1] >= ax.get_xlim()[1]:
        maxlim = ax.get_ylim()
    else:
        maxlim = ax.get_xlim()
    # Plot the regression model
    regmodel = np.polyfit(x, y, 1)
    predicted = np.polyval(regmodel, [maxlim[0], maxlim[1]])
    ax.plot([maxlim[0], maxlim[1]], predicted, '-', lw=0.5, color='red')
    # Plots stats
    if stat_blk:
        ndat = len(x)
        ymean = np.nanmean(y)
        ystdev = np.nanstd(y)
        xmean = np.nanmean(x)
        xstdev = np.nanstd(x)
        cov = np.cov(x, y)[0, 1]
        rho = stats.pearsonr(x, y)[0]
        mse = np.sqrt(((x - y) ** 2).mean())
        sor = stats.linregress(x, y)[0]
        statlist = {'ndat': ('$n = {}$'.format(ndat)),
                    'ymean': (r'$m_y = {:0.3f}$'.format(ymean)),
                    'ystdev': (r'$\sigma_y = {:0.3f}$'.format(ystdev)),
                    'xmean': (r'$m_x = {:0.3f}$'.format(xmean)),
                    'xstdev': (r'$\sigma_x = {:0.3f}$'.format(xstdev)),
                    'cov': ('$cov = {:0.3f}$'.format(cov)),
                    'rho': (r'$\rho = {:0.3f}$'.format(rho)),
                    'mse': ('$RMSE = {:0.3f}$'.format(mse)),
                    'sor': ('$SoR = {:0.3f}$'.format(sor)),
                    'newline': ''}
        # Default statistic sets
        statsets = {'minimal': ['ndat', 'rho', 'mse', 'sor'],
                    'all': ['ndat', 'newline', 'ymean', 'ystdev', 'newline', 'xmean', 'xstdev',
                            'newline', 'cov', 'rho', 'mse', 'sor'],
                    'none': None}
        txtstats, stat_xy, ha, va = get_statblk(stat_blk, statsets, statlist, stat_xy)
        if stat_ha:
            ha = stat_ha
        stat_fontsize = _set_stat_fontsize(stat_fontsize)
        ax.text(stat_xy[0], stat_xy[1], txtstats, va=va, ha=ha, fontsize=stat_fontsize,
                transform=ax.transAxes)
    # Plot the 45 degree line
    ax.plot(maxlim, maxlim, lw=0.5, color='k')
    ax.set(xlim=maxlim, ylim=maxlim)  # Reset the axis limits as the 45 degree line pushed it out

    # Export figure
    if output_file or ('pdfpages' in kwargs):
        export_image(output_file, **kwargs)

    return ax
