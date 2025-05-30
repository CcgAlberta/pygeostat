#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A collection of plotting tools to visualize bilabiate relationship between pairs of variables """
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters


@set_plot_style
def scatter_plot(x, y, wt=None, nmax=None, s=None, c=None, alpha=None, cmap=None, clim=None, cbar=False,
                cbar_label=None, stat_blk=None, stat_xy=None, stat_ha=None, stat_fontsize=None,
                roundstats=None, sigfigs=None, xlim=None, ylim=None, xlabel=None, ylabel=None, output_file=None, out_kws = None,
                title=None, grid=None, axis_xy=None, label='_nolegend_', ax=None, figsize=None,
                return_plot=False, logx=None, logy=None, **kwargs):
    '''
    Scatter plot that mimics the GSLIB scatter_plot program, providing summary statistics, kernel
    density estimate coloring, etc. NaN values are treated as null and removed from the plot and
    statistics.

    Parameters:
        x(np.ndarray or pd.Series): 1-D array with the variable to plot on the x-axis.
        y(np.ndarray or pd.Series): 1-D array with the variable to plot on the y-axis.

    Keyword arguments:
        wt(np.ndarray or pd.DataFrame): 1-D array with weights that are used in the calculation of
            displayed statistics.
        s(float or np.ndarray or pd.Series): size of each scatter point. Based on
            Parameters['plotting.scatter_plot.s'] if None.
        c(color or np.ndarray or pd.Series): color of each scatter point, as an array or valid
            Matplotlib color. Alternatively, 'KDE' may be specified to color each point according
            to its associated kernel density estimate. Based on Parameters['plotting.scatter_plot.c']
            if None.
        nmax (int): specify the maximum number of scatter points that should be displayed, which
            may be necessary due to the time-requirements of plotting many data. If specified,
            a nmax-length random sub-sample of the data is plotted. Note that this does not impact
            summary statistics.
        alpha(float): opacity of the scatter. Based on Parameters['plotting.scatter_plot.alpha'] if None.
        cmap (str): A matplotlib colormap object or a registered matplotlib
        clim (float tuple): Data minimum and maximum values
        cbar (bool): Indicate if a colorbar should be plotted or not
        cbar_label (str): Colorbar title
        stat_blk(str or list): statistics to place in the plot, which should be 'all' or
            a list that may contain ['count', 'pearson', 'spearman', 'noweightflag']. Based on
            Parameters['plotting.scatter_plot.stat_blk'] if None. Set to False to disable.
        stat_xy (float tuple): X, Y coordinates of the annotated statistics in figure
            space. Based on Parameters['plotting.scatter_plot.stat_xy'] if None.
        stat_ha (str): Horizontal alignment parameter for the annotated statistics. Can be
            ``'right'``, ``'left'``, or ``'center'``. If None, based on
            Parameters['plotting.stat_ha']
        stat_fontsize (float): the fontsize for the statistics block. If None, based on
            Parameters['plotting.stat_fontsize']. If less than 1, it is the fraction of the
            matplotlib.rcParams['font.size']. If greater than 1, it the absolute font size.
        roundstats (bool): Indicate if the statistics should be rounded to the number of digits or
            to a number of significant figures (e.g., 0.000 vs. 1.14e-5). The number of digits or
            figures used is set by the parameter ``sigfigs``. sigfigs (int): Number of significant
            figures or number of digits (depending on ``roundstats``) to display for the float
            statistics. Based on Parameters['plotting.roundstats'] and Parameters['plotting.roundstats']
            and Parameters['plotting.sigfigs'] if None.
        xlim(tuple): x-axis limits - xlim[0] to xlim[1]. Based on the data if None
        ylim(tuple): y-axis limits - ylim[0] to ylim[1]. Based on the data if None.
        xlabel(str): label of the x-axis, extracted from x if None
        ylabel(str): label of the y-axis, extracted from y if None
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
        title(str): plot title
        grid(bool): plot grid lines in each panel? Based on Parameters['plotting.grid'] if None.
        axis_xy(bool): if True, mimic a GSLIB-style scatter_plot, where only the bottom and left axes
            lines are displayed. Based on Parameters['plotting.axis_xy'] if None.
        label(str): label of scatter for legend
        ax(Matplotlib axis handle): if None, create a new figure and axis handles
        figsize(tuple): size of the figure, if creating a new one when ax = None
        logx, logy (str): permissible mpl axis scale, like `log`
        **kwargs: Optional permissible keyword arguments to pass to either: (1) matplotlib's
            scatter function

    Return:
        ax(Matplotlib axis handle)

    **Examples:**

    Basic scatter example:

    .. plot::

        import pygeostat as gs

        # Load the data
        data_file = gs.ExampleData('point3d_ind_mv')

        # Select a couple of variables
        x, y = data_file[data_file.variables[0]], data_file[data_file.variables[1]]

        # Scatter plot with default parameters
        gs.scatter_plot(x, y, figsize=(5, 5), cmap='hot')

        # Scatter plot without correlation and with a color bar:
        gs.scatter_plot(x, y, nmax=2000, stat_blk=False, cbar=True, figsize=(5, 5))

        # Scatter plot with the a constant color, transparency and all statistics
        # Also locate the statistics where they are better seen
        gs.scatter_plot(x, y, c='k', alpha=0.2, nmax=2000, stat_blk='all', stat_xy=(.95, .95),
                   figsize=(5, 5))
    '''
    # Import packages
    from scipy.stats import gaussian_kde
    from copy import deepcopy
    import pygeostat as gs
    from . utils import _set_stat_fontsize
    # Figure out the plotting axes
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    # Labels if present
    if xlabel is None:
        xlabel = gs.get_label(x)
    if ylabel is None:
        ylabel = gs.get_label(y)
    # Check the input data
    if isinstance(x, pd.DataFrame) or isinstance(x, pd.Series):
        x = x.values
    if x.ndim > 1:
        raise ValueError('x should be one-dimension!')
    if isinstance(y, pd.DataFrame) or isinstance(y, pd.Series):
        y = y.values
    if y.shape != x.shape:
        raise ValueError('x and y should be the same shape!')
    # Check the weights
    if isinstance(wt, pd.DataFrame) or isinstance(wt, pd.Series):
        wt = wt.values
    elif wt is None:
        wt = np.ones(x.shape)
    if wt.shape != x.shape:
        raise ValueError('x, y and wt should be the same shape!')
    # Remove nans if present
    idx = np.logical_and(np.isfinite(x), np.isfinite(y), np.isfinite(wt))
    x, y, wt = x[idx], y[idx], wt[idx]
    # Draw a random sub-sample if present
    xplot, yplot = deepcopy(x), deepcopy(y)
    if isinstance(nmax, int):
        if len(xplot) > nmax:
            idx1 = np.random.randint(0, len(xplot), nmax)
            xplot = xplot[idx1]
            yplot = yplot[idx1]
    else:
        idx1 = np.arange(0, len(xplot))
    # There's probably a lot of edge cases to this testing that are not yet
    # handled
    if isinstance(c, pd.DataFrame) or isinstance(c, pd.Series):
        if cbar_label is None:
            cbar_label = gs.get_label(c)
        c = c.values
    if isinstance(c, np.ndarray):
        c = c[idx]
        c = c[idx1]
    # Calculate kernel density estimate at data locations if necessary
    if c is None:
        c = Parameters['plotting.scatter_plot.c']
    kde = False
    if isinstance(c, str):
        if c.lower()[:3] == 'kde':
            pval = c.lower()[3:]
            # Points are colored based on KDE
            if logy:
                ykde = yplot.copy()
                ykde[ykde <= 0] = Parameters['plotting.log_lowerval']
                ykde = np.log(ykde)
            else:
                ykde = yplot
            if logx:
                xkde = xplot.copy()
                xkde[xkde <= 0] = Parameters['plotting.log_lowerval']
                xkde = np.log(xkde)
            else:
                xkde = xplot
            xy = np.stack((xkde, ykde), axis=1)
            kde = gaussian_kde(xy.T)
            c = kde.evaluate(xy.T)
            c = (c - min(c)) / (max(c) - min(c))
            if len(pval) > 0:
                try:
                    if pval.startswith('p'):
                        ipval = int(pval.lower()[1:])
                    else:
                        ipval = int(pval.lower())
                    assert (ipval <= 100) and (ipval >= 1)
                    ipval -= 1
                except ValueError:
                    raise ValueError('Could not interpret {} as a kde percentile!'.
                                     format(pval.lower()))
                except AssertionError:
                    raise ValueError('kde percentiles must be 1 <= p <= 100 ')
                cdfx, cdfy = gs.cdf(c, bins=101)
                clipval = np.interp(ipval / 100, cdfy, cdfx)
                c[c > clipval] = clipval
            kde = True
        else:
            cbar = False
    # Draw parameters from Parameters if necessary
    if s is None:
        s = Parameters['plotting.scatter_plot.s']
    if alpha is None:
        alpha = Parameters['plotting.scatter_plot.alpha']
    if stat_blk is None:
        stat_blk = Parameters['plotting.scatter_plot.stat_blk']
    if roundstats is None:
        roundstats = Parameters['plotting.roundstats']
    if sigfigs is None:
        sigfigs = Parameters['plotting.sigfigs']
    # Set-up some parameters
    if len(c) != xplot.shape[0]:
        cmap = False
    else:
        if cmap is None:
            cmap = Parameters['plotting.scatter_plot.cmap']
    if cmap is not False:
        clim, ticklocs, ticklabels = gs.get_contcbarargs(c, sigfigs, clim)
    if clim is None:
        clim = (None, None)
    # Set-up plot if no axis is supplied using the ImageGrid method if required or the regular way
    cax = None
    fig, ax, cax = gs.setup_plot(ax, cax=cax, cbar=cbar, figsize=figsize)
    # Scatter - let Matplotlib use the default size/color if None
    if s is None:
        if c is None:
            plot = ax.scatter(xplot, yplot, alpha=alpha, label=label, cmap=cmap,
                              vmin=clim[0], vmax=clim[1], **kwargs)
        else:
            plot = ax.scatter(xplot, yplot, c=c, alpha=alpha, label=label, cmap=cmap,
                              vmin=clim[0], vmax=clim[1], **kwargs)
    else:
        if c is None:
            plot = ax.scatter(xplot, yplot, s=s, alpha=alpha, label=label, cmap=cmap,
                              vmin=clim[0], vmax=clim[1], **kwargs)
        else:
            plot = ax.scatter(xplot, yplot, s=s, c=c, alpha=alpha, label=label, cmap=cmap,
                              vmin=clim[0], vmax=clim[1], **kwargs)
    # Setup the colorbar if required
    if cbar:
        if kde:
            if clim[0] is not None and clim[1] is not None:
                ticklocs = np.linspace(clim[0], clim[1], 3)
            else:
                ticklocs = [0, 0.5, 1]
            ticklabels = ['Low', 'Med.', 'High']
            cbar_label = 'Kernel Density Estimate'
        cbar = fig.colorbar(plot, cax=cax, ticks=ticklocs)
        # Configure the color bar
        cbar.ax.set_yticklabels(ticklabels, ha='left')
        cbar.ax.tick_params(axis='y', pad=2)
        if cbar_label is not None:
            cbar.set_label(cbar_label, ha='center', va='top', labelpad=2)
    # Set the axis extents
    if xlim is None:
        xlim = (np.min(x), np.max(x))
    if ylim is None:
        ylim = (np.min(y), np.max(y))
    if logx and xlim[0] <= 0:
        if xlim[0] == 0:
            xlim = [Parameters['plotting.log_lowerval'], ylim[1]]
        else:
            raise ValueError('ERROR: invalid clim for a log x-axis!')
    if logy and ylim[0] <= 0:
        if ylim[0] == 0:
            ylim = [Parameters['plotting.log_lowerval'], ylim[1]]
        else:
            raise ValueError('ERROR: invalid clim for a log y-axis!')
    # Set the formatting attributes
    gs.format_plot(ax, xlabel, ylabel, title, grid, axis_xy, xlim, ylim, logx, logy)
    # Setup the correlation
    if stat_blk:
        stats = ['pearson', 'spearmanr', 'count', 'noweightflag']
        # Error checking and conversion to a list of stats
        if isinstance(stat_blk, str):
            if stat_blk == 'all':
                stat_blk = stats[:-1]
            else:
                stat_blk = [stat_blk]
        elif isinstance(stat_blk, tuple):
            stat_blk = list(stat_blk)
        if isinstance(stat_blk, list):
            for stat in stat_blk:
                if stat not in stats:
                    raise ValueError('invalid stat_blk')
        else:
            raise ValueError('invalid stat_blk')
        # Build the txtstats
        txtstats = ''
        if 'count' in stat_blk:
            txtstats += r'$n = $'+str(x.shape[0])
        if 'pearson' in stat_blk:
            corr = gs.weighted_correlation(x, y, wt)
            if roundstats:
                corr = round(corr, sigfigs)
            else:
                corr = gs.round_sigfig(corr, sigfigs)
            txtstats += '\n'+r'$\rho = $'+str(corr)
        if 'spearmanr' in stat_blk:
            corr = gs.weighted_correlation_rank(x, y, wt)
            if roundstats:
                corr = round(corr, sigfigs)
            else:
                corr = gs.round_sigfig(corr, sigfigs)
            txtstats += '\n'+r'$\rho_s = $'+str(corr)
        # Note if weights were used
        if len(np.unique(wt)) > 1 and 'noweightflag' not in stat_blk:
            txtstats = txtstats + '\n\nweights used'
        # Sort the location and font size
        if stat_xy is None:
            stat_xy = Parameters['plotting.scatter_plot.stat_xy']
        if stat_ha is None:
            stat_ha = Parameters['plotting.stat_ha']
        if stat_xy[1] > 0.5:
            va = 'top'
        else:
            va = 'bottom'
        stat_fontsize = _set_stat_fontsize(stat_fontsize)
        # Draw to plot
        ax.text(stat_xy[0], stat_xy[1], txtstats, va=va, ha=stat_ha, transform=ax.transAxes,
                fontsize=stat_fontsize, linespacing=0.8)

    # Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()

    if output_file or ('pdfpages' in out_kws):
        gs.export_image(output_file, **out_kws)
        
    if return_plot:
        return ax, plot
    else:
        return ax


@set_plot_style
def scatter_plots(data, variables=None, wt=None, labels=None, nmax=None, pad=0.0, s=None, c=None,
             alpha=None, cmap=None, clim=None, cbar=True, cbar_label=None,
             stat_blk=None, stat_xy=None, stat_ha=None, stat_fontsize=None,
             roundstats=None, sigfigs=None,
             grid=None, axis_xy=None, xlim=None, ylim=None, label='_nolegend_', output_file = None, out_kws = None,
             figsize=None, **kwargs):
    '''
    Function which wraps the scatter_plot function, creating an upper matrix triangle of scatterplots
    for multiple variables.

    Parameters:
        data(np.ndarray or pd.DataFrame or gs.DataFile) : 2-D data array, which should be
            dimensioned as (ndata, nvar). Alternatively, specific variables may be selected
            with the variables argument. If a DataFile is passed and data.variables has a length
            greater than 1, those columns will be treated as the variables to plot.

    Keyword arguments:
        variables(str list): indicates the column names to treat as variables in data
        wt(np.ndarray or pd.Series or str or bool): array with weights
            that are used in the calculation of displayed statistics. Alternatively, a str may
            specify the weight column in lower. If data is a DataFile and data.wts is not None,
            then wt=True may be used to apply those weights.
        labels(tuple or nvar-list): labels for data, which are drawn from data if None
        nmax (int): specify the maximum number of scatter points that should be displayed, which
            may be necessary due to the time-requirements of plotting many data. If specified,
            a nmax-length random sub-sample of the data is plotted. Note that this does not impact
            summary statistics.
        pad(float or 2-tuple): space between each panel, which may be negative or positive. A tuple
            of (xpad, ypad) may also be used.
        align_orient(bool): align the orientation of plots in the upper and lower triangle (True),
            which causes the lower triangle plots to be flipped (x and y axes) from their
            standard symmetric orientation.
        titles(2-tuple str): titles of the lower and upper triangles (lower title, upper title)
        titlepads(2-tuple float): padding of the titles to the left of the lower triangle
            titlepads[0] and above the upper triangle (titlepads[1]). Typical required numbers
            are in the range of 0.01 to 0.5, depending on figure dimensioning.
        titlesize(int): size of the title font
        s(float or np.ndarray or pd.Series): size of each scatter point. Based on
            Parameters['plotting.scatter_plot.s'] if None.
        c(color or np.ndarray or pd.Series): color of each scatter point, as an array or valid
            Matplotlib color. Alternatively, 'KDE' may be specified to color each point according
            to its associated kernel density estimate. Based on Parameters['plotting.scatter_plot.c']
            if None.
        alpha(float): opacity of the scatter. Based on Parameters['plotting.scatter_plot.alpha'] if None.
        cmap(str): A matplotlib colormap object or a registered matplotlib
        clim(2-tuple float): Data minimum and maximum values
        cbar(bool): plot a colorbar for the color of the scatter (if variable)? (default=True)
        cbar_label(str): colorbar label(automated if KDE coloring)
        stat_blk(str or tuple): statistics to place in the plot, which should be 'all' or
            a tuple that may contain ['count', 'pearson', 'spearman']. Based on
            Parameters['plotting.scatter_plot.stat_blk'] if None. Set to False to disable.
        stat_xy(2-tuple float): X, Y coordinates of the annotated statistics in figure
            space. Based on Parameters['plotting.scatter_plot.stat_xy'] if None.
        stat_ha(str): Horizontal alignment parameter for the annotated statistics. Can be
            ``'right'``, ``'left'``, or ``'center'``. If None, based on
            Parameters['plotting.stat_ha']
        stat_fontsize(float): the fontsize for the statistics block. If None, based on
            Parameters['plotting.stat_fontsize']. If less than 1, it is the fraction of the
            matplotlib.rcParams['font.size']. If greater than 1, it the absolute font size.
        roundstats(bool): Indicate if the statistics should be rounded to the number of digits or
            to a number of significant figures (e.g., 0.000 vs. 1.14e-5). The number of digits or
            figures used is set by the parameter ``sigfigs``. sigfigs (int): Number of significant
            figures or number of digits (depending on ``roundstats``) to display for the float
            statistics. Based on Parameters['plotting.roundstats'] and Parameters['plotting.roundstats']
            and Parameters['plotting.sigfigs'] if None.
        grid(bool): plot grid lines in each panel? Based on Parameters['plotting.grid'] if None.
        axis_xy(bool): if True, mimic a GSLIB-style scatter_plot, where only the bottom and left axes
            lines are displayed. Based on Parameters['plotting.axis_xy'] if None.
        xlim(2-tuple float): x-axis limits - xlim[0] to xlim[1]. Based on the data if None
        ylim(2-tuple float): y-axis limits - ylim[0] to ylim[1]. Based on the data if None.
        label(str): label of scatter for legend
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
        figsize(2-tuple float): size of the figure, if creating a new one when ax = None
        return_handles(bool) : return figure handles? (default=False)
        **kwargs: Optional permissible keyword arguments to pass to either: (1) matplotlib's
            scatter function

    Return:
        matplotlib figure handle


    **Example:**

    Only one basic example is provided here, although all kwargs applying to the underlying scatter_plot
    function may be applied to scatter_plots.

    .. plot::

        import pygeostat as gs

        # Load the data, which registers the variables attribute
        data_file = gs.ExampleData('point3d_ind_mv')

        # Plot with the default KDE coloring
        fig = gs.scatter_plots(data_file, nmax=1000, stat_xy=(0.95, 0.95), pad=(-1, -1), s=10,
                          figsize=(10, 10))
    '''
    import pandas as pd
    import pygeostat as gs
    # Parse the data, variables and wt inputs, returning appropriate inputs
    data, wt, labels = _handle_variables_wt(data, variables, wt, labels)
    nvar = data.shape[1]
    # Iterate over the pairs
    fig, axes = plt.subplots(nvar-1, nvar-1, figsize=figsize)
    for i in np.arange(0, nvar-1):
        for j in np.arange(1, nvar):
                if i < j:
                    _, plot = scatter_plot(data.iloc[:, j], data.iloc[:, i], wt=wt, s=s, c=c,
                                      alpha=alpha, clim=clim, cmap=cmap, cbar=False, nmax=nmax,
                                      cbar_label=False, stat_blk=stat_blk, stat_xy=stat_xy,
                                      stat_fontsize=stat_fontsize, return_plot=True,
                                      stat_ha=stat_ha, roundstats=roundstats, sigfigs=sigfigs,
                                      xlim=xlim, ylim=ylim, ax=axes[i][j-1], xlabel=labels[j],
                                      ylabel=labels[i], grid=grid, axis_xy=axis_xy, **kwargs)
                    if i == j-1:
                        _tickoff(axes[i][j-1], xtickoff=False, ytickoff=False)
                    else:
                        _tickoff(axes[i][j-1], xtickoff=True, ytickoff=True)
                else:
                    axes[i][j-1].axis('off')
    try:
        fig.tight_layout(h_pad=pad[1], w_pad=pad[0])
    except:
        fig.tight_layout(h_pad=pad, w_pad=pad)
    # Figure out if KDE was used
    kde = False
    if c is None:
        c = Parameters['plotting.scatter_plot.c']
        if isinstance(c, str):
            if c.lower() == 'kde':
                kde = True
    if (not kde and not isinstance(c, pd.DataFrame) and not isinstance(c, np.ndarray) and
            not isinstance(c, pd.Series)):
        cbar = False
    # Colorbar
    if cbar:
        cbar_ax = fig.add_axes([0.2, .15, .03, .25])
        if kde:
            cbar = fig.colorbar(plot, cax=cbar_ax, ticks=[0, .5, 1])
            cbar.ax.set_yticklabels(['Low', 'Med.', 'High'])
            cbar.set_label('Kernel Density Estimate', ha='center', va='top', labelpad=2)
                #cbar.ax.set_title('KDE')
        else:
            cbar = fig.colorbar(plot, cax=cbar_ax)
            try:
                cbar_label = gs.get_label(c)
            except:
                pass
            if cbar_label is not None:
                cbar.set_label(cbar_label, ha='center', va='top', labelpad=2)

    # Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()
        
    if output_file or ('pdfpages' in out_kws):
        gs.export_image(output_file, **out_kws)
    return fig


@set_plot_style
def scatter_plots_lu(lower, upper, lower_variables=None, upper_variables=None, lowwt=None, uppwt=None,
                    lowlabels=None, upplabels=None, nmax=None,
                    pad=0.0, align_orient=False, titles=None, titlepads=None, titlesize=None,
                    s=None, c=None, alpha=None, cbar=True,
                    cbar_label=None, cmap=None, clim=None, stat_blk=None, stat_xy=None,
                    stat_ha=None, stat_fontsize=None, roundstats=None, sigfigs=None,
                    xlim=None, ylim=None, label='_nolegend_', output_file = None, out_kws = None,
                    grid=True, axis_xy=None, figsize=None, return_handle=False, **kwargs):
    '''
    Function which wraps the scatter_plot function, creating an upper/lower matrix triangle of
    scatterplots for comparing the scatter of multiple variables in two data sets.

    Parameters:
        lower(np.ndarray or pd.DataFrame or gs.DataFile): 2-D data array, which should be
            dimensioned as (ndata, nvar). Alternatively, specific variables may be selected
            with the variables argument. If a DataFile is passed and data.variables has a length
            greater than 1, those columns will be treated as the variables to plot.
            This data is plotted in the lower triangle.
        upper(np.ndarray or pd.DataFrame or gs.DataFile): see the description for lower, although
            this data is plotted in the upper triangle.

    Keyword arguments:
        lower_variables(nvar-tuple str): indicates the column names to treat as variables in lower
        upper_variables(nvar-tuple str): indicates the column names to treat as variables in upper
        lowwt(np.ndarray or pd.Series or str or bool): array with weights that are used in the
            calculation of displayed statistics for the lower data. Alternatively, a str may
            specify the weight column in lower. If lower is a DataFile and lower.wt is not None,
            then wt=True may be used to apply those weights.
        uppwt(np.ndarray or pd.DataFrame or str or bool): see the description for
            lowwt, although these weights are applied to upper.
        lowlabels(nvar-tuple str): labels for lower, which are drawn from lower if None
        upplabels(nvar-tuple str): labels for upper, which are drawn from upper if None
        nmax (int): specify the maximum number of scatter points that should be displayed, which
            may be necessary due to the time-requirements of plotting many data. If specified,
            a nmax-length random sub-sample of the data is plotted. Note that this does not impact
            summary statistics.
        pad(float or 2-tuple): space between each panel, which may be negative or positive. A tuple
            of (xpad, ypad) may also be used.
        align_orient(bool): align the orientation of plots in the upper and lower triangle (True),
            which causes the lower triangle plots to be flipped (x and y axes) from their
            standard symmetric orientation.
        titles(2-tuple str): titles of the lower and upper triangles (lower title, upper title)
        titlepads(2-tuple float): padding of the titles to the left of the lower triangle
            titlepads[0] and above the upper triangle (titlepads[1]). Typical required numbers
            are in the range of 0.01 to 0.5, depending on figure dimensioning.
        titlesize(int): size of the title font
        s(float or np.ndarray or pd.Series): size of each scatter point. Based on
            Parameters['plotting.scatter_plot.s'] if None.
        c(color or np.ndarray or pd.Series): color of each scatter point, as an array or valid
            Matplotlib color. Alternatively, 'KDE' may be specified to color each point according
            to its associated kernel density estimate. Based on Parameters['plotting.scatter_plot.c']
            if None.
        alpha(float): opacity of the scatter. Based on Parameters['plotting.scatter_plot.alpha'] if None.
        cmap(str): A matplotlib colormap object or a registered matplotlib
        clim(2-tuple float): Data minimum and maximum values
        cbar(bool): plot a colorbar for the color of the scatter (if variable)? (default=True)
        cbar_label(str): colorbar label(automated if KDE coloring)
        stat_blk(str or tuple): statistics to place in the plot, which should be 'all' or
            a tuple that may contain ['count', 'pearson', 'spearman']. Based on
            Parameters['plotting.scatter_plot.stat_blk'] if None. Set to False to disable.
        stat_xy(2-tuple float): X, Y coordinates of the annotated statistics in figure
            space. Based on Parameters['plotting.scatter_plot.stat_xy'] if None.
        stat_ha(str): Horizontal alignment parameter for the annotated statistics. Can be
            ``'right'``, ``'left'``, or ``'center'``. If None, based on
            Parameters['plotting.stat_ha']
        stat_fontsize(float): the fontsize for the statistics block. If None, based on
            Parameters['plotting.stat_fontsize']. If less than 1, it is the fraction of the
            matplotlib.rcParams['font.size']. If greater than 1, it the absolute font size.
        roundstats(bool): Indicate if the statistics should be rounded to the number of digits or
            to a number of significant figures (e.g., 0.000 vs. 1.14e-5). The number of digits or
            figures used is set by the parameter ``sigfigs``. sigfigs (int): Number of significant
            figures or number of digits (depending on ``roundstats``) to display for the float
            statistics. Based on Parameters['plotting.roundstats'] and Parameters['plotting.roundstats']
            and Parameters['plotting.sigfigs'] if None.
        grid(bool): plot grid lines in each panel? Based on Parameters['plotting.grid'] if None.
        axis_xy(bool): if True, mimic a GSLIB-style scatter_plot, where only the bottom and left axes
            lines are displayed. Based on Parameters['plotting.axis_xy'] if None.
        xlim(2-tuple float): x-axis limits - xlim[0] to xlim[1]. Based on the data if None
        ylim(2-tuple float): y-axis limits - ylim[0] to ylim[1]. Based on the data if None.
        label(str): label of scatter for legend
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
        figsize(2-tuple float): size of the figure, if creating a new one when ax = None
        return_handles(bool) : return figure handles? (default=False)
        **kwargs: Optional permissible keyword arguments to pass to either: (1) matplotlib's
            scatter function

    Return:
        matplotlib figure handle

    **Examples:**

    Plot with varying orientations that provide correct symmetry (above) and ease of comparison
    (below). Here, the data is treated as both the data and a realization (first two arguments)
    for the sake of demonstration.

    .. plot::

        import pygeostat as gs
        import numpy as np

        # Load the data, which registers the variables attribute
        data_file1 = gs.ExampleData('point3d_ind_mv')
        data_file2 = gs.ExampleData('point3d_ind_mv')
        mask = np.random.rand(len(data_file2))<0.3
        data_file2.data = data_file2.data[mask]

        # Plot with the standard orientation
        fig = gs.scatter_plots_lu(data_file1, data_file2, titles=('Data', 'Realization'), s=10, nmax=1000,
                             stat_xy=(0.95, 0.95), pad=(-1, -1), figsize=(10, 10))

        # Plot with aligned orientation to ease comparison
        fig = gs.scatter_plots_lu(data_file1, data_file2, titles=('Data', 'Realization'), s=10, nmax=1000,
                             stat_xy=(0.95, 0.95), pad=(-1, -1), figsize=(10, 10), cmap='jet',
                             align_orient=True)

    '''
    import pygeostat as gs
    import matplotlib as mpl
    # Parse the data, variables and wt inputs, returning appropriate inputs
    lower, lowwt, lowlabels = _handle_variables_wt(lower, lower_variables, lowwt, lowlabels)
    upper, uppwt, upplabels = _handle_variables_wt(upper, upper_variables, uppwt, upplabels)
    nvar = upper.shape[1]
    if lower.shape[1] != nvar:
        raise ValueError('upper and lower were coerced into differing number of variables!')
    # Iterate over the pairs
    fig, axes = plt.subplots(nvar, nvar, figsize=figsize)
    for i in range(nvar):
        axes[i][i].axis('off')
        for j in range(nvar):
            if i < j:
                _, plot = scatter_plot(upper.iloc[:, j], upper.iloc[:, i], wt=uppwt, s=s, c=c, nmax=nmax,
                                  alpha=alpha, clim=clim, cmap=cmap, cbar=False, cbar_label=False,
                                  stat_blk=stat_blk, stat_xy=stat_xy, stat_fontsize=stat_fontsize,
                                  stat_ha=stat_ha, roundstats=roundstats, sigfigs=sigfigs,
                                  xlim=xlim, ylim=ylim, ax=axes[i][j], xlabel=upplabels[j],
                                  ylabel=upplabels[i], grid=grid, axis_xy=axis_xy,
                                  return_plot=True, **kwargs)
                if i + 1 == j:
                    _tickoff(axes[i][j], xtickoff=False, ytickoff=False)
                else:
                    _tickoff(axes[i][j], xtickoff=True, ytickoff=True)
            elif i > j and align_orient:
                _, plot = scatter_plot(lower.iloc[:, i], lower.iloc[:, j], wt=lowwt, s=s, c=c, nmax=nmax,
                                  alpha=alpha, clim=clim, cmap=cmap, cbar=False, cbar_label=False,
                                  stat_blk=stat_blk, stat_xy=stat_xy, stat_fontsize=stat_fontsize,
                                  stat_ha=stat_ha, roundstats=roundstats, sigfigs=sigfigs,
                                  xlim=xlim, ylim=ylim, ax=axes[i][j], xlabel=lowlabels[i],
                                  ylabel=lowlabels[j], grid=grid, axis_xy=axis_xy,
                                  return_plot=True, **kwargs)
                if j == 0 and i == nvar-1:
                    _tickoff(axes[i][j], xtickoff=False, ytickoff=False)
                elif i == nvar-1:
                    _tickoff(axes[i][j], xtickoff=False, ytickoff=True)
                elif j == 0:
                    _tickoff(axes[i][j], xtickoff=True, ytickoff=False)
                else:
                    _tickoff(axes[i][j], xtickoff=True, ytickoff=True)
            elif i > j and not align_orient:
                _, plot = scatter_plot(lower.iloc[:, j], lower.iloc[:, i], wt=lowwt, s=s, c=c, nmax=nmax,
                                  alpha=alpha, clim=clim, cmap=cmap, cbar=False, cbar_label=False,
                                  stat_blk=stat_blk, stat_xy=stat_xy, stat_fontsize=stat_fontsize,
                                  stat_ha=stat_ha, roundstats=roundstats, sigfigs=sigfigs,
                                  xlim=xlim, ylim=ylim, ax=axes[i][j], xlabel=lowlabels[j],
                                  ylabel=lowlabels[i], grid=grid, axis_xy=axis_xy,
                                  return_plot=True, **kwargs)
                if j == 0 and i == nvar-1:
                    _tickoff(axes[i][j], xtickoff=False, ytickoff=False)
                elif i == nvar-1:
                    _tickoff(axes[i][j], xtickoff=False, ytickoff=True)
                elif j == 0:
                    _tickoff(axes[i][j], xtickoff=True, ytickoff=False)
                else:
                    _tickoff(axes[i][j], xtickoff=True, ytickoff=True)
    try:
        fig.tight_layout(h_pad=pad[1], w_pad=pad[0])
    except:
        fig.tight_layout(h_pad=pad, w_pad=pad)
    fig.subplots_adjust(top=.95, right=.95, left=.07)
    if titles is not None:
        if len(titles) != 2:
            raise ValueError('titles should be a 2-list of strings!')
        if titlepads is not None:
            if titlepads[0] is None:
                titlepads[0] = 3.*fig.dpi
            if titlepads[1] is None:
                titlepads[1] = 0.0
        else:
            titlepads = (0.08*fig.dpi, 0.01)
        if titlesize is None:
            titlesize = mpl.rcParams['font.size']
        gs.supaxislabel('y', titles[0], label_prop={'weight': 'bold', 'fontsize': titlesize},
                        fig=fig, labelpad=titlepads[0])
        fig.suptitle(titles[1], weight='bold', fontsize=titlesize, y=0.98+titlepads[1])

    # Figure out if KDE was used
    kde = False
    if c is None:
        c = Parameters['plotting.scatter_plot.c']
        if isinstance(c, str):
            if c.lower() == 'kde':
                kde = True
    if (not kde and not isinstance(c, pd.DataFrame) and not isinstance(c, np.ndarray) and
            not isinstance(c, pd.Series)):
        cbar = False
    # Colorbar
    if cbar:
        fig.subplots_adjust(bottom=.15)
        #ax = fig.add_axes([0.07, .15, .88, .8])
        cbar_ax = fig.add_axes([0.1, 0.04, 0.8, 0.02])
        if kde:
            cbar = fig.colorbar(plot, cax=cbar_ax, ticks=[0, .5, 1], orientation='horizontal')
            cbar.ax.set_xticklabels(['Low', 'Med.', 'High'])
            cbar.ax.set_title('Kernel Density Estimate')
                #cbar.ax.set_title('KDE')
        else:
            cbar = fig.colorbar(plot, cax=cbar_ax, orientation='horizontal')
            try:
                cbar_label = gs.get_label(c)
            except:
                pass
            if cbar_label is not None:
                cbar.ax.set_title(cbar_label)
    else:
        fig.subplots_adjust(bottom=.05)

    # Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()
        
    if output_file or ('pdfpages' in out_kws):
        gs.export_image(output_file, **out_kws)
        
    return fig


def _handle_variables_wt(data, variables, wt, labels):
    '''Given data, variables, wt and labels input, return data as a DataFrame of (ndata, nvar)
        dimension and wt of (ndata) dimension.'''
    import pygeostat as gs
    # Weight handling
    if isinstance(wt, str):
        if isinstance(data, gs.DataFile) or isinstance(data, pd.DataFrame):
            wt = data[wt].values
        else:
            raise ValueError(('wt as column specifier is only valid if data is a DataFile'
                              ' or DataFrame'))
    elif isinstance(wt, bool):
        if isinstance(data, gs.DataFile):
            if wt:
                if data.wts is None:
                    raise ValueError('wt=True is only valid if data.wts is not None')
                wt = data[data.wts].values
        else:
            raise ValueError('wt as a boolean is only valid if data is a DataFile')
    elif wt is not None:
        raise ValueError('invalid wt type!')
    if wt is not None:
        if wt.ndim > 1:
            if wt.shape[0] > 1:
                raise ValueError('wt must be 1D!')
            else:
                wt = wt.flatten()
    # Variable handling
    if isinstance(variables, list) or isinstance(variables, tuple):
        if isinstance(data, gs.DataFile) or isinstance(data, pd.DataFrame):
            data = data[variables]
        else:
            raise ValueError('variables is only valid if data is a DataFile or DataFrame')
    elif isinstance(data, gs.DataFile):
        if isinstance(data.variables, list):
            data = data[data.variables]
        else:
            data = data.data
    if not isinstance(data, pd.DataFrame):
        try:
            data = pd.DataFrame(data)
        except:
            raise ValueError('could not coerce provided data into a pandas DataFrame!')
        nvar = data.shape[1]
        data.columns = ['Var'+str(i+1) for i in range(nvar)]
    else:
        nvar = data.shape[1]
    if data.shape[1] < 3:
        raise ValueError('nvar < 3 is invalid - use scatter_plot for plotting two variables!')
    # Check the weights now that variables are assembled
    if wt is not None:
        if len(wt) != data.shape[0]:
            raise ValueError('wt does not have the same number of obervations as data!')

    if labels is None:
        labels = data.columns
    return data, wt, labels


def _tickoff(ax, xtickoff, ytickoff):
    '''Remove the xtick and/or ytick labels from the an axis handle'''
    if xtickoff:
        ax.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False)
        ax.set_xlabel('')
    if ytickoff:
        ax.tick_params(
            axis='y',
            which='both',
            left=False,
            right=False,
            labelleft=False)
        ax.set_ylabel('')
