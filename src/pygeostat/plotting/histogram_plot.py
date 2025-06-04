#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A basic histogram plotting routine using matplotlib"""

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters

@set_plot_style
def histogram_plot(data, var=None, weights=None, cat=None, catdict=None, bins=None, icdf=False, lower=None,
                    upper=None, ax=None, figsize=None, xlim=None, ylim=None, title=None, xlabel=None,
                    stat_blk=None, stat_xy=None, stat_ha=None, roundstats=None, sigfigs=None, color=None,
                    edgecolor=None, edgeweights=None, grid=None, axis_xy=None, label_count=False,
                    rotateticks=None, plot_style=None, custom_style=None, output_file=None, out_kws=None,
                    stat_fontsize=None, stat_linespacing=None, logx=False, **kwargs):
    """
    Generates a matplotlib style histogram with summary statistics. Trimming is now only applied
    to NaN values (Pygeostat null standard).

    The only required required parameter is ``data``. If ``xlabel`` is left to its default value of
    ``None`` and the input data is contained in a pandas dataframe or series, the column
    information will be used to label the x-axis.

    Two statistics block sets are available: ``'all'`` and the default ``'minimal'``. The
    statistics block can be customized to a user defined list and order. Available statistics are
    as follows:

    >>> ['count', 'mean', 'stdev', 'cvar', 'max', 'upquart', 'median', 'lowquart', 'min',
    ...  'p10', 'p90']

    The way in which the values within the statistics block are rounded and displayed can be
    controlled using the parameters ``roundstats`` and ``sigfigs``.

    Please review the documentation of the :func:`gs.set_style()
    <pygeostat.plotting.set_style.set_style>` and :func:`gs.export_image()
    <pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
    their use in this function can be understood.

    Parameters:
        data (np.ndarray, pd.DataFrame/Series, or gs.DataFile): data array, which must be 1D
            unless var is provided. The exception being a DataFile, if data.variables
            is a single name.
        var (str): name of the variable in data, which is required if data is not 1D.
        weights (np.ndarray, pd.DataFrame/Series, or gs.DataFile or str): 1D array of declustering
             weights for the data. Alternatively the declustering weights name in var. If data
             is a DataFile, it may be string in data.columns, or True to use data.weights
             (if data.weights is not None).
        cat (bool or str): either a cat column in data.data, or if True uses data.cat if data.cat
            is not None
        catdict (dict or bool): overrides bins. If a categorical variable is being plotted, provide
            a dictionary where keys are numeric (categorical codes) and values are their associated
            labels (categorical names). The bins will be set so that the left edge (and associated
            label) of each bar is inclusive to each category. May also be set to True, if data is
            a DataFile and data.catdict is initialized.
        bins (int or list): Number of bins to use, or a list of bins
        icdf (bool): Indicator to plot a CDF or not
        lower (float): Lower limit for histogram
        upper (float): Upper limit for histogram
        ax (mpl.axis): Matplotlib axis to plot the figure
        figsize (tuple): Figure size (width, height)
        xlim (float tuple): Minimum and maximum limits of data along the x axis
        ylim (float tuple): Minimum and maximum limits of data along the y axis
        title (str): Title for the plot
        xlabel (str): X-axis label
        stat_blk (bool): Indicate if statistics are plotted or not
        stat_xy (float tuple): X, Y coordinates of the annotated statistics in figure
            space. Based on Parameters['plotting.histogram_plot.stat_xy'] if a histogram and
            Parameters['plotting.histogram_plot.stat_xy'] if a CDF, which defaults to the top right when
            a PDF is plotted and the bottom right if a CDF is plotted.
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
            statistics
        color (str or int or dict): Any permissible matplotlib color or a integer which is used to draw
            a color from the pygeostat color pallet ``pallet_pastel``> May also be a dictionary of colors,
            which are used for each bar (useful for categories). colors.keys() must align with bins[:-1]
            if a dictionary is passed. Drawn from Parameters['plotting.cmap_cat'] if catdict is used
            and their keys align.
        edgecolor (str): Any permissible matplotlib color for the edge of a histogram bar
        grid(bool): plots the major grid lines if True. Based on Parameters['plotting.grid']
            if None.
        axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
            visible) if axis_xy is True. Based on Parameters['plotting.axis_xy'] if None.
        label_count (bool): label the number of samples found for each category in catdict. Does
            nothing if no catdict is found
        rotateticks (bool tuple): Indicate if the axis tick labels should be rotated (x, y)
        plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
            to turn it off
        custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
        **kwargs: Optional permissible keyword arguments to pass to either: (1) matplotlib's hist
            function if a PDF is plotted or (2) matplotlib's plot function if a CDF is plotted.

    Returns:
        ax (ax): matplotlib Axes object with the histogram

    **Examples:**

    A simple call:

    .. plot::

        import pygeostat as gs
        # load some data
        dfl = gs.ExampleData("point3d_ind_mv")
        # plot the histogram_plot
        gs.histogram_plot(dfl, var="Phi", bins=30)

    |

    Change the colour, number of significant figures displayed in the statistics, and pass some
    keyword arguments to matplotlibs hist function:

    .. plot::

        import pygeostat as gs
        # load some data
        dfl = gs.ExampleData("point3d_ind_mv")
        # plot the histogram_plot
        gs.histogram_plot(dfl, var="Phi", color='#c2e1e5', sigfigs=5, log=True, density=True)

    |

    Plot a CDF while also displaying all available statistics, which have been shifted up:

    .. plot::

        import pygeostat as gs
        # load some data
        dfl = gs.ExampleData("point3d_ind_mv")
        # plot the histogram_plot
        gs.histogram_plot(dfl, var="Phi", icdf=True, stat_blk='all', stat_xy=(1, 0.75))
        # Change the CDF line colour by grabbing the 3rd colour from the color pallet
        # ``cat_vibrant`` and increase its width by passing a keyword argument to matplotlib's
        # plot function. Also define a custom statistics block:
        gs.histogram_plot(dfl, var="Phi", icdf=True, color=3, lw=3.5, stat_blk=['count','upquart'])

    |

    Generate histograms of Phi considering the categories:

    .. plot::

        import pygeostat as gs
        # load some data
        dfl = gs.ExampleData("point3d_ind_mv")
        cats = [1, 2, 3, 4, 5]
        colors = gs.catcmapfromcontinuous("Spectral", 5).colors
        # build the required cat dictionaries
        dfl.catdict = {c: "RT {:02d}".format(c) for c in cats}
        colordict =  {c: colors[i] for i, c in enumerate(cats)}
        # plot the histogram_plot
        f, axs = plt.subplots(2, 1, figsize=(8, 6))
        for var, ax in zip(["Phi", "Sw"], axs):
            gs.histogram_plot(dfl, var=var, cat=True, color=colordict, bins=40, figsize=(8, 4), ax=ax,
                       xlabel=False, title=var)

    |

    Generate cdf subplots considering the categories:

    .. plot::

        import pygeostat as gs
        # load some data
        dfl = gs.ExampleData("point3d_ind_mv")
        cats = [1, 2, 3, 4, 5]
        colors = gs.catcmapfromcontinuous("Spectral", 5).colors
        # build the required cat dictionaries
        dfl.catdict = {c: "RT {:02d}".format(c) for c in cats}
        colordict =  {c: colors[i] for i, c in enumerate(cats)}
        # plot the histogram_plot
        f, axs = plt.subplots(2, 2, figsize=(12, 9))
        axs=axs.flatten()
        for var, ax in zip(dfl.variables, axs):
            gs.histogram_plot(dfl, var=var, icdf=True, cat=True, color=colordict, ax=ax)

    Recreate the `Proportion` class plot

    .. plot::

        import pygeostat as gs
        # load some data
        dfl = gs.ExampleData("point3d_ind_mv")
        cats = [1, 2, 3, 4, 5]
        colors = gs.catcmapfromcontinuous("Spectral", 5).colors
        # build the required cat dictionaries
        dfl.catdict = {c: "RT {:02d}".format(c) for c in cats}
        colordict =  {c: colors[i] for i, c in enumerate(cats)}
        # plot the histogram_plot
        ax = gs.histogram_plot(dfl, cat=True, color=colordict, figsize=(7, 4), rotateticks=(45, 0),
                        label_count=True)

    """
    import pygeostat as gs
    from . utils import format_plot, _set_stat_fontsize, _format_grid, _format_tick_labels, setup_plot, catcmapfromcontinuous
    from . cmaps import _cat_pastel_data, _cat_vibrant_data
    import copy
    # Now converting to a numpy array, as encountering some odd pandas performance, and there's
    # no major disadvantagve to application of a numpy in this context to my knowledge - RMB
    # If a list is passed convert it to a series so that trimming can take place
    # weights
    if isinstance(weights, str):
        if isinstance(data, pd.DataFrame) or isinstance(data, gs.DataFile):
            weights = data[weights]
    elif isinstance(weights, bool):
        if weights:
            if isinstance(data, gs.DataFile):
                if data.weights is None:
                    raise ValueError('weights=True but data.weights is None!')
                elif isinstance(data.weights, list):
                    raise ValueError('weights=True but data.weights is a list!')
                weights = data[data.weights].values
            else:
                raise ValueError('weights=True is only valid if data is a DataFile!')
        else:
            weights = None
    if isinstance(weights, pd.Series) or isinstance(weights, pd.DataFrame):
        weights = weights.values
    # cats for continuous histogram_plots
    if isinstance(cat, str):
        if isinstance(data, pd.DataFrame) or isinstance(data, gs.DataFile):
            cat = data[cat]
    elif isinstance(cat, bool):
        if cat:
            if isinstance(data, gs.DataFile):
                if data.cat is None:
                    raise ValueError('cat=True but data.cat is None!')
                cat = data[data.cat].values
                if catdict is None and data.catdict is None:
                    raise ValueError("pass a `catdict` when setting `cat`")
                else:
                    catdict = data.catdict
            else:
                raise ValueError('cat=True is only valid if data is a DataFile!')
        else:
            cat = None
    if isinstance(cat, pd.Series) or isinstance(cat, pd.DataFrame):
        cat = cat.values
    # Handle categorical dictionary
    if isinstance(catdict, bool):
        if catdict:
            if not isinstance(data, gs.DataFile):
                raise ValueError('catdict as a bool is only valid if data is a DataFile!')
            if data.catdict is None:
                raise ValueError('catdict as a bool is only valid if data is not None!')
            catdict = data.catdict
    # Variable
    # Handle data that is 2-D and/or a DataFile
    if isinstance(var, str):
        if isinstance(data, pd.DataFrame) or isinstance(data, gs.DataFile):
            if isinstance(cat, str):
                cat = data[cat]
            data = data[var]
        else:
            raise ValueError('var as a string is only valid if data is a DataFile or DataFrame!')
    elif isinstance(data, gs.DataFile):
        if isinstance(data.variables, str):
            data = data[data.variables]
        elif cat is not None:
            if isinstance(cat, str):
                data = data[cat]
            elif var is None and isinstance(cat, (np.ndarray, list)):
                data = cat
        elif len(data.columns) == 1:
            data = data.data
        else:
            raise ValueError('Could not coerce data (DataFile) into a 1D dataset!')
    # Get the xlabel if possible before converting to a numpy array
    if isinstance(data, pd.Series) or isinstance(data, pd.DataFrame):
        if xlabel is None:
            xlabel = gs.get_label(data)
        data = data.values
    elif isinstance(data, list):
        data = np.array(data)
    if isinstance(cat, (pd.Series, pd.DataFrame)):
        cat = cat.values
    # Should be numpy by now...
    if data.ndim > 1:
        if data.shape[1] > 1:
            raise ValueError('Could not coerce data into a 1D dataset!')
        else:
            data = data.flatten()
    # Handle Null values if needed
    idx = np.isnan(data)
    nullcnt = np.sum(idx)
    if nullcnt > 0:
        data = data[~idx]
        if weights is not None:
            weights = weights[~idx]
        if cat is not None:
            cat = cat[~idx]
    # Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()
    # Set-up plot if no axis is supplied
    _, ax, _ = setup_plot(ax, figsize=figsize, aspect=False)
    # Infer some default parameters
    if weights is None:
        weights = np.ones(len(data)) / len(data)
    else:
        weights = weights / np.sum(weights)
    # Some quick error checks
    assert(np.all(weights) >= 0.0), 'weights less than 0 not valid'
    # Categories
    if isinstance(catdict, dict) and var is None:
        if not all([isinstance(float(i), float) for i in catdict.keys()]):
            raise ValueError('if catdict is dict., all keys should be an int/float!')
        # The bins are set to begin at the start of each category
        # bins go from 0.5 to (icat + 1) + 0.5
        # label is centered at (icat + 1)
        bins = np.arange(len(catdict) + 1) + 0.5
    if color is None and isinstance(catdict, dict):
        # Color each bin by the category color?
        if isinstance(Parameters['plotting.cmap_cat'], dict):
            temp = Parameters['plotting.cmap_cat']
            if list(sorted(temp.keys())) == list(sorted(catdict.keys())):
                color = temp
        else:
            color = catcmapfromcontinuous(Parameters["plotting.cmap"], len(catdict)).colors
    if isinstance(color, dict):
        if list(sorted(color.keys())) != list(sorted(catdict.keys())):
            raise ValueError(('if color is a dictionary, keys must align with '
                              'bins[:-1]! Consider using a single color.'))
        temp = color
        color = []
        for _, v in sorted(temp.items()):
            color.append(v)
    # Color setup
    if isinstance(color, int):
        # Grab a color from ``cat_vibrant`` if an integer is passed
        color = _cat_pastel_data[color % len(_cat_vibrant_data)]
    if not icdf:
        if color is None:
            color = Parameters['plotting.histogram_plot.facecolor']
        if edgecolor is None:
            edgecolor = Parameters['plotting.histogram_plot.edgecolor']
        if edgeweights is None:
            if "lw" in kwargs:
                edgeweights = kwargs.pop("lw")
            else:
                edgeweights = Parameters["plotting.histogram_plot.edgeweight"]
    else:
        if color is None and icdf:
            color = Parameters['plotting.histogram_plot.cdfcolor']
    plotdata = copy.deepcopy(data)
    plotweights = copy.deepcopy(weights)
    if xlim is not None:
        plotdata[data < xlim[0]] = xlim[0]
        plotdata[data > xlim[1]] = xlim[1]
    # Main plotting
    if icdf:

        def singlecdf(ax, data, weights, lower, upper, bins, color, label=None, **kwargs):
            """ local function to plot a single cdf """
            cdf_x, cdfvals = gs.cdf(data, weights=weights, lower=lower, upper=upper, bins=bins)
            # Matplotlib is a memory hog if to many points are used. Limit the number of points the CDF
            # is build with to 1000. The tails are given extra attention to make sure they are defined
            # nicely.
            if len(cdf_x) > 1000:
                cdfinterp = scipy.interpolate.interp1d(x=cdfvals, y=cdf_x)
                cdfvals = np.concatenate([np.arange(cdfvals.min(), 0.1, 0.001),
                                          np.arange(0.1, 0.9, 0.01),
                                          np.arange(0.9, cdfvals.max(), 0.001)])
                cdf_x = []
                for val in cdfvals:
                    cdf_x.append(cdfinterp(val))
                cdf_x = np.array(cdf_x)
            fig = ax.plot(cdf_x, cdfvals, color=color, label=label, **kwargs)
            return fig

        if catdict is not None:
            if var is not None:
                stat_blk = False
                for icat, c in enumerate(catdict):
                    clr = color[icat]
                    catidx = cat == c
                    fig = singlecdf(ax, plotdata[catidx], plotweights[catidx], lower, upper, bins, clr,
                                    label=catdict[c], **kwargs)
            else:
                raise ValueError("`icdf=True` and `catdict` only makes sense with a `var` defined")
        else:
            fig = singlecdf(ax, plotdata, plotweights, lower, upper, bins, color, **kwargs)
        if ylim is None:
            ylim = (0, 1.0)
    else:
        if bins is None:
            bins = Parameters['plotting.histogram_plot.histbins']
        label = kwargs.pop("label", None)
        if bins is None:
            if len(plotdata) < 200:
                bins = 20
            elif len(plotdata) < 500:
                bins = 25
            else:
                bins = 30
        if logx:
            if catdict is not None:
                raise ValueError('Cannot have logx with catdict!')
            if xlim is None:
                minv = np.log10(max(plotdata.min(), 1e-10))
                maxv = np.log10(plotdata.max())
            else:
                minv = np.log10(max(xlim[0], 1e-10))
                maxv = np.log10(xlim[1])
            if np.isnan([minv, maxv]).any():
                raise ValueError('ERROR converting your data to log base! are there negatives?')
            bins = np.logspace(minv, maxv, bins)
        if catdict is not None:
            if var is None:
                for icat, cat in enumerate(catdict):
                    plotdata[data == cat] = icat + 1
                histclr = None
            else:
                # generate lists of data per cat
                plotdata = [plotdata[cat == c] for c in catdict]
                plotweights = [weights[cat == c] for c in catdict]
                label = list(catdict.values())
                histtype = kwargs.pop("histtype", "stepfilled")
                stat_blk = False
                if "stacked" not in kwargs:
                    kwargs["stacked"] = True
                histclr = color
        histtype = kwargs.pop("histtype", "bar")
        if not isinstance(color, list):
            ax.hist(plotdata, bins, weights=plotweights, color=color, edgecolor=edgecolor,
                    histtype=histtype, label=label, lw=edgeweights, **kwargs)
        else:
            _, _, patches = ax.hist(plotdata, bins, weights=plotweights, histtype=histtype,
                                    color=histclr, edgecolor=edgecolor,
                                    label=label, lw=edgeweights, **kwargs)
            try:
                for patch, clr in zip(patches, color):
                    patch.set_facecolor(clr)
            except (AttributeError, ValueError):
                pass
        if catdict is not None and label_count:
            nd = len(data)
            for icat, cat in enumerate(catdict):
                count = np.count_nonzero(data == cat)
                pcat = (weights * (data == cat).astype(float)).sum()
                ax.text(icat + 1, pcat, count, ha="center", va="bottom")
    # Summary stats
    if stat_blk is None:
        stat_blk = Parameters['plotting.histogram_plot.stat_blk']
    if stat_xy is None:
        if icdf:
            stat_xy = Parameters['plotting.histogram_plot.stat_xy_cdf']
        else:
            stat_xy = Parameters['plotting.histogram_plot.stat_xy']
    if stat_blk:
        if sigfigs is None:
            sigfigs = Parameters['plotting.sigfigs']
        if roundstats is None:
            roundstats = Parameters['plotting.roundstats']
        if stat_ha is None:
            stat_ha = Parameters['plotting.stat_ha']
        if stat_linespacing is None:
            stat_linespacing = Parameters['plotting.stat_linespacing']
        if stat_linespacing is None:
            stat_linespacing = 1.0
        # Force no bins and upper/lower for median
        cdf_x, cdfvals = gs.cdf(data, weights=weights)
        # Currently defined statistics, possible to add more quite simply
        if np.mean(data) == 0:
            cdata = float("nan")
        elif roundstats:
            cdata = round((np.std(data) / np.mean(data)), sigfigs)
        else:
            cdata = gs.round_sigfig((np.std(data) / np.mean(data)), sigfigs)
        if roundstats:
            mean = round(gs.weighted_mean(data, weights), sigfigs)
            median = round(gs.percentile_from_cdf(cdf_x, cdfvals, 50.0), sigfigs)
            stdev = round(np.sqrt(gs.weighted_variance(data, weights)), sigfigs)
            minval = round(np.min(data), sigfigs)
            maxval = round(np.max(data), sigfigs)
            upquart = round(np.percentile(data, 75), sigfigs)
            lowquart = round(np.percentile(data, 25), sigfigs)
            p10 = round(np.percentile(data, 10), sigfigs)
            p90 = round(np.percentile(data, 90), sigfigs)
        else:
            mean = gs.round_sigfig(gs.weighted_mean(data, weights), sigfigs)
            median = gs.round_sigfig(gs.percentile_from_cdf(cdf_x, cdfvals, 50.0), sigfigs)
            stdev = gs.round_sigfig(np.sqrt(gs.weighted_variance(data, weights)), sigfigs)
            minval = gs.round_sigfig(np.min(data), sigfigs)
            maxval = gs.round_sigfig(np.max(data), sigfigs)
            upquart = gs.round_sigfig(np.percentile(data, 75), sigfigs)
            lowquart = gs.round_sigfig(np.percentile(data, 25), sigfigs)
            p10 = gs.round_sigfig(np.percentile(data, 10), sigfigs)
            p90 = gs.round_sigfig(np.percentile(data, 90), sigfigs)
        statistics = {'mean': (r'$m = %g$' % mean),
                      'median': (r'$x_{{50}} = %g$' % median),
                      'count': ('$n = %i$' % len(data)),
                      'count_trimmed': ('$n_{trim} = %i$' % nullcnt),
                      'stdev': (r'$\sigma = %g$' % stdev),
                      'cvar': ('$CV = %g$' % cdata),
                      'min': ('$x_{{min}} = %g$' % minval),
                      'max': ('$x_{{max}} = %g$' % maxval),
                      'upquart': ('$x_{{75}} = %g$' % upquart),
                      'lowquart': ('$x_{{25}} = %g$' % lowquart),
                      'p10': ('$x_{{10}} = %g$' % p10),
                      'p90': ('$x_{{90}} = %g$' % p90)}
        # Default statistic sets
        if stat_blk == 'varlabel' and 'label' in kwargs:
            statistics['varlabel'] = kwargs['label']
        statsets = {'minimal': ['count', 'mean', 'median', 'stdev'],
                    'all': ['count', 'mean', 'stdev', 'cvar', 'max', 'upquart', 'median',
                            'lowquart', 'min'],
                    'varlabel': ['varlabel', 'count', 'mean', 'stdev', 'cvar', 'max', 'upquart',
                                 'median', 'lowquart', 'min'],
                    'none': None}
        # Use a default statistic set
        if isinstance(stat_blk, bool) and stat_blk:
            stat_blk = 'all'
        if isinstance(stat_blk, str):
            if stat_blk in statsets:
                stat_blk = statsets[stat_blk]
            else:
                print('WARNING: stats value of: "' + stat_blk + '" does not exist - '
                      'default to no stats')
                stat_blk = None
        # Use a supplied statistic set, but check for bad ones
        else:
            badstats = [s for s in stat_blk if s not in statistics]
            stat_blk = [s for s in stat_blk if s in statistics]
            for badstat in badstats:
                print('WARNING: stats value of: "' + badstat + '" does not exist - '
                      'It was removed from summary statistics list')
        # Form the stats string
        if stat_blk:
            if nullcnt != 0:
                stat_blk.insert(stat_blk.index('count') + 1, 'count_trimmed')
            stat_blk = [statistics[s] for s in stat_blk]
            txtstats = '\n'.join(stat_blk)
            if len(np.unique(weights)) > 1:
                txtstats = txtstats + '\n\nweights used'
            if stat_xy[1] > 0.5:
                va = 'top'
            else:
                va = 'bottom'
            # Set the stat_fontsize
            stat_fontsize = _set_stat_fontsize(stat_fontsize)
            ax.text(stat_xy[0], stat_xy[1], txtstats, va=va, ha=stat_ha, transform=ax.transAxes,
                    fontsize=stat_fontsize, linespacing=stat_linespacing)
    # Label as required
    if icdf:
        ylabel = 'Cumulative Distribution Function'
    elif 'density' in kwargs:
        ylabel = 'Probability Density Function (PDF)'
    else:
        ylabel = 'Frequency'
    ax = format_plot(ax, xlabel, ylabel, title, axis_xy=axis_xy, xlim=xlim, ylim=ylim, logx=logx)
    if catdict is not None and var is None:
        ticlocs = [i + 1 for i in range(len(catdict.keys()))]
        ax.set_xticks(ticlocs)
        ax.set_xticklabels(catdict.values())
        ax.set_xlim(0.25, len(catdict) + 0.75)
    elif catdict is not None and var is not None:
        ax.legend()
    _format_tick_labels(ax, rotateticks)
    # format_plot doesn't handle some specialized axis_xy and grid requirements
    # for histogram_plot...
    if icdf:
        # Ensure that we have top spline, in case it was removed above
        ax.spines['top'].set_visible(True)
        _format_grid(ax, grid, below=False)
    else:
        # The grid should be below for a histogram
        _format_grid(ax, grid, below=True)
    # Export figure
    if output_file or ('pdfpages' in out_kws):
        gs.export_image(output_file, **out_kws)
    return ax
