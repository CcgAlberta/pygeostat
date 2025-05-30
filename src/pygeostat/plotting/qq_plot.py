#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A basic quantile by quantile plot to compare two probability distribution"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


from . set_style import set_plot_style


@set_plot_style
def qq_plot(data, reference_data, data_weight=None, reference_weight=None, limits=None, npoints=0,
          log_scale=None, ax=None, figsize=None, title=None, xlabel=None,
          ylabel=None, s=None, percent=True, color='k', grid=None, axis_xy=None, plot_style=None,
          custom_style=None, output_file=None, out_kws=None, line=True, ntickbins=5, **kwargs):
    """
    Plot a QQ plot between the reference. Pretty much the probplt but with 2 datasets and plotting
    the quantiles between them

    Parameters:
        data: Tidy (long-form) 1D data where a single column of the variable exists with each row
            is an observation. A pandas dataframe/series or numpy array can be passed.
        reference_data: Tidy (long-form) 1D data or a valid scipy.stats distribution (e.g. "norm"). A
            pandas dataframe/series or numpy array can be passed.
        data_weight: 1D dataframe, series, or numpy array of declustering weights for the data.
        reference_weight: 1D dataframe, series, or numpy array of declustering weights for the data.
        lower (float): Lower trimming limits
        upper (float): Upper trimming limits
        limits (tuple): the min and max value of the axes
        ax (mpl.axis): Matplotlib axis to plot the figure
        log_scale (bool): yes or no to log_scale
        npoints (int): set to 0 to use all points
        figsize (tuple): Figure size (width, height)
        xlim (float tuple): Minimum and maximum limits of data along the x axis
        title (str): Title for the plot
        xlabel (str): X-axis label. A default value of ``None`` will try and grab a label from the
            passed ``data``. Pass ``False`` to not have an xlabel.
        s (int): Size of points
        color (str or int): Any permissible matplotlib color or a integer which is used to draw a
            color from the pygeostat color pallet ``pallet_pastel`` (useful for iteration)
        grid(bool): plots the major grid lines if True. Based on gsParams['plotting.grid']
            if None.
        axis_xy (bool): converts the axis to GSLIB-style axis visibility (only left and bottom
            visible) if axis_xy is True. Based on gsParams['plotting.axis_xy'] if None.
        plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
            to turn it off
        custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
        line (bool): Plot the reference 1:1 line
        ntickbins (int or tuple): modify the number of ticks. Only works if log_scale == ``True``
        **kwargs: Optional permissible keyword arguments to pass to matplotlib's scatter function

    Returns:
        ax (ax): matplotlib Axes object with the histogram

    **Examples:**

    A simple call:
    
    .. plot::

        import pygeostat as gs
        import numpy as np
        # load some data
        gs.qq_plot(np.random.randn(1000),np.random.randn(1000))

    """

    # Import the rest of the packages
    import numpy as np
    import pandas as pd
    from .utils import get_label, format_plot
    from .export_image import export_image
    from ..statistics.cdf import cdf
    import matplotlib.pyplot as plt

    if xlabel is None:
        xlabel = get_label(data)
    if ylabel is None:
        ylabel = get_label(reference_data)
    # Coerce the passed data and wt into a single pandas dataframe
    
    try:
        data = pd.DataFrame(data=data)
        data.reset_index(inplace=True, drop=True)
        if data.shape[1] != 1:
            raise ValueError("The passed `data` is not 1-D")
    except:
        raise ValueError("Please ensure the passed `data` can be coerced into a pandas dataframe"
                         " (i.e., data = pd.DataFrame(data=data)")
    if data_weight is not None:
        try:
            data_weight = pd.DataFrame(data=data_weight, columns=["wt"])
            if data_weight.shape[1] != 1:
                raise ValueError("The passed `data_weight` is not 1-D")
        except:
            raise ValueError("Please ensure the passed `data_weight` can be coerced into a pandas"
                             " dataframe (i.e., data_weight = pd.DataFrame(wt=data_weight)")
    else:
        data_weight = pd.DataFrame(data=np.ones((len(data), 1)))
    data = pd.concat([data, data_weight], axis=1)
    data.columns = ['data', 'data_weight']

    # Do the same for the reference data
    try:
        reference_data = pd.DataFrame(data=reference_data)
        reference_data.reset_index(inplace=True, drop=True)
        if reference_data.shape[1] != 1:
            raise ValueError("The passed `reference_data` is not 1-D")
    except:
        raise ValueError("Please ensure the passed `reference_data` can be coerced into a pandas dataframe"
                         " (i.e., reference_data = pd.DataFrame(data=reference_data)")
    if reference_weight is not None:
        try:
            reference_weight = pd.DataFrame(data=reference_weight)
            if reference_weight.shape[1] != 1:
                raise ValueError("The passed `reference_weight` is not 1-D")
        except:
            raise ValueError("Please ensure the passed `reference_weight` can be coerced into a pandas"
                             " dataframe (i.e., reference_weight = pd.DataFrame(reference_weight=reference_weight)")
    else:
        reference_weight = pd.DataFrame(data=np.ones((len(reference_data), 1)))
    reference_data = pd.concat([reference_data, reference_weight], axis=1)
    reference_data.columns = ['reference_data', 'reference_weight']

    # Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()

    # build the cdf's using all bins
    print(len(data))
    mp1, p1 = cdf(data['data'], weights=data['data_weight'])
    mp2, p2 = cdf(reference_data['reference_data'], weights=reference_data['reference_weight'])

    if npoints <= 0:
        npoints = min([len(data), len(reference_data), 1e10])

    # generate and interpolate the equally spaced set of quantiles
    probs = np.linspace(1 / npoints, 1, npoints)
    x = np.interp(probs, p2, mp2)
    y = np.interp(probs, p1, mp1)

    # Set-up plot if no axis is supplied
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    # Plot the figure
    ax.scatter(x, y, s=s, c=color, lw=0, **kwargs)

    if limits is None:
        limits = np.minimum(np.min(mp1), np.min(mp2)), np.maximum(np.max(mp1), np.max(mp2))

    # apply the log_scale
    if log_scale:
        ax.set_xscale('log')
        ax.set_yscale('log')

    if line:
        ax.plot(limits, limits, color='r', lw=0.5)

    # deal with the labels and formatting
    if xlabel is None:
        xlabel = '$ref_{quantiles}$'
    if ylabel is None:
        ylabel = '$data_{quantiles}$'
    ax = format_plot(ax, xlabel=xlabel, ylabel=ylabel, title=title, grid=grid, axis_xy=axis_xy,
                   xlim=limits, ylim=limits)
    # Export figure
    if output_file or ('pdfpages' in out_kws):
        export_image(output_file, **out_kws)

    return ax
