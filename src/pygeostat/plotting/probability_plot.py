#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A basic probability plotting routine reminiscent of probability_plot from gslib"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from . set_style import set_plot_style


@set_plot_style
def probability_plot(data, wt=None, lower=None, upper=None, logscale=True, ax=None, figsize=None,
                    xlim=None, ylim=(0.01, 99.99), nyticks=15, title=None, xlabel=None, s=2, color='k',
                    grid=True, axis_xy=None, line=False, plot_style=None, custom_style=None, output_file=None,
                    out_kws=None, **kwargs):
    """
    Create either a normal or a lognormal probability plot. This plot displays all the data
    values on a chart that illustrates the general distribution shape and the behavior of the
    extreme values. Please review the documentation of the :func:`gs.set_style()
    <pygeostat.plotting.set_style.set_style>` and :func:`gs.export_image()
    <pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
    their use in this function can be understood.

    This function requires the python package probscale. It can be installed be executing the
    following code in the command prompt:

        >>> pip install probscale

    Parameters:
        data: Tidy (long-form) 1D data where a single column of the variable exists with each row
            is an observation. A pandas dataframe/series or numpy array can be passed.
        wt: 1D dataframe, series, or numpy array of declustering weights for the data.
        lower (float): Lower trimming limits
        upper (float): Upper trimming limits
        ax (mpl.axis): Matplotlib axis to plot the figure
        figsize (tuple): Figure size (width, height)
        xlim (float tuple): Minimum and maximum limits of data along the x axis
        ylim (float tuple): Minimum and maximum limits of data along the y axis, e.g.(0.001, 99.999)
        nyticks (int): the number of ticks on the y axis to show. Currently disabled due to altered
            matplotlib functionality.
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
        **kwargs: Optional permissible keyword arguments to pass to matplotlib's scatter function

    Returns:
        ax (ax): matplotlib Axes object with the histogram

    **Examples:**

        A simple call:

        .. plot::

            import pygeostat as gs
            data = gs.ExampleData('oilsands')
            gs.probability_plot(data.data['Fines'], logscale=False)

    """
    # Import the rest of the packages
    import numpy as np
    import pandas as pd
    from .export_image import export_image
    from .. statistics.cdf import cdf
    import matplotlib.pyplot as plt
    from .utils import format_plot, get_label
    try:
        import probscale
    except ImportError:
        raise ImportError('Looks like probscale is not installed, use:\n' +
                          '>>> pip install probscale\n')
    # Ensure valid VTK option
    # Coerce the passed data and wt into a single pandas dataframe
    if xlabel is None:
        xlabel = get_label(data)
    try:
        data = pd.DataFrame(data=data)
        data.reset_index(inplace=True, drop=True)
        if data.shape[1] != 1:
            raise ValueError("The passed `data` is not 1-D")
    except:
        raise ValueError("Please ensure the passed `data` can be coerced into a pandas dataframe"
                         " (i.e., data = pd.DataFrame(data=data)")
    if wt is not None:
        try:
            wt = pd.DataFrame(data=wt)
            if wt.shape[1] != 1:
                raise ValueError("The passed `wt` is not 1-D")
        except:
            raise ValueError("Please ensure the passed `wt` can be coerced into a pandas"
                             " dataframe (i.e., wt = pd.DataFrame(wt=wt)")
    else:
        wt = pd.DataFrame(data=np.ones((len(data), 1)))
    data = pd.concat([data, wt], axis=1)
    data.columns = ['data', 'wt']
    # Trim the data
    if lower is None:
        lower = data['data'].min()
    if upper is None:
        upper = data['data'].max()
    data.query("%s <= data <= %s" % (lower, upper))
    # Sanity checks
    assert((data['wt']).all() >= 0.0), 'Weights less than 0 not valid'
    # Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()
    # Calculate the cdf
    cdf_x, cdfvals = cdf(data['data'], weights=data['wt'])
    cdfvals = cdfvals * 100
    # Set-up plot if no axis is supplied
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    # Plot the figure
    if line:
        ax.plot(cdf_x, cdfvals, color=color, **kwargs)
    else:
        ax.scatter(cdf_x, cdfvals, s=s, c=color, lw=0, **kwargs)
    # Set some figure aesthetics
    if logscale:
        ax.set_xscale('log')
    ax.set_yscale('prob')
    ax.locator_params(axis='y')
    ax = format_plot(ax, xlabel=xlabel, ylabel='Cumulative Probability', title=title,
                   grid=grid, axis_xy=axis_xy, xlim=xlim, ylim=ylim)
    # Export figure
    if output_file or ('pdfpages' in out_kws):
        export_image(output_file, **out_kws)

    return ax
