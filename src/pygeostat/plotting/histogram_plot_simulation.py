#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import glob
import os
import types
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters

@set_plot_style
def histogram_plot_simulation(simulated_data, reference_data, reference_variable=None, reference_weight=None, reference_n_sample=None,
                            simulated_column=None, griddef=None, nreal=None,
                            n_subsample=None, simulated_limits=False, ax=None,
                            figsize=None, xlim=None, title=None, xlabel=None, stat_blk='all',
                            stat_xy=(0.95, 0.05), reference_color=None, simulation_color=None, alpha=None, lw=1,
                            plot_style=None, custom_style=None, output_file=None, out_kws=None, sim_kws=None,
                            **kwargs):
    """
    histogram_plot_simulation emulates the pygeostat histogram_plot program as a means of checking histogram
    reproduction of simulated realizations to the original histogram. The use of python generators
    is a very flexible and easy means of instructing this plotting function as to what to plot.

    The function accepts five types of simulated input passed to the ``simulated_data`` argument:

        #. 1-D array like data (numpy or pandas) containing 1 or more realizations of simulated
           data.
        #. 2-D array like data (numpy or pandas) with each column being a realization and each row
           being an observation.
        #. List containing location(s) of realization file(s).
        #. String containing the location of a folder containing realization files. All files in
           the folder are read in this case.Can contain
        #. String with a wild card search (e.g., './data/realizations/*.out')
        #. Python generator object that yields a 1-D numpy array.

    The function accepts two types of reference input passed to the ``reference_data`` argument:

        #. Array like data containing the reference variable
        #. String containing the location of the reference data file (e.g., './data/data.out')

    This function uses pygeostat for plotting and numpy to calculate statistics.

    The only parameters required are ``reference_data`` and ``simulated_data``. If files are to be read or a 1-D
    array is passed, the parameters ``griddef`` and ``nreal`` are required. ``simulated_column`` is required
    for reading files as well. It is assumed that an equal number of realizations are within each
    file if multiple file locations are passed. Sub-sampling of datafiles can be completed by
    passing the parameter ``n_subsample``. If a file location is passed to ``reference_data``, the parameters
    ``reference_variable`` and ``reference_n_sample`` are required. All other arguments are optional or determined
    automatically if left at their default values. If ``xlabel`` is left to its default value of
    ``None``, the column information will be used to label the axes if present. Three keyword
    dictionaries can be defined. (1) ``sim_kws`` will be passed to pygeostat histogram_plot used for
    plotting realizations (2) ``out_kws`` will be passed to the pygeostat exportfig function and
    (3) ``**kwargs`` will be passed to the pygeostat histogram_plot used to plot the reference data.


    Two statistics block sets are available: ``'minimal'`` and the default ``'all'``. The
    statistics block can be customized to a user defined list and order. Available statistics are
    as follows:

    >>> ['nreal', 'realavg', 'realavgstd', 'realstd', 'realstdstd', 'ndat', 'refavg', 'refstd']

    Please review the documentation of the :func:`gs.set_plot_style() <pygeostat.plotting.set_plot_style>` and
    :func:`gs.export_image() <pygeostat.plotting.export_image>` functions for details on their
    parameters so that their use in this function can be understood.

    Parameters:
        simulated_data: Input simulation data
        reference_data: Input reference data

    Keyword Arguments:
        reference_variable (int, str): Required if sub-sampling reference data. The column containing the data
            to be sub-sampled
        reference_weight: 1D dataframe, series, or numpy array of declustering weights for the data. Can also
            be a string of the column in the reference_data if reference_data is a string, or a bool if reference_data.weights
            is a string
        reference_n_sample (int): Required if sub-sampling reference data. The number of data within the
            reference data file to sample from
        griddef (GridDef): A pygeostat class GridDef created using :class:`gs.GridDef
            <pygeostat.data.grid_definition.GridDef>`
        simulated_column (int): column number in the simulated data file
        nreal (int): Required if sub-sampling simulation data. The total number of realizations
            that are being plotted. If a HDF5 file is passed, this parameter can be used to limit
            the amount of realizations plotted (i.e., the first ``nreal`` realizations)
        n_subsample (int): Required if sub-sampling is used. The number of sub-samples to draw.
        ax (mpl.axis): Matplotlib axis to plot the figure
        figsize (tuple): Figure size (width, height)
        xlim (float tuple): Minimum and maximum limits of data along the x axis
        title (str): Title for the plot
        xlabel (str): X-axis label
        stat_blk (str or list): Indicate what preset statistics block to write or a specific list
        stat_xy (str or float tuple): X, Y coordinates of the annotated statistics in figure
            space. The default coordinates specify the bottom right corner of the text block
        reference_color (str): Colour of original histogram
        simulation_color (str): Colour of simulation histograms
        alpha (float): Transparency for realization variograms (0 = Transparent, 1 = Opaque)
        lw (float): Line width in points. The width provided in this parameter is used for the
            reference variogram, half of the value is used for the realization variograms.
        plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
            to turn it off
        custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
        sim_kws: Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.histogram_plot() <pygeostat.plotting.histogram_plot.histogram_plot>` for plotting realization
            histograms and by extension, matplotlib's plot function if the keyword passed is not
            used by :func:`gs.histogram_plot() <pygeostat.plotting.histogram_plot.histogram_plot>`
        **kwargs: Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.histogram_plot() <pygeostat.plotting.histogram_plot.histogram_plot>` for plotting the reference
            histogram and by extension, matplotlib's plot function if the keyword passed is not
            used by :func:`gs.histogram_plot() <pygeostat.plotting.histogram_plot.histogram_plot>`

    Returns:
        ax (ax): matplotlib Axes object with the histogram reproduction plot

    **Examples:**

    .. plot:: 
    
        import pygeostat as gs
        import pandas as pd

        # Global setting using Parameters
        gs.Parameters['data.griddef'] = gs.GridDef([10,1,0.5, 10,1,0.5, 2,1,0.5])
        gs.Parameters['data.nreal'] = nreal = 100
        size = gs.Parameters['data.griddef'].count();

        reference_data = pd.DataFrame({'value':np.random.normal(0, 1, size = size)})

        # Create example simulated data
        simulated_data =pd.DataFrame(columns=['value'])
        for i, offset in enumerate(np.random.rand(nreal)*0.04 - 0.02):
            simulated_data= simulated_data.append(pd.DataFrame({'value':np.random.normal(offset, 1, size = size)}))

        gs.histogram_plot_simulation(simulated_data, reference_data, reference_variable='value',
                                title='Histogram Reproduction', grid=True)
    """
    # -----------------------------------------------------------------------
    #  Sanity checks, file type determination, and try loading fortran
    # -----------------------------------------------------------------------
    import pygeostat as gs
    from . utils import format_plot, _set_stat_fontsize
    # Figure out what type of input simulated_data is
    subsamp = False
    ndim = False
    generator = False
    folder = False
    wildcard = False
    array = False
    filelist = False
    if isinstance(simulated_data, types.GeneratorType):
        nreal = 0
        generator = True
        iterator = simulated_data
    elif isinstance(simulated_data, str) and ('*' in simulated_data):
        wildcard = True
    elif isinstance(simulated_data, str) and os.path.isdir(simulated_data):
        folder = True
    elif isinstance(simulated_data, list):
        filelist = True
    elif any([isinstance(simulated_data, pd.DataFrame), isinstance(simulated_data, np.ndarray), isinstance(simulated_data, pd.Series)]):
        array = True
        if subsamp:
            raise ValueError("Sub-sampling won't work if the data is already in memory")
    else:
        raise ValueError("The passed `simulated_data` is not a valid input format")
    if n_subsample is not None and isinstance(n_subsample, (int, float)):
        subsamp = True
    # Make sure the required parameters are passed
    if griddef is None:
        griddef = Parameters['data.griddef']
        if griddef is None:
            raise ValueError("A gs.GridDef must be passed when reading from files")

    if nreal is None:
            nreal = Parameters['data.nreal']
            if nreal is None:
                raise ValueError("The number of realizations to be read must be specified when"
                                 " reading from files")

    if any([folder, wildcard, filelist]):
        
        if simulated_column is None:
            raise ValueError("The column in the files that contains the simulation data must be"
                             " specified")
    elif array:
        if (len(simulated_data.shape) == 1) or (simulated_data.shape[1] == 1):
            if not isinstance(griddef, gs.GridDef):
                raise ValueError("If a 1-D array is passed, a gs.GridDef must be passed to"
                                 " `griddef`")
            if nreal is None:
                raise ValueError("The number of realizations must be passed if dealing with a"
                                 " 1-D array")
    # Figure out what type of input reference_data is
    if isinstance(reference_data, str) and (n_subsample is not False):
        refsubsamp = True
    else:
        refsubsamp = False
    # Try to load the subsample function if it required
    if subsamp or refsubsamp:
        if not isinstance(griddef, gs.GridDef):
            raise ValueError("A gs.GridDef is required for subsampling.")
        try:
            from pygeostat.fortran.subsample import subsample
        except:
            raise ImportError("The fortran subroutine subsample could not be loaded, please ensure"
                              " it has been compiled correctly.")

    # -----------------------------------------------------------------------
    #  Handle data input
    # -----------------------------------------------------------------------
    # Set-up variables
    realavg, realavgstd, realstd, realstdstd, refavg, refstd = ([] for i in range(6))
    # Handle pd and np input
    if array:
        if isinstance(simulated_data, pd.DataFrame):
            simulated_data = simulated_data.values
        if (len(simulated_data.shape) == 2) and (simulated_data.shape[1] > 1):
            ndim = 2
        else:
            ndim = 1
        # Handle 1-D arrays
        if ndim == 1:
            ncell = griddef.count()
        # Handle 2-D arrays
        if ndim == 2:
            nreal = simulated_data.shape[1]
    # Handle folder and wildcard searches
    if folder:
        if simulated_data[-1] != '/':
            simulated_data = simulated_data + '/'
        simulated_data = simulated_data + '*'
        files = []
        for filepath in glob.glob(simulated_data):
            files.append(filepath)
    if wildcard:
        files = []
        for filepath in glob.glob(simulated_data):
            files.append(filepath)
    if filelist:
        files = simulated_data

    if any([folder, wildcard, filelist]):
        ncell = griddef.count()
        ndim = 1
        # Check and make sure the number of files and the nreal value passed makes sense
        if nreal % len(files) != 0:
            raise ValueError(" The number of realizations passed is not divisible by the number"
                             " of files passed/found. Please make sure there are the same number"
                             " of realizations in each file and that the sum of them match the"
                             " nreal argument passed.")
        # Read the data
        simulated_data = []
        for file in files:
            data = gs.DataFile(file).data.iloc[:, simulated_column - 1]
            if simulated_limits:
                if isinstance(simulated_limits, (int, float)):
                    data = data.loc[data > simulated_limits]
                elif len(simulated_limits) == 2:
                    data = data.loc[data.between(simulated_limits[0], simulated_limits[1])]
                else:
                    raise ValueError('simulated_limits must be a value or tuple')
            simulated_data.extend(data)
        simulated_data = np.array(simulated_data)
        if simulated_limits:
            ncell = len(data)
    # Handle sub-sampling if required
    if subsamp:
        ncell = griddef.count()
        if array:
            if ndim == 1:
                simulated_data = np.reshape(simulated_data, (nreal, griddef.count())).T
                ndim = 2
            newarr = np.zeros((n_subsample, nreal))
            for ireal in range(nreal):
                ridx = np.random.permutation(ncell)[:n_subsample]
                newarr[:, ireal] = simulated_data[ridx, ireal]
            simulated_data = newarr
        else:
            # Sub-sample all of the files and combine into a single numpy array
            files = simulated_data
            file_nreal = int(nreal / len(files))
            simulated_data = []
            for fl in files:
                dump = subsample(fl, simulated_column, ncell, n_subsample, file_nreal, gs.rseed())
                dump = np.transpose(dump)
                simulated_data.extend(dump)
            simulated_data = np.array(simulated_data)
            simulated_data = np.transpose(simulated_data)
            nreal = simulated_data.shape[1]

    # Create a generator for the realizations
    if generator is False:
        def _itersimulated_data():
            for i in range(0, nreal):
                if subsamp or ndim == 2:
                    real = simulated_data[:, i]
                elif ndim == 1:
                    real = simulated_data[(i * ncell):(((i + 1) * ncell) - 1)]
                yield real
        iterator = _itersimulated_data()

    # -----------------------------------------------------------------------
    #  Plot Figure
    # -----------------------------------------------------------------------
    # Set figure style parameters

    # Handle dictionary defaults
    if sim_kws is None:
        sim_kws = dict()
    if out_kws is None:
        out_kws = dict()
    # Set-up figure
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    if simulation_color is None:
        simulation_color = Parameters['plotting.histogram_plot_simulation.simclr']
    if alpha is None:
        alpha = Parameters['plotting.histogram_plot_simulation.alpha']
    # Plot realization histograms
    for real in iterator:
        # Append intermediate realization dist statistics to variables
        if stat_blk:
            realavg.append(np.nanmean(real))
            realstd.append(np.nanstd(real))
        if generator:
            nreal += 1
        gs.histogram_plot(real, ax=ax, icdf=True, lw=(lw / 2), stat_blk=False, color=simulation_color, alpha=alpha, plot_style=False, **sim_kws)
    # Calculate more realization dist statistics if required
    if stat_blk:
        realavgstd = np.std(realavg)
        realavg = np.mean(realavg)
        realstdstd = np.std(realstd)
        realstd = np.mean(realstd)
    # Sub-sample original distribution if needed

    if isinstance(reference_data, str):
        if subsamp:
            reference_data = subsample(reference_data, reference_variable, reference_n_sample, n_subsample, 1, gs.rseed())
            reference_data = reference_data[:, 0]
        else:
            reference_data = gs.DataFile(reference_data)
            if isinstance(reference_variable, str):
                reference_variable = reference_data.gscol(reference_variable) - 1
            if isinstance(reference_weight, str):
                reference_weight = reference_data.gscol(reference_weight) - 1
                reference_weight = reference_data.data.values[:, reference_weight]
            elif isinstance(reference_weight, bool):
                if reference_weight:
                    if isinstance(reference_data.weights, str):
                        reference_weight = reference_data[reference_data.weights]
                    else:
                        raise ValueError('reference_weight=True is only valid if reference_data.weights is a string!')
            reference_data = reference_data.data.values[:, reference_variable]
    elif isinstance(reference_data, gs.DataFile):
        if isinstance(reference_weight, str):
            reference_weight = reference_data[reference_weight].values
        elif isinstance(reference_weight, bool):
            if reference_weight:
                if isinstance(reference_data.weights, str):
                    reference_weight = reference_data[reference_data.weights].values
                else:
                    raise ValueError('reference_weight=True is only valid if reference_data.weights is a string!')
        if isinstance(reference_variable, str):
            reference_data = reference_data[reference_variable].values
        elif isinstance(reference_data.variables, str):
            reference_data = reference_data[reference_data.variables].values
        elif len(list(reference_data.columns)) == 1:
            reference_data = reference_data.data.values
        else:
            raise ValueError('could not coerce reference_data into a 1D array!')

    # Plot the reference histogram
    if reference_color is None:
        reference_color = Parameters['plotting.histogram_plot_simulation.refclr']
    if not isinstance(reference_data, bool) and (reference_data is not None):
        gs.histogram_plot(reference_data, weights=reference_weight, ax=ax, icdf=True, stat_blk=False, color=reference_color,
                   lw=lw, plot_style=False, **kwargs)
    # Calculate reference dist statistics if required
    if stat_blk:
        if not isinstance(reference_data, bool) and (reference_data is not None):
            refavg = gs.weighted_mean(reference_data, reference_weight)
            refstd = np.sqrt(gs.weighted_variance(reference_data, reference_weight))
            ndat = len(reference_data)
        else:
            refavg = np.nan
            refstd = np.nan
            ndat = np.nan
    # Configure plot
    if xlabel is None:
        xlabel = gs.get_label(reference_data)
        if xlabel is None:
            xlabel = gs.get_label(simulated_data)
    # axis_xy and grid are applied by format_plot based on the current Parameters setting
    # no kwarg in this function for now since it's already loaded
    ax = format_plot(ax, xlabel, 'Cumulative Frequency', title, xlim=xlim)
    # Ensure that we have top spline, in case it was removed above
    ax.spines['top'].set_visible(True)
    # Plot statistics block
    if stat_blk:
        statlist = {'nreal': (r'$n_{real} = %0.0f$' % nreal),
                    'realavg': (r'$m_{real} = %0.3f$' % realavg),
                    'realavgstd': (r'$\sigma_{m_{real}} = %0.3f$' % realavgstd),
                    'realstd': (r'$\sigma_{real} = %0.3f$' % realstd),
                    'realstdstd': (r'$\sigma_{\sigma_{real}} = %0.3f$' % realstdstd),
                    'ndat': ('$n_{ref} = %0.0f$' % ndat),
                    'refavg': (r'$m_{ref} = %0.3f$' % refavg),
                    'refstd': (r'$\sigma_{ref} = %0.3f$' % refstd)}
        statsets = {'all': ['nreal', 'realavg', 'realavgstd', 'realstd', 'realstdstd', 'ndat',
                            'refavg', 'refstd'],
                    'minimal': ['nreal']}
        if subsamp:
            statlist['n_subsample'] = '$n_{subsample} = %0.0f$' % n_subsample
            statsets['all'].append('n_subsample')
        txtstats, stat_xy, ha, va = gs.get_statblk(stat_blk, statsets, statlist, stat_xy)
        stat_fontsize = _set_stat_fontsize(None)
        ax.text(stat_xy[0], stat_xy[1], txtstats, va=va, ha=ha, fontsize=stat_fontsize,
                transform=ax.transAxes)
    # Export figure
    if output_file or ('pdfpages' in out_kws):
        gs.export_image(output_file, **out_kws)

    return ax
