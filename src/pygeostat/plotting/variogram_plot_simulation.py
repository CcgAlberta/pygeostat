#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Varplot wrapper to allow variogram reproduction plotting"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import matplotlib as mpl
import matplotlib.pyplot as plt

from . set_style import set_plot_style
from .. pygeostat_parameters import Parameters

@set_plot_style
def variogram_plot_simulation(simulated_data, simulated_var_id, reference_data=False, reference_var_id=1, sill=1, ax=None, figsize=None,
                            xlim=None, ylim=None, trimylim=True, title=None, variable=None, xlabel=None,
                            unit=None, ylabel=None, refclr=None, simclr=None, alpha=None, ls='-', lw=1,
                            lw_real=None, legend_label=None, plot_style=None, custom_style=None, output_file=None,
                            out_kws=None, **kwargs):
    """
    variogram_plot_simulation provides a means of checking variogram reproduction of simulation realizations to the
    input variogram model. A variogram reference system is used to know which variograms to plot
    from a output GSLIB varsim file. Please refer to :func:`gs.get_uniquevarids()
    <pygeostat.plotting.variogram_utils.get_uniquevarids>` for information to understand how it operates.

    The function uses pygeostats :func:`gs.variogram_plot() <pygeostat.plotting.variogram_plot.variogram_plot>` to plot
    both simulation, and if present, reference variograms.

    The only required parameters are ``simulated_data`` and ``simulated_var_id``. The function only accepts
    ``simulated_data`` as a pandas dataframe. ``simulated_var_id`` is the corresponding  ID value derived from
    :func:`gs.get_uniquevarids() <pygeostat.plotting.variogram_utils.get_uniquevarids>`, to the function
    knows which  variograms to plot. All other parameters are optional or have default values. The
    line width of the simulation plots are hard coded to be half of the size input with the ``lw``
    parameter.

    Please review the documentation of the :func:`gs.set_style()
    <pygeostat.plotting.set_style.set_style>` and :func:`gs.export_image()
    <pygeostat.plotting.export_image.export_image>` functions for details on their parameters so that
    their use in this function can be understood.

    Parameters:
        simulated_data (pd.DataFrame): Dataframe containing the following columns: 'Variogram Index',
            'Variogram Number', 'Calculation Azimuth', and 'Calculation Dip'. Variogram distance
            data is also needed as must have one of the following column headers: 'h', 'Lag
            Distance', or 'Distance'
        simulated_var_id (int): A 'Variogram ID' derived from :func:`gs.get_uniquevarids()
            <pygeostat.plotting.variogram_utils.get_uniquevarids>`
        reference_data (pd.DataFrame or bool): Dataframe containing the reference variogram value,
            variogram distance, and variogram index (if required) data as columns. The dataframe
            must contain the correct column IDs. The column header containing the variogram
            distance can be: 'h', 'Lag Distance', or 'Distance.' The column header containing the
            variogram values can be: 'vario', 'Variogram Value', or 'Variogram'. If reference data
            does not need to be plotted, set ``reference_data`` to ``False``
        reference_var_id (int): Point to which reference variogram you would like to plot if there are
            multiple variogram within your reference dataframe. The dataframe must contain the
            correct column ID. The column header containing the variogram index values can be:
            'Variogram Index' or 'Index'
        sill (float): Value to plot a horizontal line representing the variograms sill
        ax (mpl.axis): Matplotlib axis to plot the figure
        figsize (tuple): Figure size (width, height)
        xlim (float tuple): Minimum and maximum limits of data along the x axis
        ylim (float tuple): Minimum and maximum limits of data along the y axis
        trimylim (bool): Indicate if realization plots should cease plotting once they pass the
            limits specified by ``ylim``.
        title (str): Title for the plot. Can use ``dirs`` from :func:`gs.get_uniquevarids()
            <pygeostat.plotting.variogram_utils.get_uniquevarids>` if direction data is desired in the
            title.
        variable (str): By default, titles are generated based on the data provided, which include
            the variable number and direction information. If you would like to keep the direction
            data but update the variogram name manually, use this parameter to give that variable
            name
        xlabel (str): X-axis label
        yalabl (str): Y-axis label
        unit (str): Distance units used for lag distance. Only used if the keyword parameter
            ``xlabel`` is left to its default value of ``None``.
        refclr (str): Any Matplotlib color for the reference variogram
        simclr (str): Any Matplotlib color for the realization variograms
        alpha (float): Transparency for realization variograms (0 = Transparent, 1 = Opaque)
        ls (float): A valid Matplotlib line style for both reference and realization variograms
        lw (float): Line width in points. The width provided in this parameter is used for the
            reference variogram.
        lw_real (float): Line width in points. If no value is passed, half the value of ``lw`` is
            used for the realization variograms
        legend_label (str or bool): A string containing the label that will be attached to the
            reference variogram if it exists. If the value `True` is passed, only the realizations
            will have a label.
        plot_style (str): Use a predefined set of matplotlib plotting parameters as specified by
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Use ``False`` or ``None``
            to turn it off
        custom_style (dict): Alter some of the predefined parameters in the ``plot_style`` selected.
        output_file (str): Output figure file name and location
        out_kws (dict): Optional dictionary of permissible keyword arguments to pass to
            :func:`gs.export_image() <pygeostat.plotting.export_image.export_image>`
        **kwargs: Optional permissible keyword arguments to pass to :func:`gs.variogram_plot()
            <pygeostat.plotting.variogram_plot.variogram_plot>` and by extension, matplotlib's plot function if the
            keyword passed is not used by :func:`gs.variogram_plot() <pygeostat.plotting.variogram_plot.variogram_plot>`

    Returns:
        ax (ax): matplotlib Axes object with the variogram

    Examples:
        Import a variogram model manually into python using pygeostat functions:

        >>> #Import the reference variogram models into python
        >>> var1 = gs.VargModel(vargstr='''3 0
        >>> 1 0.46 0 0 0
        >>> 16 16 16
        >>> 1 0.06 0 0 0
        >>> 32 32 32
        >>> 1 0.48 0 0 0
        >>> 64 64 64''')
        >>> #Generate a dataframe that has lag distances values for the above model
        >>> model1 = var1.model(azm=0, dip=0, nlags=200, lagdist=0.5, returnstr=False)

        Load output from varsim:

        >>> varsimulated_data = gs.DataFile('varsim_reals.out', readfl=True)

        Look at what variograms are within the variogram data:

        >>> gs.get_uniquevarids(varsimulated_data.data)
        Variogram ID: 1 ... Variable: 1, Azimuth 90, Dip 0
        Variogram ID: 2 ... Variable: 2, Azimuth 90, Dip 0
        Variogram ID: 3 ... Variable: 3, Azimuth 90, Dip 0

        A simple call using the variogram model loaded and the set of realizations associated with
        Variogram ID 1:

        >>> gs.variogram_plot_simulation(simulated_data=varsimulated_data.data, simulated_var_id=1, reference_data=model1)

        .. image:: ./figures/variogram_plot_simulation1_150.png

        |

        Use a custom variogram title, and fix the x-axis limits:

        >>> gs.variogram_plot_simulation(simulated_data=varsimulated_data.data, simulated_var_id=1, reference_data=model1, xlim=(0,150),
        ...              title='PVar1 Variogram Reproduction')

        .. image:: ./figures/variogram_plot_simulation2_150.png

        |

        To plot all the variables at once, we first need a dictionary that has model data loaded
        within it. Create a dictionary using the already loaded model and files generated from
        varcalc for the other variables:

        >>> models={}
        >>> models[1] = model1
        >>> models[2] = gs.DataFile('varsim_model_var2.out', readfl=True).data
        >>> models[3] = gs.DataFile('varsim_model_var3.out', readfl=True).data

        Now loop through the variables, use a variable name to update the title produced:

        >>> variables=['PVar1', 'PVar2', 'PVar3']
        >>> for i in range(1,4):
        >>>     var=variables[i-1]
        >>>     gs.variogram_plot_simulation(simulated_data=varsimulated_data.data, simulated_var_id=i, reference_data=models[i], variable=var,
        ...                  xlim=(0,150))

        .. image:: ./figures/variogram_plot_simulation3_150.png
        .. image:: ./figures/variogram_plot_simulation4_150.png
        .. image:: ./figures/variogram_plot_simulation5_150.png

    """
    import pygeostat as gs
    from .variogram_utils import get_uniquevarids
    from . utils import format_plot, setup_plot
    from .export_image import export_image

    # Handle dictionary defaults
    if out_kws is None:
        out_kws = dict()
    # Get the enumerate data for the datafiles that exist
    simulated_var_ids, simdirs = get_uniquevarids(simulated_data, mode='ref', source='varsim')
    if not isinstance(reference_data, bool):
        if 'Variogram Number' in reference_data.columns:
            reference_var_ids, _ = get_uniquevarids(reference_data, mode='ref')
            refindex = reference_var_ids[reference_var_id]
        else:
            refindex = None
    # Set-up plot if no axis is supplied
    fig, ax, cax = setup_plot(ax, figsize=figsize, aspect =False)

    #  Set-up some plotting defaults
    if lw_real is None:
        lw_real = lw / 2
    if ylim is None and sill == 1:
        ylim = (0, 1.2)
    # Plot simulation variograms
    if trimylim:
        if isinstance(ylim, tuple):
            simulated_data = gs.trimylim(simulated_data, ylim[1])
    if simclr is None:
        simclr = Parameters['plotting.variogram_plot_simulation.simclr']
    if alpha is None:
        alpha = Parameters['plotting.variogram_plot_simulation.alpha']
    for i, index in enumerate(simulated_var_ids[simulated_var_id]):
        if legend_label:
            if i == (len(simulated_var_ids[simulated_var_id]) - 1):
                label = 'Realizations'
            else:
                label = False
        else:
            label = False
        gs.variogram_plot(data=simulated_data, index=index, sill=False, experimental=False, ax=ax, color=simclr,
                  ylabel=False, xlabel=False, title=False, lw=lw_real, alpha=alpha, label=label,
                  plot_style=None, **kwargs)
    # Plot the reference variogram
    if refclr is None:
        refclr = Parameters['plotting.variogram_plot_simulation.refclr']
    if not isinstance(reference_data, bool):
        if legend_label is True:
            legend_label = False
        gs.variogram_plot(data=reference_data, sill=False, experimental=False, index=refindex, ax=ax, color=refclr,
                  ylabel=False, xlabel=False, title=False, lw=lw, label=legend_label,
                  plot_style=None,
                  **kwargs)
    # Plot figure and axis labels
    if unit is None:
        unit = Parameters['plotting.unit']
    if unit is None or unit == '':
        unit = ''
    else:
        unit = ' ({})'.format(unit)
    if xlabel is None:
        ax.set_xlabel(Parameters['plotting.lagname'] + unit)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel is None:
        ax.set_ylabel(r'$\gamma$   ', fontsize=mpl.rcParams['font.size'] *
                      Parameters['plotting.gammasize'], rotation=0)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title is None:
        ndirs = len(simdirs)
        if variable is None:
            var = "Variable %s" % int((simulated_var_id / ndirs) + 0.5)
        else:
            var = variable
        direction = (simulated_var_id % ndirs) - 1
        if direction < 0:
            direction = ndirs - 1
        azm = simdirs[direction][0]
        dip = simdirs[direction][1]
        ax.set_title("%s (Azm: %s°, Dip: %s°)" % (var, azm, dip))
    if title:
        ax.set_title(title)
    # Set plot limits
    if not isinstance(sill, bool):
        ax.plot(ax.get_xlim(), (sill, sill), 'k-', ls='-', lw=1.5)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if legend_label:
        ax.legend()
    # Export figure
    if output_file or ('pdfpages' in out_kws):
        export_image(output_file, **out_kws)
    #  Return axis
    return ax
