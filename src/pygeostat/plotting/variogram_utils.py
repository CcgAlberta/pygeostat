#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A collection of utility functions for variogram plotting"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import pandas as pd
import numpy as np


def get_uniquevarids(data, mode='print', source=None):
    """
    Enumerate the unique variograms within a output gslib variogram files based on the number of
    directions and variables found. Prints the unique variograms found to the screen so a varid can
    be determined for :func:`gs.varplt() <pygeostat.plotting.varplt.varplt>` or
    :func:`gs.varpltsim() <pygeostat.plotting.varpltsim.varpltsim>`. Alternatively, a dictionary
    containing the Variogram IDs and their corresponding 'Variogram Index' values and a list
    containing the directions found can be returned.

    Files from varcalc, varmodel, and varsim are accepted.

    .. note:
        This function relies on the arrangement of data within the output files to determine some
        parameters. If it is having some difficulties, the parameter ``source`` can be used to
        possibly solve any issues.

    Parameters:
        data (pd.DataFrame):  Dataframe containing the output data from gslib variogram programs.
            The columns 'Variogram Index', 'Variogram Number', 'Calculation Azimuth', and
            'Calculation Dip' must be unedited; specifically their headers.
        mode (str): Specify the required output from the function. If ``'print'`` is passed,
            the enumerated variogram IDs will be printed with their corresponding variable and
            direction information. If ``'ref'`` is passed, a dictionary is build for each variogram
            ID containing their corresponding 'Variogram Index' values and a direction
            list``dirs``. source (str): Manually set what the data source is. Permissible input
            is: ``'varcalc'``, ``'varmodel'``, and ``'varsim'``.

    Returns:
        varids (dict): Returned if ``mode`` is set to ``'ref'``. A dictionary is build containing
        each variogram ID and their corresponding 'Variogram Index' value or values if multiple
        realizations are found.

    Returns:
        dirs (list): Returned if ``mode`` is set to ``'ref'``. A list of directions found within
        the variogram data

    Examples:
        Load variogram output into python:

        >>> varsimdat = gs.DataFile('varsim_reals.out', readfl=True)

        Check the variogram IDs that will be used by varpltsim:

        >>> gs.get_uniquevarids(varsimdat.data)
            Variogram ID: 1 ... Variable: 1, Azimuth 90, Dip 0
            Variogram ID: 2 ... Variable: 2, Azimuth 90, Dip 0
            Variogram ID: 3 ... Variable: 3, Azimuth 90, Dip 0

        Return the varids and print the index values attached to each Variogram ID:

        >>> varids, dirs = gs.get_uniquevarids(varsimdat.data, mode='ref')
        >>> print(varids[1])
        [1, 4, 7, ..., 292, 295, 298]

        The list ``varids`` can now be used to repeadidly plot :func:`gs.varplt()
        <pygeostat.plotting.varplt.varplt>`, which is what :func:`gs.varpltsim()
        <pygeostat.plotting.varpltsim.varpltsim>` does.

        Have a look at the list ``dirs`` which was generated above:

        >>> print(dirs)
        [(90, 0)]

        The output ``dirs`` can be useful in situations such as generating plot titles.

    """
    # Infer and extract the required columns
    if 'Variogram Index' in data.columns:
        vario_index = data['Variogram Index']
    else:
        print("Error: The column 'Variogram Index' could not be found")
        return
    if 'Variogram Number' in data.columns:
        vario_num = data['Variogram Number']
    else:
        print("Error: The column 'Variogram Number' could not be found")
        return
    if 'Calculation Azimuth' not in data.columns and 'Calculation Dip' not in data.columns:
        print("Error: One or more of the columns 'Calculation Azimuth' and 'Calculation Dip' could"
              "not be found")
        return
    # Get the number of directions
    dirs = []
    for direction, _ in data.groupby(['Calculation Azimuth', 'Calculation Dip'], sort=False):
        dirs.append(direction)
    ndir = len(dirs)
    # Get the number of lags per variogram
    nlegs = len(vario_index[vario_index == 1])
    # Get the number of variables
    nvarios = len(vario_index.unique())
    nvar = len(vario_num.unique()) // ndir
    if source == 'varmodel' or nvarios == 1 or nvarios == ndir:
        nvar = 1
        nreal = 1
    elif source == 'varsim':
        nreal = vario_num[vario_num == 1].count() // nlegs
    else:
        nreal = vario_num[vario_num == 1].count() // nlegs // ndir
    # Check that the number of variograms in the file match what should be there.
    calcnvarios = nreal * nvar * ndir
    if calcnvarios != nvarios:
        print("Error: The number of directions, variables, and realizations do not match the",
              "number of\n       variograms found in the file")
        return
    # Build a reference list of variogram IDs and the 'Variogram Index' values that correspond to
    # that variogram ID
    varids = {}
    i = 1
    for j in range(1, nvar + 1):
        for k in range(0, ndir):
            varids[i] = list(range(i, nvarios + 1)[::(nvar * ndir)])
            if mode == 'print':
                print("Variogram ID: %s ... Variable: %s, Azimuth %s, Dip %s"
                      % (i, j, dirs[k][0], dirs[k][1]))
            i += 1
    # Return some variables if in the right mode
    if mode == 'ref':
        return varids, dirs


def trimylim(simdat, ylim):
    """
    Remove realization experimental variogram points beyond a variogram value so they won't be
    plotted.

    Parameters:
        simdat (pd.DataFrame):
        ylim (float): Value to trim points after
    """
    # Get the column containing the index
    if 'Variogram Index' in simdat.columns:
        index_name = 'Variogram Index'
    elif 'Index' in simdat.columns:
        index_name = 'Index'
    trimdata = []
    for index in np.unique(simdat[index_name]):
        # Figure out where to start trimming
        data = simdat[simdat[index_name] == index]
        indices = data[data.rolling(center=False, window=5).mean()['Variogram Value'] > ylim]
        indices = indices.index.tolist()
        # Trim if required
        if len(indices) > 0:
            minix = min(indices)
            # Dump the rows after the ceiling is reached
            trimdata.append(data.iloc[:minix, :])
        else:
            # Dump the data even if it doesn't need trimming
            trimdata.append(data)
    # Rebuild the trimmed DataFrame
    simdat = pd.concat(trimdata)

    return simdat
