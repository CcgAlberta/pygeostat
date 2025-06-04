#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''A collection of utility functions for data prosessing and manipulation'''

import os
import math
import copy
import numpy as np
import pandas as pd
from ..data import GridDef as GridDef
from ..data import DataFile as DataFile

def round_sigfig(value, sigfigs):
    """
    Round a float or integer to a specified number of significant figures. Also handles effectively
    zero, infinity, and negative infinity values.\n
    From: http://stackoverflow.com/questions/3410976/

    Parameters:
        value (int or float): Value that requires rounding
        sigfigs (int): Number of significant figures to round the value to

    Returns:
        new_value (int or float): Rounded value

    Example:
        >>> gs.round_sigfig(-0.00032161, 3)
        >>> -0.00322

    .. codeauthor:: pygeostat development team 2015-10-13
    """
    if value == 0 or math.isnan(value):
        return value
    else:
        new_value = round(value, -int(math.floor(math.log10(abs(value)))) + (sigfigs - 1))
    if abs(new_value) < 1e-15 and abs(new_value) >= 0:
        new_value = 0.00
    elif new_value > 1e+15:
        new_value = float("inf")
    elif new_value < -1e+15:
        new_value = float("-inf")

    return new_value


def fileheader(datafl, mute=False):
    """
    Read a GSLIB file from python and return the header information. Useful for large files.

    .. codeauthor:: pygeostat development team 2016-02-15
    """
    # Get the header and print out if required
    nvar = 999
    variables = []
    with open(datafl) as file:
        for i, line in enumerate(file):
            line = line.replace('\n', '')
            if i == 0:
                if not mute:
                    print(line)
                header = line
            elif i == 1:
                if not mute:
                    print(line)
                line = line.split()
                nvar = int(line[0])
            elif i <= nvar + 1:
                if not mute:
                    print(line)
                variables.append(line)
            else:
                break
            i += 1
    # Return the data if required
    if mute:
        return header, nvar, variables


def corrmatstr(corrmat, fmt):
    """
    Converts a correlation matrix that is currently a numpy matrix or a pandas dataframe, into a
    space delimited string. Correlation matrix strings are required in the parameter files of CCG
    programs such as USGSIM and supersec.

    Currently, this function is hard coded to return two formats as specified by the ``fmt``
    argument, one for ``'usgsim'`` and one for ``'supersec'``. ``'usgsim'`` returns the full
    correlation matrix while ``'supersec'`` returns only the upper triangle of the matrix, without
    the diagonal values.

    Parameters:
        corrmat: Correlation matrix as either a pandas dataframe (pd.DataFrame) or numpy matrix
            (np.ndarray).
        fmt (str): Indicate which format to return. Accepts only one of ``['usgsim', 'supersec']``

    Returns:
        corrstr (str): Correlation matrix as a space delimited string.

    .. codeauthor:: pygeostat development team 2016-03-15
    """
    # Sanity Check
    if fmt not in ['supersec', 'usgsim']:
        return ValueError("One of 'supersec' or 'usgsim' must be passed")
    # Convert the correlation matrix to a string
    if not isinstance(corrmat, pd.DataFrame):
        corrmat = pd.DataFrame(data=corrmat)
    if fmt == 'supersec':
        if isinstance(corrmat, pd.DataFrame):
            corrmat = corrmat.as_matrix()
        ind = np.triu_indices(len(corrmat), 1)
        i = 0
        corrstr = ""
        for j in range(len(corrmat) - 1, 0, -1):
            for k in range(j):
                corrstr = corrstr + "%s " % corrmat[ind][i]
                i += 1
            corrstr = corrstr + "\n"
    elif (fmt == 'usgsim') or fmt:
        corrstr = ""
        for i in range(0, len(corrmat)):
            string = (' '.join(map(str, corrmat.ix[i, :])))
            corrstr = corrstr + string + '\n'

    return corrstr


def slice_grid(data, griddef, orient, slice_number, slice_thickness=None, tmin=None):
    """
    Slice a 3-D grid.

    Parameters:
        data: 1-D array or a tidy long-form dataframe with a single column containing the variable
            in question and each row is an observation
        griddef (GridDef): A pygeostat GridDef class created using
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`
        orient (str): Orientation to slice data. ``'xy'``, ``'xz'``, ``'yz'`` are the only accepted
            values
        slice_number (int): Grid cell location along the axis not plotted to take the slice of data to
            plot
        slice_thickness (int): Number of slices around the selected slice_number to average the attributes
        tmin (float): A minimum threshold value to filter/key out the data

    Returns:
        view (np.ndarray): 1-D array of the sliced data

    .. codeauthor:: pygeostat development team 2014-04-19
    """
    if not isinstance(slice_number, int):
        slice_number = int(slice_number)
    if data.shape != (griddef.nz, griddef.ny, griddef.nx):
        data = data.reshape((griddef.nz, griddef.ny, griddef.nx))
    if orient == 'xy':
        view = data[slice_number, :, :]
    elif orient == 'xz':
        view = data[:, slice_number, :]
    elif orient == 'yz':
        view = data[:, :, slice_number]
    else:
        raise ValueError("Not a valid orientation! {}".format(orient))
    # new functionality to get the average of a set of slices, honoring null values if present
    if slice_thickness is not None:
        slicerange = np.arange(slice_number - slice_thickness, slice_number + slice_thickness)
        nslice = len(slicerange)
        slices = np.zeros((nslice, *view.shape), dtype=float)
        for isl, slno in enumerate(slicerange):
            if orient == 'xy':
                slices[isl, :, :] = data[slno, :, :]
            elif orient == 'xz':
                slices[isl, :, :] = data[:, slno, :]
            elif orient == 'yz':
                slices[isl, :, :] = data[:, :, slno]
        if tmin is None:
            view = np.mean(slices, axis=0)
        else:
            slices = np.ma.masked_array(slices, mask=(slices <= tmin))
            view = np.ma.mean(slices, axis=0)
    return view


def slicescatter(data, orient, slice_number, slicetol, griddef=None, x=None, y=None, z=None):
    """
    Slice scattered data based on a GSLIB style grid definition.

    Parameters:
        data (pd.DataFrame or gs.DataFile): Dataframe where each column is a variable and each row
            is an observation. Must contain the coordinate columns required depending on the value
            of ``orient``. If a :class:`gs.DataFile <pygeostat.data.data.DataFile>` class is
            passed, its attribute ``griddef``, ``x``, ``y``, and ``z`` will be extracted.
        var (str): Column header of variable under investigation
        orient (str): Orientation to slice data. ``'xy'``, ``'xz'``, ``'yz'`` are the only accepted
            values
        slice_number (int): Grid cell location along the axis not plotted to take the slice of data to
            plot
        slicetol (float): Slice tolerance to plot point data (i.e. plot +/- ``slicetol`` from the
            center of the slice). Any negative value plots all data. Default is to plot all data.
        griddef (gs.GridDef): A pygeostat GridDef class created using
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`. Required if the
            attribute cannot be retrieved from ``data`` if it is a
            :class:`gs.DataFile <pygeostat.data.data.DataFile>` class.
        x (str): Column header of x-coordinate. Required if the attribute cannot be retrieved from
            ``data`` if it is a :class:`gs.DataFile <pygeostat.data.data.DataFile>` class.
        y (str): Column header of x-coordinate. Required if the attribute cannot be retrieved from
            ``data`` if it is a :class:`gs.DataFile <pygeostat.data.data.DataFile>` class.
        z (str): Column header of x-coordinate. Required if the attribute cannot be retrieved from
            ``data`` if it is a :class:`gs.DataFile <pygeostat.data.data.DataFile>` class.

    Returns:
        pointview (pd.DataFrane): pd.DataFrame of the sliced data
    .. codeauthor:: pygeostat development team - 2016-04-11
    """
    # Handle default parameters
    if isinstance(data, DataFile):
        if not griddef:
            if data.griddef:
                griddef = data.griddef
        if not x:
            x = data.x
        if not y:
            y = data.y
        if not z:
            z = data.z
        data = data.data
    if not x:
        if 'x' in data.columns:
            x = 'x'
    if not y:
        if 'y' in data.columns:
            y = 'y'
    if not z:
        if 'z' in data.columns:
            z = 'z'
    # Sanity checks
    if not any([isinstance(data, pd.DataFrame), isinstance(data, DataFile)]):
        raise ValueError("The parameter data must be either a gs.DataFile or a pd.DataFrame")
    if orient == 'xy':
            if any([x is None, y is None]):
                raise ValueError("`orient` is set to `'xy'` yet the column ID for `x` and `y` has"
                                 " not been set and could not be automatically retrieved.")
    elif orient == 'xz':
        if any([x is None, z is None]):
            raise ValueError("`orient` is set to `'xz'` yet the column ID for `x` and `z` has"
                             " not been set and could not be automatically retrieved.")
    elif orient == 'yz':
        if any([y is None, z is None]):
            raise ValueError("`orient` is set to `'yz'` yet the column ID for `y` and `z` has"
                             " not been set and could not be automatically retrieved.")
    else:
        raise ValueError("An `orient` setting of %s is not permissible" % orient)
    if not isinstance(griddef, GridDef):
        raise ValueError("A valid gs.GridDef was not passed")
    # Slice the data
    if orient == 'xy':
        if griddef.nz > 1:
            pointview = data[data[z] >= (griddef.get_slice_coordinate('xy', slice_number) - slicetol *
                                         griddef.zsiz)]
            pointview = pointview[pointview[z] <=
                                  (griddef.get_slice_coordinate('xy', slice_number) + slicetol * griddef.zsiz)]
    elif orient == 'xz':
        if griddef.ny > 1:
            pointview = data[data[y] >= (griddef.get_slice_coordinate('xz', slice_number) - slicetol *
                                         griddef.ysiz)]
            pointview = pointview[pointview[y] <=
                                  (griddef.get_slice_coordinate('xz', slice_number) + slicetol * griddef.ysiz)]
    elif orient == 'yz':
        if griddef.nx > 1:
            pointview = data[data[x] >= (griddef.get_slice_coordinate('yz', slice_number) - slicetol *
                                         griddef.xsiz)]
            pointview = pointview[pointview[x] <=
                                  (griddef.get_slice_coordinate('yz', slice_number) + slicetol * griddef.xsiz)]

    return pointview


def fixpath(path):
    """
    Convert a file path to an absolute path if required and make sure there are only forward
    slashes.

    If copying the path directly from windows explorer or something that will produce a path like
    that, make sure to indicate to python that the string is raw. This is done by placing a ``r``
    in front of the string. For example:

    >>> string = r"A string with backslashes \\ \\ \\ \\"

    Example:
        Make sure to place an ``r`` in front of the string so funny things don't happen. A simple
        call:

        >>> gs.fixpath(r"D:\\Data\\data.dat")

    .. codeauthor:: pygeostat development team - 2016-02-07
    """
    # Convert any relative paths to absolute paths if needed
    if not os.path.isabs(path):
        path = os.path.abspath(path)
    # Fix slashes if required
    path = path.replace('\\\\\\', '\\')
    path = path.replace('\\\\', '\\')
    path = path.replace('\\', '/')

    return path


def is_numeric(s):
    """Returns true if a value can be converted to a floating point number"""
    try:
        float(s)
        return True
    except ValueError:
        return False
    except TypeError:
        print('ERROR: Must have null value in working dictionary')
        return False


def ensure_dir(f):
    """Function to make sure that directory(s) exists and if not, create it"""
    import os, copy
    if isinstance(f, str):
        dirs = [f]
    else:
        dirs = copy.copy(f)
    for f in dirs:
        if not os.path.exists(f):
            os.makedirs(f)

def ensure_path(path):
    """
    Function ensures that all folders in a given path or list of paths are created if they
    do not exist
    """
    if isinstance(path, str):
        path = [path]
    path = [p.replace('\\', '/') for p in path]
    for p in path:
        incremental_path = ''
        dirlist = p.split('/')
        for d in dirlist:
            incremental_path += d + '/'
            if ':' not in d and d != '':
                ensure_dir(incremental_path)

def nearest_eucdist(x, y=None, z=None):
    """
    Calculate the euclidean distance to the nearest sample for each sample.

    Parameters:
        x (np.array): Array of the coordinate in the x direction

    Keyword Arguments:
        y (np.array): Array of the coordinate in the y direction
        z (np.array): Array of the coordinate in the z direction

    Returns:
        dist (np.array): Array of the euclidean distance to the nearest sample for each
        sample

    .. codeauthor:: pygeostat development team - 2016-07-28
    """
    import scipy
    nsamp = len(x)
    coords = []
    if y is not None and z is not None:
        for coord in zip(x, y, z):
            coords.append(coord)
    if y is not None and z is None:
        for coord in zip(x, y):
            coords.append(coord)
    if z is not None and y is None:
        for coord in zip(x, z):
            coords.append(coord)
    distmat = scipy.spatial.distance.cdist(coords, coords, 'euclidean')
    distmat[np.tril_indices(nsamp, k=0, m=nsamp)] = np.nan
    nearest = []
    for i in range(nsamp - 1):
        dist = distmat[i, :]
        nearest.append(np.nanmin(dist))
    nearest = np.array(nearest)

    return nearest
