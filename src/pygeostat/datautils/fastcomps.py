#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''fastcomposites.py: Faster compositing with length saving for usage in modeling'''

__author__ = 'pygeostat development team'
__date__ = '2016-12-17'
__version__ = '1.000'

import numpy as np
import pandas as pd
from ..data import DataFile
from numba import jit


def get_composites(composites, datafl, vartypes='continuous', null=-999.0):
    '''Returns a pandas DataFrame with upscaled composites.

    Parameters:
        datafl (DataFile): pygslib DataFile with dh, ifrom, ito and at least 1 variable to upscale
        composites (DataFrame): pandas DataFrame with datafl.dh, datafl.ifrom, datafl.ito
        corresponding to the composite locations

    Keyword Args:
        vartypes (str OR dict): 'continuous', 'categorical' or a dictionary of variables like:
        {'Cu':'continuous','Facies','categorial'}
        null: the null value (values less than or equal to this will not be used)

    Returns:
        upscaled (DataFrame): Pandas DataFrame with values of upscaled variable

    .. codeauthor: pygeostat development team 2014-04-03'''

    # Get a list of drill holes for composites
    dhids = sorted(list(set(composites[datafl.dh])))

    # Check that all drillholes in the composite list are in the datafl
    datafl_dhids = list(set(datafl.data[datafl.dh]))
    for dhid in dhids:
        assert(dhid in datafl_dhids), 'ERROR: {} not in datafl'.format(dhid)

    # Determine the number of composites
    ncomposites = len(composites[datafl.dh])

    # Create a dictionary of variable upscaling methods if not supplied
    if isinstance(vartypes, str):
        vartypedict = {}
        for var in datafl.data.columns:
            if var not in [datafl.dh, datafl.ifrom, datafl.ito]:
                vartypedict[var] = vartypes.lower()
    else:
        vartypedict = vartypes

    # Determine number of variables
    nvar = len(vartypedict)

    # Convert categorical to numbers
    cat_to_numeric_dict = {}
    for var, vartype in vartypedict.items():
        if vartype.startswith('cat'):
            catlist = [cat for cat in list(set(datafl.data[var])) if cat != null]
            catdict = dict((cat, i) for cat, i in zip(catlist, range(1,len(catlist)+1)))
            catdict[null] = null
            cat_to_numeric_dict[var] = catdict
            datafl.data[var] = datafl.data[var].apply(catdict.get)


    # Group data by drill hole ID
    grouped_composites = composites.groupby([datafl.dh], squeeze=True, as_index=False)
    grouped_data = datafl.data.groupby([datafl.dh], squeeze=True, as_index=False)

    upscaled = composites.copy()

    for var, vartype in vartypedict.items():
        comp_arrs = []
        for dhid, compdf in grouped_composites:
            try:
                datadf = grouped_data.get_group(dhid)
            except KeyError:
                raise Exception('No sample data for {}'.format(dhid))

            # Composited array needs to be N x 4 (from, to, variable value, comp length)
            comp_arr = np.zeros((len(compdf), 4))
            comp_arr[:,:-2] = compdf[[datafl.ifrom, datafl.ito]].values

            # Data array only contains from, to, variable value
            data_arr = datadf[[datafl.ifrom, datafl.ito, var]].values

            # Composite and get the value and length
            comp_arr[:,2:] = np.apply_along_axis(comp_value_and_length,
                                      1,
                                      comp_arr,
                                      data_arr,
                                      vartype,
                                      null)
            comp_arrs.append(comp_arr)
        comp_arr = np.concatenate(comp_arrs)
        upscaled[var] = comp_arr[:,2]
        upscaled[var+'_length'] = comp_arr[:,3]

    # Map categorical variables back
    for var, vartype in vartypedict.items():
        if vartype.startswith('cat'):
            revdict = {v:k for k,v in cat_to_numeric_dict[var].items()}
            datafl.data[var] = datafl.data[var].apply(revdict.get)
            upscaled[var] = upscaled[var].apply(revdict.get)

    return(upscaled)

@jit(nopython=True)
def row_length_in_interval(datarow, compfrom, compto, null):
    'Finds the length in the interval along the drill hole - tda'
    # Is the variable value valid?
    if datarow[2] == null:
        return 0.0
    # Case 1 - entire sample interval in composite
    if (datarow[0] >= compfrom) and (datarow[1] <= compto):
        return datarow[1] - datarow[0]
    # Case 2 - partially bottom of interval
    elif ((datarow[0] >= compfrom) and (datarow[0] < compto) and
          (datarow[1] > compto)):

        return compto - datarow[0]
    # Case 3 - partially in top of interval
    elif ((datarow[0] < compfrom) and (datarow[1] > compfrom) and
          (datarow[1] <= compto)):

        return datarow[1] - compfrom
    # Case 4 - composite bracketed by sample interval
    elif (datarow[0] <= compfrom) and (datarow[1] >= compto):
        return compto - compfrom
    # Case 5 - no length in interval
    else:
        return 0.0

@jit
def comp_value_and_length(comprow, data_arr, vartype, null):
    '''Return the composited value and length'''
    lengths = np.apply_along_axis(row_length_in_interval,
                                  1,
                                  data_arr,
                                  comprow[0],
                                  comprow[1],
                                  null)
    total_length = np.sum(lengths)
    if np.sum(lengths) <= 0.0:
        return((null, total_length))
    if vartype.startswith('con'):
        return((np.dot(data_arr[:,2], lengths)/total_length), total_length)
    else:
        unique_cats = list(set(data_arr[lengths > 0.0, 2]))
        if len(unique_cats) == 1:
            return((unique_cats[0], total_length))
        else:
            unique_cat_wts = np.zeros(len(unique_cats))
            # Note: list comprehension is not used to permit usage of JIT
            for i, cat in enumerate(unique_cats):
                unique_cat_wts[i] = np.sum(lengths[data_arr[:,2] == cat])
            unique_cat_wts = unique_cat_wts.tolist()
        return((unique_cats[unique_cat_wts.index(max(unique_cat_wts))], total_length))
