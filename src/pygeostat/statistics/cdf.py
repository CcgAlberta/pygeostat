#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""cdf.py: CDF functions"""

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np


def cdf(var, weights=None, lower=None, upper=None, bins=None):
    """
    Calculates an empirical CDF using the standard method of the midpoint of a
    histogram bin. Assumes that the data is already trimmed or that iwt=0 for
    variables which are not valid.

    If 'bins' are provided, then a binned histogram approach is used. Otherwise,
    the CDF is constructed using all points (as done in GSLIB).

    Notes:
    'lower' and 'upper' limits for the CDF may be supplied and will be returned
    appropriately

    Parameters:
        var (array): Array passed to the cdf function
        weights (str): column where the weights are stored
        Lower (float): Lower limit
        upper (float): Upper Limit
        bins (int): Number of bins to use

    Returns:
        midpoints: np.ndarray
            array of bin midpoints
        cdf: np.ndarray
            cdf values for each midpoint value
    """
    # Weight normalization
    if weights is None:
        weights = np.ones(len(var)) / len(var)
    else:
        weights = weights / np.sum(weights)
    # Conversion to numpy arrays
    if not isinstance(weights, np.ndarray):
        weights = np.array(weights)
    if not isinstance(var, np.ndarray):
        var = np.array(var)
    else:
        # Make sure the shape is alright if a column array was sent
        if var.shape[0] > 1:
            var = var.transpose()
    # Binned experimental CDF
    if bins is not None:
        if lower is None:
            lower = np.min(var)
        if upper is None:
            upper = np.max(var)
        cdf, bin_edges = np.histogram(var, weights=weights, bins=bins, range=[lower, upper])
        midpoints = np.zeros(len(bin_edges) - 1)
        for idx, upper_edge in enumerate(bin_edges[1:]):
            midpoints[idx] = 0.5 * (upper_edge + bin_edges[idx])
        cdf = np.cumsum(cdf)
    # GSLIB-style all-data experimental CDF
    else:
        order = var.argsort()
        midpoints = var[order]
        cdf = np.cumsum(weights[order])
        cdf = cdf - cdf[0] / 2.0
    # Add lower and upper values if desired
    if lower is not None:
        cdf = np.append([0.0], cdf)
        midpoints = np.append([lower], midpoints)
    if upper is not None:
        cdf = np.append(cdf, [1.0])
        midpoints = np.append(midpoints, [upper])

    return(midpoints, cdf)


def percentile_from_cdf(cdf_x, cdf, percentile):
    """
    Given 'x' values of a CDF and corresponding CDF values, find a given
    percentile. Percentile may be a single value or an array-like and must be in
    [0, 100] or CDF bounds
    """
    # Handle multiple percentiles
    percentile = np.array(percentile)
    # Assumes percentiles are in [0, 100]
    percentile = percentile / 100.0
    assert((percentile >= cdf[0]).all() and (percentile <= cdf[-1]).all()), \
        'Percentile {} must be in cdf bounds'.format(percentile)
    # Piece-wise linear interpolation
    xvals = np.interp(percentile, cdf, cdf_x)
    # Return float if 1-length array
    try:
        if len(xvals) == 1:
            xvals = xvals[0]
    except TypeError:
        pass
    return xvals


def z_percentile(z, cdf_x, cdf):
    """
    Given `'cdf_x`` values of a CDF and corresponding ``cdf`` values, find the percetile of a
    given value ``z``. Percentile may be a single value or an array-like and must be in [0, 100]
    or CDF bounds.
    """
    # Sanity check
    if not ((z >= cdf_x[0]).all() and (z <= cdf_x[-1]).all()):
        raise ValueError("The value `z` must be within the bounds of the array `cdf_x`")
    # Piece-wise linear interpolation
    pvals = np.interp(z, cdf_x, cdf)
    # Return float if 1-length array
    try:
        if len(pvals) == 1:
            pvals = pvals[0]
    except TypeError:
        pass
    return pvals


def variance_from_cdf(cdfx, cdfy, nsample=1000):
    """
    Compute the variance by randomly sampling the cdf by brute force
    """
    return stdev_from_cdf(cdfx, cdfy, nsample) ** 2


def stdev_from_cdf(cdfx, cdfy, nsample=1000):
    """
    Compute the stdsev of the cdf from a n-sample sized random sample from the cdf
    """
    randsamp = np.random.rand(nsample)
    rand_z = np.interp(randsamp, cdfy, cdfx)
    return np.std(rand_z)


def build_indicator_cdf(prob_ik, zvals):
    """
    Build the X-Y data required to plot a categorical cdf

    Parameters:
        prob_ik: np.ndarray
            the p-vals corresponding to the cutoffs
        zvals: np.ndarray
            the corresponding z-value specifying the cutoffs

    Returns:
        points, midpoints: np.ndarray
            the x and y coordinates of (1) the cutoffs and (2) the midpoints for each cutoff

    """
    if prob_ik[0] != 0:
        raise ValueError('Must pass a licit CDF to this function, e.g. [0, 1]')
    else:
        points = []
        ipts = []
    for i, x in enumerate(zvals):
        if i == 0:
            points.append([x, 0])
        else:
            points.append([x, prob_ik[i - 1]])
        points.append([x, prob_ik[i]])
        ipts.append([x, (points[-1][1] + points[-2][1]) / 2])
    points = np.array(points)
    ipts = np.array(ipts)
    return points, ipts