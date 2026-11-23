#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various tools for calculating statistics
"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np


def weighted_mean(var, wts):
    """Calculates the weighted mean"""
    return np.average(var, weights=wts)


def weighted_variance(var, wts):
    """Calculates the weighted variance"""
    return np.average((var - weighted_mean(var, wts))**2, weights=wts)


def weighted_skew(var, wts):
    """Calculates the weighted skewness"""
    return (np.average((var - weighted_mean(var, wts))**3, weights=wts) /
            weighted_variance(var, wts)**(1.5))

def weighted_kurtosis(var, wts):
    """Calculates the weighted skewness"""
    return (np.average((var - weighted_mean(var, wts))**4, weights=wts) /
            weighted_variance(var, wts)**(2))

def weighted_covariance(x, y, wt):
    """Calculates the weighted covariance"""
    return (np.average((x - weighted_mean(x, wt)) *
                       (y - weighted_mean(y, wt)), weights=wt))


def weighted_correlation(x, y, wt):
    """Calculates the weighted correlation"""
    return (weighted_covariance(x, y, wt) /
            (np.sqrt(weighted_variance(x, wt)) * np.sqrt(weighted_variance(y, wt))))


def weighted_correlation_rank(x, y, wt):
    """Calculatees the weighted spearman rank correlation"""
    from scipy.stats import rankdata
    x = rankdata(x)
    y = rankdata(y)
    return weighted_correlation(x, y, wt)

def near_positive_definite(input_matrix):
    """This function uses R to calculate the nearest positive definite matrix
    within python. An installation of R with the library "Matrix" is required.
    The module rpy2 is also needed

    The only requirement is an input matrix. Can be either a pandas dataframe or
    numpy-array.

    Parameters:
        input_matrix:  input numpy array or pandas dataframe, not numpy matrix

    Returns:
        (np.array):  Nearest positive definite matrix as a numpy-array
    """
    import pandas as pd
    import numpy as np
    from ..utility.logging import printerr
    import pygeostat as gs
    # Try and load r2py
    try:
        import rpy2.robjects as robjects
        from rpy2.robjects import r
        from rpy2.robjects.packages import importr
    except ImportError:
        printerr(("near_positive_definite could not be loaded. Please install the r2py library"
                     " and the software R with the library 'Matrix' to enable it. Installation"
                     " instructions can be found within pygeostat's documentation."),
                    errtype='error')
        return
    # Convert input matrix to a numpy array if it is a pd.DataFrame
    if isinstance(input_matrix, pd.DataFrame):
        input_matrix = input_matrix.as_matrix()
    # Determine matrix shape
    dim = input_matrix.shape
    # Call matrix R library
    matcalc = importr("Matrix")
    # Convert numpy array to RObject then to R matrix
    pdmat = robjects.FloatVector(input_matrix.reshape((input_matrix.size)))
    pdmat = robjects.r.matrix(pdmat, nrow=dim[0], ncol=dim[1], byrow=True)
    # Calculate nearest positive definite matrix
    pdmat = matcalc.near_positive_definite(pdmat)
    # Convert calculated matrix to python string
    pdmat = pdmat[0]  # Extract near_positive_definite matrix from R object
    pdmat = r.toString(pdmat)  # Convert R binary object to a string
    pdmat = pdmat.r_repr()  # Convert R string to python string
    pdmat = pdmat.replace('"', "")  # Clean up string
    pdmat = pdmat.replace(' ', "")  # Clean up string
    # Convert comma delimited string to list then to np array
    pdmat = [float(x) for x in pdmat.split(',')]
    pdmat = np.array(pdmat)
    pdmat = np.reshape(pdmat, dim)  # Restore near_positive_definite matrix to the original input shape

    return pdmat


def accsim(truth, reals, pinc=None, reals_gaussian=None, data=None, weights=None,
           global_variance=1.0):
    """
    Calculates the proportion of locations where the true value falls within symmetric p-PI
    intervals when completing a jackknife study. A portion of the data is excluded from the
    conditioning dataset and the excluded sample locations simulated values are then checked.

    The summary measures are Clayton Deutsch's accuracy (A), precision (P) and goodness (G), all in
    percent. With ``p`` the probability intervals and ``t`` the fraction of the true values in each:

    * ``A = (1 - 2 * mean(max(0, p - t))) * 100``: between 0 and 100; intervals on or above the 1:1
      line add nothing, so only those below it lower A
    * ``P = (1 - mean(local variance) / global variance) * 100``, in Gaussian units. It is negative
      if the mean local variance is larger than the global variance, and is not clipped. The true
      values are not used
    * ``G = 2 / (1/A + 1/P)``, the harmonic average of A and P, between them when both are positive

    .. seealso::

        Pyrcz, M. J., & Deutsch, C. V. (2014). Geostatistical Reservoir Modeling (2nd ed.). New
        York, NY: Oxford University Press, p. 350-351.

        Deutsch, C. V. (2026). Refresh on the Accuracy Plot and Measures of Performance of
        Probabilistic Prediction. Paper 40, CCG Annual Report 28.

    Arguments:
        truth: Tidy (long-form) 1D data where a single column containing the true values.
            A pandas dataframe/series or numpy array can be passed
        reals: Tidy (long-form) 2D data where a single column contains values from a single
            realizations and each row contains the simulated values from a single truth location.
            A pandas dataframe or numpy matrix can be passed

    Keyword Arguments:
        pinc (float): Increments between the probability intervals within (0, 1). If None, 0.1
            (9 intervals, 10% to 90%) for fewer than 300 true values and 0.05 (19 intervals,
            5% to 95%) otherwise
        reals_gaussian: The realizations in Gaussian units (before the back transform), the same
            shape as ``reals``. Used for P and G
        data: The data values that the normal score transform of the simulation was made from. If
            ``reals_gaussian`` is None, the same transform is applied to ``reals`` to get Gaussian
            units for P and G
        weights: Declustering weights of ``data``, if the transform used them
        global_variance (float): The variance of the transformed data, 1 for a normal score transform

    Returns:
        propavg (pd.DataFrame): Dataframe with the calculated probability intervals and the
        fraction within the interval

    Returns:
        sumstats (dict): Dictionary containing the average variance (U), mean squared error (MSE),
        accuracy (acc, A), precision (pre, P) and goodness (goo, G). P and G are NaN without
        ``reals_gaussian`` or ``data``

    """
    import pandas as pd
    truth = _truth_1d(truth)
    if isinstance(reals, pd.DataFrame):
        reals = reals.values
    elif not isinstance(reals, np.ndarray):
        raise ValueError("The argument `reals` must be a pd.DataFrame or np.matrix")
    if pinc is None:
        pinc = default_probability_increment(len(truth))
    if not global_variance > 0:
        raise ValueError("`global_variance` must be positive")
    gaussian_variances = None
    if reals_gaussian is not None:
        if len(reals_gaussian) != len(reals) or any(
                np.shape(g) != np.shape(r) for g, r in zip(reals_gaussian, reals)):
            raise ValueError("`reals_gaussian` must have the shape of `reals`, {}".format(reals.shape))
        gaussian_variances = [np.var(np.asarray(g, dtype=float)) for g in reals_gaussian]
    elif data is not None:
        table = _normal_score_table(data, weights)
        # Row by row, so locations with different numbers of realizations work too
        gaussian_variances = [np.var(np.interp(np.asarray(r, dtype=float), *table)) for r in reals]
    pints, propindic, variances = _interval_responses(truth, reals, pinc)
    return _accuracy_summary(pints, propindic, variances, gaussian_variances, global_variance)


def default_probability_increment(ndata):
    """The probability increment for ``ndata`` true values: 9 intervals (10% to 90%) for fewer than
    a few hundred data (taken as 300), 19 intervals (5% to 95%) for more"""
    return 0.1 if ndata < 300 else 0.05


def _normal_score_table(data, weights=None):
    """
    The normal score transform of ``data`` as a table (values, scores), for np.interp: the data
    sorted, each with the normal quantile of the midpoint of its (weighted) cumulative probability,
    and the mean score for tied values. Values beyond the data take the lowest or highest score.
    """
    from scipy.stats import norm
    data = np.asarray(data, dtype=float).ravel()
    weights = np.ones(len(data)) if weights is None else np.asarray(weights, dtype=float).ravel()
    keep = ~np.isnan(data)
    data, weights = data[keep], weights[keep]
    order = np.argsort(data, kind='stable')
    data, weights = data[order], weights[order]
    cumulative = np.cumsum(weights)
    scores = norm.ppf((cumulative - weights / 2) / cumulative[-1])
    values, first = np.unique(data, return_index=True)
    scores = np.array([scores[a:b].mean() for a, b in zip(first, list(first[1:]) + [len(data)])])
    return values, scores


def accmik(truth, thresholds, mikprobs, pinc=None):
    """
    Similar to accsim but accepting mik distributions instead

    Parameters
    ----------
        truth: np.ndarray
            Tidy (long-form) 1D data where a single column containing the true values.
            A pandas dataframe/series or numpy array can be passed
        thresholds: np.ndarray
            1D array of thresholds where each CDF is defined by these thresholds and the probability
            given in the mikprobs array for each location.
        mikprobs: np.ndarray
            Tidy (long-form) 2D data where a single column contains values from a single
            MIK cutoff and each row contains the simulated values for the corresponding single
            truth location. A pandas dataframe or numpy matrix can be passed
        pinc: float
            Increments between the probability intervals to calculate within (0, 1)

    Returns
    -------
        propavg: pd.DataFrame
            Dataframe with the calculated probability intervals and the fraction within the interval
        sumstats: dict
            Dictionary containing the average variance (U), mean squared error (MSE), accuracy
            measure (acc), precision measure (pre), and a goodness measure (goo)

    """
    import pandas as pd
    truth = _truth_1d(truth)
    if isinstance(mikprobs, pd.DataFrame):
        mikprobs = mikprobs.values
    elif not isinstance(mikprobs, np.ndarray):
        raise ValueError("The argument `mikprobs` must be a pd.DataFrame or np.matrix")
    if pinc is None:
        pinc = default_probability_increment(len(truth))
    pints, propindic, variances = _interval_responses(truth, mikprobs, pinc, cdf_x=thresholds)
    return _accuracy_summary(pints, propindic, variances)


def _truth_1d(truth):
    """
    Return ``truth`` as a 1D array. Accepts a 1D array, a pd.Series or a single column.
    """
    import pandas as pd
    if isinstance(truth, (pd.Series, pd.DataFrame)):
        truth = truth.values
    elif not isinstance(truth, np.ndarray):
        raise ValueError("The argument `truth` must be a pd.DataFrame, pd.Series, or np.matrix")
    if truth.ndim == 2 and truth.shape[1] == 1:
        truth = truth[:, 0]
    if truth.ndim != 1:
        raise ValueError("The argument `truth` must be a single column, got shape {}".format(
            truth.shape))
    return truth


def _accuracy_summary(pints, propindic, variances, gaussian_variances=None, global_variance=1.0):
    """
    Fraction of the true values inside each probability interval, and the summary statistics
    returned by accsim and accmik: the accuracy A, precision P and goodness G (percent) of
    Clayton Deutsch. P and G need the local variances in Gaussian units (``gaussian_variances``).
    """
    import pandas as pd
    fracin = np.array([np.mean(propindic[pint]) for pint in pints])
    propavg = pd.DataFrame({'ProbInt': pints, 'FracIn': fracin})
    excess = fracin - pints
    accuracy = 100 * (1 - 2 * np.mean(np.maximum(0, -excess)))
    precision = np.nan
    goodness = np.nan
    if gaussian_variances is not None:
        precision = 100 * (1 - np.mean(gaussian_variances) / global_variance)
        # The harmonic average is not defined for a negative P, and is 0 when A or P is 0
        if precision == 0 or accuracy == 0:
            goodness = 0.0
        elif precision > 0:
            goodness = 2 / (1 / accuracy + 1 / precision)
    sumstats = {'avgvar': np.mean(variances),
                'mse': np.mean(excess ** 2),
                'acc': accuracy, 'pre': precision, 'goo': goodness}
    return propavg, sumstats


def _variance_piecewise_cdf(cdf_x, cdf_y):
    """
    Exact variance of the distribution whose CDF linearly interpolates ``cdf_y`` (from 0 to 1)
    at ``cdf_x``: each step between thresholds holds its probability spread uniformly.
    """
    cdf_x = np.asarray(cdf_x, dtype=float)
    mass = np.diff(cdf_y)
    lower, upper = cdf_x[:-1], cdf_x[1:]
    mean = np.sum(mass * (lower + upper) / 2)
    lower, upper = lower - mean, upper - mean  # centred, so large coordinates keep precision
    return np.sum(mass * (lower ** 2 + lower * upper + upper ** 2) / 3)


def _interval_responses(truth, reals, pinc, cdf_x=None):
    """
    When cdf_x is None, reals contains the simulated values from which a cdf should be computed.
    Otherwise the ``'reals'`` contains the distribution F(cdf_x) values for each location (nloc,
    nquant)

    Shared by accsim and accmik.

    Parameters:

        truth: np.ndarray
            tidy 1D array of truth values
        reals: np.ndarray
            tidy 2D array of `reals` where if cdf_x is None these are the realizations from which a cdf
            is built, otherwise cdf_x defines the z-values and each row of reals contains the
            corresponding probabilites defining the local cdf's
        pinc: float
            the incremement of the probability intervals
        cdf_x: np.ndarray, optional
            contains the z-values when ``reals`` contains the F(z)

    Returns:
        pints: np.ndarray
            a range of pinc spaced probability intervals
        propindic: dict
            the dictionary used in accsim and accmik functions
        variances: list
            the list of variances
    """
    from .cdf import cdf, z_percentile
    if not isinstance(truth, np.ndarray):
        truth = np.array(truth)
    isjagged = False
    if isinstance(reals, list) or reals.dtype == "O":
        # assume reals is `jagged` (nested lists), each location has a different # simulated vals
        isjagged = True
    elif not isinstance(reals, np.ndarray):
        reals = np.array(reals)
    if truth.shape[0] != reals.shape[0]:
        raise ValueError('`truth` and `reals` must have the same dimension along the first axis!')
    # initializse the variables
    pints = np.arange(pinc, 1, pinc)
    propindic = {pint: [] for pint in pints}
    variances = []
    if cdf_x is not None and (reals[:, 0] != 0).any():
        # Pad the lower tail one threshold step down, where every CDF is 0
        reals = np.c_[np.zeros(truth.shape[0]), reals]
        cdf_x = np.insert(cdf_x, [0], cdf_x[0] - (cdf_x[1] - cdf_x[0]))
    if cdf_x is not None and (reals[:, -1] != 1).any():
        # Pad the upper tail one threshold step up, where every CDF is 1
        reals = np.c_[reals, np.ones(truth.shape[0])]
        cdf_x = np.append(cdf_x, cdf_x[-1] + (cdf_x[-1] - cdf_x[-2]))
    # Calculate the indicator responses and local variances
    for i in range(truth.shape[0]):
        if cdf_x is None:
            if isjagged:
                ecdf = cdf(reals[i])  # each element in reals is a list of sim vals
                v = np.var(reals[i])
            else:
                ecdf = cdf(reals[i, :])  # reals is a 2D array with standard size that can be sliced
                v = np.var(reals[i, :])
            variances.append(v)
        else:
            ecdf = (cdf_x, reals[i, :])
            variances.append(_variance_piecewise_cdf(ecdf[0], ecdf[1]))
        for pint in pints:
            if ecdf[0][0] <= truth[i] <= ecdf[0][-1]:
                p = z_percentile(truth[i], ecdf[0], ecdf[1])
                plower = 0.5 - (pint / 2)
                pupper = 0.5 + (pint / 2)
                if plower <= p <= pupper:
                    indic = 1
                else:
                    indic = 0
            else:
                indic = 0
            propindic[pint].append(indic)
    return pints, propindic, variances
