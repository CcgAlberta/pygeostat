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


def accsim(truth, reals, pinc=0.05):
    """
    Calculates the proportion of locations where the true value falls within symmetric p-PI
    intervals when completing a jackknife study. A portion of the data is excluded from the
    conditioning dataset and the excluded sample locations simulated values are then checked.

    .. seealso::

        Pyrcz, M. J., & Deutsch, C. V. (2014). Geostatistical Reservoir Modeling (2nd ed.). New
        York, NY: Oxford University Press, p. 350-351.

    Arguments:
        truth: Tidy (long-form) 1D data where a single column containing the true values.
            A pandas dataframe/series or numpy array can be passed
        reals: Tidy (long-form) 2D data where a single column contains values from a single
            realizations and each row contains the simulated values from a single truth location.
            A pandas dataframe or numpy matrix can be passed

    Keyword Arguments:
        pinc (float): Increments between the probability intervals to calculate within (0, 1)

    Returns:
        propavg (pd.DataFrame): Dataframe with the calculated probability intervals and the
        fraction within the interval

    Returns:
        sumstats (dict): Dictionary containing the average variance (U), mean squared error (MSE),
        accuracy measure (acc), precision measure (pre), and a goodness measure (goo)

    """
    import pandas as pd
    import pygeostat as gs
    # Handle input
    if isinstance(truth, pd.Series):
        truth = truth.values
    elif isinstance(truth, pd.DataFrame):
        truth = truth.values
    elif not isinstance(truth, np.ndarray):
        raise ValueError("The argument `truth` must be a pd.DataFrame, pd.Series, or np.matrix")
    if isinstance(truth, np.ndarray) and len(truth.shape) == 1:
        truth = np.reshape(truth, (truth.shape[0], 1))
    if isinstance(reals, pd.DataFrame):
        reals = reals.values
    elif not isinstance(reals, np.ndarray):
        raise ValueError("The argument `reals` must be a pd.DataFrame or np.matrix")
    try:
        data = np.concatenate((truth, reals), axis=1)
        data = pd.DataFrame(data=data)
    except:
        raise ValueError("The `truth` and `reals` data could not be coerced into a pd.DataFrame")
    # Initialize some variables
    pints = np.arange(pinc, 1, pinc)
    propindic = dict([pint, []] for pint in pints)
    variances = []
    acc = dict([pint, 0] for pint in pints)
    pre = dict([pint, 0] for pint in pints)
    goo = dict([pint, 0] for pint in pints)
    # Calculate the indicator responses and local variances
    for i, values in data.iterrows():
        cdf = gs.cdf(values[1:].values)
        variances.append(np.var(values[1:].values))
        for pint in pints:
            if cdf[0][0] <= values[0] <= cdf[0][-1]:
                p = gs.z_percentile(values[0], cdf[0], cdf[1])
                plower = 0.5 - (pint / 2)
                pupper = 0.5 + (pint / 2)
                if plower <= p <= pupper:
                    indic = 1
                else:
                    indic = 0
            else:
                indic = 0
            propindic[pint].append(indic)
    # Calculate the average proportions and average variance
    propavg = []
    for pint in pints:
        avg = np.average(propindic[pint])
        propavg.append([pint, avg])
    propavg = pd.DataFrame(propavg, columns=['ProbInt', 'FracIn'])
    # Calculate the summary statistics
    avgvar = np.average(variances)
    mse = ((propavg['ProbInt'].values - propavg['FracIn'].values) ** 2).mean()
    acc = 0
    pre = 0
    goo = 0
    for i, values in propavg.iterrows():
        if values.iloc[1] >= values.iloc[0]:
            acc = acc + 1
            pre = pre + (values.iloc[1] - values.iloc[0])
            goo = goo + (values.iloc[1] - values.iloc[0])
        else:
            goo = goo + (2 * (values.iloc[0] - values.iloc[1]))
    acc = acc / len(propavg)
    pre = 1 - ((2 * pre) / len(propavg))
    goo = 1 - (goo / len(propavg))
    sumstats = {'avgvar': avgvar, 'mse': mse, 'acc': acc, 'pre': pre, 'goo': goo}

    return propavg, sumstats

def accmik(truth, thresholds, mikprobs, pinc=0.05):
    """
    Similar to accsim but accepting mik distributions instead

    Mostly pulled from accsim

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
    # Handle input
    if isinstance(truth, pd.Series):
        truth = truth.values
    elif isinstance(truth, pd.DataFrame):
        truth = truth.values
    elif not isinstance(truth, np.ndarray):
        raise ValueError("The argument `truth` must be a pd.DataFrame, pd.Series, or np.matrix")
    if isinstance(truth, np.ndarray) and len(truth.shape) == 1:
        truth = np.reshape(truth, (truth.shape[0], 1))
    if isinstance(mikprobs, pd.DataFrame):
        mikprobs = mikprobs.values
    elif not isinstance(mikprobs, np.ndarray):
        raise ValueError("The argument `mikprobs` must be a pd.DataFrame or np.matrix")
    # Initialize some variables
    pints, propindic, variances = _interval_responses(truth, mikprobs, pinc, cdf_x=thresholds)
    # Calculate the average proportions and average variance
    propavg = []
    for pint in pints:
        avg = np.average(propindic[pint])
        propavg.append([pint, avg])
    propavg = pd.DataFrame(propavg, columns=['ProbInt', 'FracIn'])
    # Calculate the summary statistics
    avgvar = np.average(variances)
    mse = ((propavg['ProbInt'].values - propavg['FracIn'].values) ** 2).mean()
    acc = 0
    pre = 0
    goo = 0
    for i, values in propavg.iterrows():
        if values[1] >= values[0]:
            acc = acc + 1
            pre = pre + (values[1] - values[0])
            goo = goo + (values[1] - values[0])
        else:
            goo = goo + (2 * (values[0] - values[1]))
    acc = acc / len(propavg)
    pre = 1 - ((2 * pre) / len(propavg))
    goo = 1 - (goo / len(propavg))
    sumstats = {'avgvar': avgvar, 'mse': mse, 'acc': acc, 'pre': pre, 'goo': goo}
    return propavg, sumstats


def _interval_responses(truth, reals, pinc, cdf_x=None):
    """
    When cdf_x is None, reals contains the simulated values from which a cdf should be computed.
    Otherwise the ``'reals'`` contains the distribution F(cdf_x) values for each location (nloc,
    nquant)

    Mostly pulled from the original accsim

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
    from .cdf import variance_from_cdf
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
    if cdf_x is not None and reals[0, 0] != 0:
        reals = np.c_[np.ones(truth.shape[0]), reals]
        cdf_x = np.insert(cdf_x, [0], cdf_x[0] - (cdf_x[1] - cdf_x[0]))
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
            variances.append(variance_from_cdf(ecdf[0], ecdf[1]))
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
