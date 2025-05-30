#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Kernel Density Estimation (KDE) calculation module """
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from scipy.stats import gaussian_kde
from statsmodels.nonparametric.kde import KDEUnivariate
from statsmodels.nonparametric.kernel_density import KDEMultivariate
import numpy as np

def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):

    """Kernel Density Estimation based on different packages and different kernels
    Note that scipy weights its bandwidth by the covariance of the
    input data.  To make the results comparable to the other methods,
    we divide the bandwidth by the sample standard deviation here.
    """

    if bandwidth > 0:
        kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    else:
        kde = gaussian_kde(x, **kwargs)

    return kde.evaluate(x_grid)


def kde_statsmodels_u(x, x_grid, bandwidth=0.2, **kwargs):
    """Univariate Kernel Density Estimation with Statsmodels"""
    kde = KDEUnivariate(x)
    kde.fit(bw=bandwidth, **kwargs)
    return kde.evaluate(x_grid)


def kde_statsmodels_m(x, x_grid, bandwidth=0.2, **kwargs):
    """Multivariate Kernel Density Estimation with Statsmodels"""
    kde = KDEMultivariate(x, bw=bandwidth * np.ones_like(x),
                          var_type='c', **kwargs)
    return kde.pdf(x_grid)


def kde_sklearn(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    try:
        from sklearn.neighbors import KernelDensity
    except ImportError:
        raise ImportError("The plotting library Scikit-learn was not found, please install it to to enable"
                            " kde_sklern. It can be installed via Anaconda by executing the following"
                            " command in either cygwin or windows command prompt:                     "
                            "conda install Scikit-learn")
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)

# The Scipy KDE implementation contains only the common Gaussian Kernel. Statsmodels contains seven kernels,
# while Scikit-learn contains six kernels, each of which can be used with one of about a dozen distance metrics,
# resulting in a very flexible range of effective kernel shapes.

kde_funcs = [kde_statsmodels_u, kde_statsmodels_m, kde_scipy, kde_sklearn]
kde_funcnames = ['Statsmodels-U', 'Statsmodels-M', 'Scipy', 'Scikit-learn']
