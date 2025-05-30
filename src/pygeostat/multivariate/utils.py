#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Various tools to aid in multivariate geostatistics.
"""
__author__ = 'pygeostat development team'
__date__ = '2016-03-16'
__version__ = '1.000'

import numpy as np
import pandas as pd


def mds(data, variables=None):
    """
    Python implementation of the MDS coordinates calculated by the CCG program corrmat when set to
    'ordination' mode.

    MDS coordinates are calculate from the correlation matrix

    .. seealso::

        1. Deutsch, M. V, & Deutsch, C. V. (2013). A Program to Calculate and Display Correlation
           Matrices with Relevant Analysis. Edmonton AB. Retrieved from http://www.ccgalberta.com

    Parameters:
        data: Tidy (long-form) 2-D data where each column is a variable and each row is an
            observation. A pandas dataframe or numpy array may be passed.
        variables (list): Variables from the pd.DataFrame passed with ``data`` to calculate
            coordinates for

    Returns:
        coords (pd.DataFrame): The 3-D MDS coordinates calculated

    .. codeauthor:: pygeostat development team 2016-05-30
    """
    # Calculate the correlation matrix
    if isinstance(data, pd.DataFrame) and (variables is not None):
        corrmat = data[variables].corr()
    else:
        corrmat = np.corrcoef(data, rowvar=0)
    nvar = corrmat.shape[0]
    # Calculate the MDS coordinates as done in the program 'corrmat'
    proxmat = 1 - corrmat
    imat = np.identity(nvar)
    jmat = imat - (1 / nvar)
    bmat = np.dot(jmat, proxmat)
    bmat = -0.5 * np.dot(bmat, jmat)
    eigval, eigvec = np.linalg.eig(bmat)
    eigval.sort()
    eigval = eigval[::-1]
    eigvalmat = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i == j:
                eigvalmat[i, j] = np.sqrt(abs(eigval[i]))
    coords = np.dot(eigvec[:, :3], eigvalmat)
    coords[:, 0] = coords[:, 0] * -1
    coords[:, 2] = coords[:, 2] * -1
    # Save the coordinates as a pd.DataFrame and return
    coords = pd.DataFrame(data=coords, columns=['X1', 'X2', 'X3'])

    return coords
