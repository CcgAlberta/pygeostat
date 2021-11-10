# !/usr/bin/env python
#  -*- coding: utf-8 -*-
# public
def subsample(datafl, col, ncell, nsub, nreal, rseed):
    """
    When datasets become too large to efficiently use, a sub-sample maybe extracted and used in
    ieu. A Fisher-Yates shuffle is implemented.This subroutine was designed with the intent of
    wrapping it for usein python. Motivated from Jarred L. Deutsch's use of Fisher-Yates shuffle in
    gslib program histpltsim.

    Assumes that the data file being read is in GSLIB format. Outputs an 1D array in the case that
    only one realization is sub-sampled, or a2D array in the case where multiple realizations are
    sub-sampled.Each realizations sub-sample is within a unique column.

    Parameters:
        datafl (str): A single input datafile with the variable being sampled
        col (int): Column containing data to sub-sample
        ncell (int): The number of cells within a single realization, can also be interpreted as the
            number of data to read, then sub-sample
        nsub (int): Number of sub-samples to output
        nreal (int): Number of realizations to sub-sample
        rseed (int): A seed value to pass to the random number generator

    Returns:
        subsamp (np.array): Output array with a realization in each column


    Examples:
        First, the function needs to be loaded into Python:

        >>> from pygeostat.fortran.subsample import subsample

        A simple call:

        >>> subsamp = subsample(datafl='./sgsim.out', col=1, ncell=griddef.count(), nsub=5000,
        ...                     nreal=100, rseed=gs.rseed())

    .. codeauthor:: Warren E. Black - 2015-09-22
    """
    from pygeostat.fortran.subsample import subsample
    subsamp = subsample(datafl,col,ncell,nsub,nreal,rseed)

    return subsamp
