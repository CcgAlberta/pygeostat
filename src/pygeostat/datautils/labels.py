#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''labels.py: Contains utilities for creating index columns and labels'''
__author__ = 'Tyler Acorn'
__date__ = '2015-10-08'
__version__ = '1.001'

import numpy as np


def insert_real_idx(data, num_real=0, bindex=True, real_column='Realization',
                    bi_column='BlockIndex'):
    '''
    This will insert realization index columns. By default it will use the griddef associated
    with the file.

    Parameters:
        num_real (int): If you do not have a griddef associated with the file you can tell it how
            many realizations there are
        bindex (bool): True or False for adding a block index
        real_column (str): Set the name of the column used for the Realizations Index
        bi_column (str): Set the name of the column used for the Block Index

    **Process**

        If there are already a "real_column" or "bi_column" columns it will overwrite the values in
        these columns If the "real_column" and "bi_column" columns aren't in the dataframe it will
        insert these columns at the front.

    .. codeauthor:: pygeostat development team 2015-09-30
    '''

    from ..data import DataFile
    # ---------------
    # Insert columns at front of the dataframe if they do not already exist
    # ---------------
    if real_column not in data.data.columns:
        data.data.insert(loc=0, column=real_column, value=0)
        if bindex is True:
            if bi_column not in data.columns and bindex is True:
                data.data.insert(loc=1, column=bi_column, value=0)
    # ---------------
    # Calculate Number of Realizations and the size of the realizations
    # ---------------
    if data.griddef is None:
        if num_real != 0:
            real_size = len(data.data.index) / num_real
            # The next step I make sure we have enough rows in our dataframe. is this necessary?
    elif data.griddef is not None:
        num_real = len(data.data.index) / data.griddef.count()
        real_size = data.griddef.count()
    else:
        print('ERROR: No griddef assoicated with file and num_real is still 0')
    temp_labels = make_labels('R', (int(num_real)), padding=(len(str(int(num_real)))))
    ridx = np.repeat(temp_labels, real_size)
    data.data[real_column] = ridx
    if bindex is True:
        temp_labels = make_labels('B', (int(real_size)), padding=(len(str(int(real_size)))))
        bidx = np.tile(temp_labels, num_real)
        data.data[bi_column] = bidx


def make_labels(prefix, num, padding=0):
    '''
    Returns a series of lables combining a prefix and a number with leading zeros

    Parameters:
        prefix (str): any letter(s) that you want as the prefix (for example B for blockindex)
        num (int): The number of labels you want.
        padding (int): if given an integer value will pad the numbers with zeros until the prefix +
            the numbers equal the length of the padding value

    Return:
        Series: This will return a series with "n" number of labels starting from 1

    Note:
        Barrowed from website http://pandas.pydata.org/pandas-docs/stable/advanced.html#advanced

    Examples:
        Creating an array of labels

        >>> label = gs.datautils.make_labels('R', 3, padding=3)
        >>> label
        >>> [R001, R002, R003]

    .. codeauthor:: pygeostat development team 2015-09-21
    '''

    return ["%s%s" % (prefix, str(i).zfill(padding)) for i in range(1, (num + 1))]
