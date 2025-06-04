#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''utils.py: Contains utility functions for providing status alerts in a script'''
__author__ = 'pygeostat development team'
__date__ = '2015'
__version__ = '1.000'

import textwrap


def printerr(text, width=80, errtype=None):
    """
    Small utility to print custom errors with proper indentation and text wrapping. The only error
    types coded are `error` and `warning`.

    Parameters:
        text (str): String of text without the preceding error flag that will be formated
        width (int): The maximum length of wrapped lines
        errtype (str): Indicate which type of error to format the string as. The default value of
            `None` will only text wrap the string to the specified `width`.

    .. codeauthor:: pygeostat development team 2015-11-01
    """
    import pygeostat as gs

    if isinstance(errtype, str):
        errtype.lower()
    if errtype == 'error':
        text = 'ERROR: ' + text
        subsequent_indent = "       "
    elif errtype == 'warning':
        text = 'WARNING: ' + text
        subsequent_indent = "         "
    else:
        subsequent_indent = ""
    print(textwrap.fill(text, width=width, subsequent_indent=subsequent_indent))


def log_progress(sequence, size=None, name='Items'):
    """
    Parameters:
        sequence: The thing that is being iterated over

    Example:
        >>> from pygeostat import log_progress
        >>> from time import sleep
        >>> for i in log_progress(range(200)):
        ...     sleep(0.1)

    .. image:: ./figures/log_progress.gif

    """
    try:
        from tqdm import tqdm, tqdm_notebook
    except:
        print("pip install tqdm to install the progressbar functionality !")
    try:
        get_ipython()
        return tqdm_notebook(sequence, desc=name, total=size)
    except:
        return tqdm(sequence, name, total=size, ascii=True)