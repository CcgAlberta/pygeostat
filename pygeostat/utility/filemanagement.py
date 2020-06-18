#!/usr/bin/env python
# -*- coding: utf-8 -*

'''filemanagement.py: Contains basic utility functions for creating/removing files and directories'''
from __future__ import absolute_import, division, print_function

__author__ = 'pygeostat development team'
__date__ = '2018'
__version__ = '1.000'


def mkdir(dirnames, verbose=False):
    '''
    Make directory(s).

    Parameters:
        dirnames(str or list): directory(s) to create
        verbose(bool): print to screen if an error is encountered,
            including if the dirnames already exists
    '''
    from os import makedirs
    if isinstance(dirnames, str):
        dirnames = [dirnames]
    for dirname in dirnames:
        try:
            makedirs(dirname)
        except:
            if verbose:
                print('{} not created - likely exists already!'.format(dirname))


def rmdir(dirnames, verbose=False):
    '''
    Remove directory(s).

    Parameters:
        dirnames(str or list): directory(s) to remove
        verbose(bool): print to screen if an error is encountered,
            including if the dirnames has already been removed
    '''
    from shutil import rmtree
    if isinstance(dirnames, str):
        dirnames = [dirnames]
    for dirname in dirnames:
        try:
            rmtree(dirname)
        except:
            if verbose:
                print('{} not removed - likely does not exist!'.format(dirname))


def rmfile(filenames, verbose=False):
    '''
    Remove file(s).

    Parameters:
        filenames(str): file to remove
        verbose(bool): print to screen if an error is encountered,
            including if the file has already been removed
    '''
    from os import remove
    if isinstance(filenames, str):
        filenames = [filenames]
    for filename in filenames:
        try:
            remove(filename)
        except:
            if verbose:
                print('{} not removed - likely does not exist!'.format(filename))
