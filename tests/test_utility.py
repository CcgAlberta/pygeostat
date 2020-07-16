#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

__author__ = 'pygeostat development team'
__date__ = '2020-01-04'
__version__ = '1.0.0'

import os, sys
try:
    import pygeostat as gs
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname( __file__ ), r'..')))
    import pygeostat as gs

import unittest
import warnings
import subprocess


class GetExecutableTest(unittest.TestCase):


    '''
    Test suite for loading executable files from a protected repository
    '''

    def setUp(self):

        self.wrong_token = 'wrong_token'


    def test_token_error(self):

        '''
        Test if the proper error handleing is in place for wrong access token
        '''

        with self.assertRaises(Exception):
            gs.get_executable(access_token=self.wrong_token)


if __name__ == '__main__':
	subprocess.call([sys.executable, '-m', 'unittest', str(__file__), '-v'])