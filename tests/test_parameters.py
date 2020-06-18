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


class PygeostatParametersTest(unittest.TestCase):
    '''
     Test suite for pygeostat global parameters
    '''

    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.params = gs.Parameters
        self.key = 'data.tmin'
        self.value = -789
        self.wrong_value = 'wrong'

    def test_description_exist(self):
        '''
            A test to check that all parameters have a description
        '''
        for name, value in self.params.items():
            try:
                self.params.describe(name, verbose=False)
            except KeyError:
                self.fail('The provided parameter, "{}" does not have a description'.format(name))

    def test_update_key(self):
        '''
        A collection of test to check the update behavior for pygeostat Parameters object
        '''

        # Correct assignment
        self.params.update({self.key:self.value})
        
        self.assertEqual(self.params[self.key], self.value)

        # Wrong Key
        # with self.assertRaises(KeyError):
        #     self.params.update({'FakeParameters':0})

        # Correct key but wrong value
        with self.assertRaises(ValueError):
            self.params.update({self.key:self.wrong_value})


    def test_set_default_parameters(self):
        '''
        Test the capability of writing pygeostat default parameters into a config file under the user folder
        '''
        try:
            self.params.set_systemdefault()
        except Exception as ex:
            self.fail('Unable to set system defaults \n{}'.format(str(ex)))

    def test_get_default_parameters(self):
        '''
        test getting the default parameters
        '''

        try:
            self.params.get_systemdefault()
        except Exception as ex:
            self.fail('Unable to set system defaults \n{}'.format(str(ex)))



if __name__ == '__main__':
	subprocess.call([sys.executable, '-m', 'unittest', str(__file__), '-v'])