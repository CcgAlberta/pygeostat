#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

__author__ = 'pygeostat development team'
__date__ = '2020-01-04'
__version__ = '1.0.0'

import os, sys
try:
    import pygeostat as gs
except (ImportError, ModuleNotFoundError):
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname( __file__ ), r'..')))
    import pygeostat as gs
import unittest
import warnings
import subprocess


class PygeostatPlotStyleTest(unittest.TestCase):
    '''
     Test suite for pygeostat plot style settings
    '''

    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.plot_style = gs.PlotStyle


    def test_update_key(self):
        '''
        A collection of test to check the update behavior for pygeostat Parameters object
        '''

        # Correct assignment
        self.plot_style.update({'axes.grid': True})
        
        self.assertEqual(self.plot_style['axes.grid'], True)

        # Wrong Key
        with self.assertRaises(KeyError):
            self.plot_style.update({'FakeParameters':0})
    

    def test_set_default_plot_style(self):
        '''
        Test the capability of writing pygeostat plot styles into a config file under the user folder
        '''
        try:
            self.plot_style.set_systemdefault()
        except Exception as ex:
            self.fail('Unable to set system defaults \n{}'.format(str(ex)))

    def test_get_default_plot_style(self):
        '''
        test getting the default plot style
        '''

        try:
            self.plot_style.get_systemdefault()
        except Exception as ex:
            self.fail('Unable to set system defaults \n{}'.format(str(ex)))

if __name__ == '__main__':
	subprocess.call([sys.executable, '-m', 'unittest', str(__file__), '-v'])