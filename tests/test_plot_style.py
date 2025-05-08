#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import warnings

from pygeostat.plotting.set_style import PlotStyle


# Fixture component used in the tests
@pytest.fixture(scope="module")
def plot_style():
    """
    Fixture providing the PlotStyle object from gs for tests

    Returns the PlotStyle class to test its static methods and properties
    """

    warnings.simplefilter('ignore', category=ImportWarning)
    return PlotStyle


def test_update_key(plot_style):
    """Test to update plot style parameters with valid keys."""

    # Original so it can be restored later
    original_value = plot_style.get('axes.grid', None)

    try:
        # Update with valid key
        plot_style.update({'axes.grid': True})
        assert plot_style['axes.grid'] is True
    finally:
        # Restore orignal value to avoid affect rest of tests
        if original_value is not None:
            plot_style.update({'axes.grid': original_value})




#class PygeostatPlotStyleTest(unittest.TestCase):
#    '''
#     Test suite for pygeostat plot style settings
#    '''
#
#    def setUp(self):
#        warnings.simplefilter('ignore', category=ImportWarning)
#        self.plot_style = gs.PlotStyle
#
#
#    def test_update_key(self):
#        '''
#        A collection of test to check the update behavior for pygeostat Parameters object
#        '''
#
#        # Correct assignment
#        self.plot_style.update({'axes.grid': True})
#        
#        self.assertEqual(self.plot_style['axes.grid'], True)
#
#        # Wrong Key
#        with self.assertRaises(KeyError):
#            self.plot_style.update({'FakeParameters':0})
#    
#
#    def test_set_default_plot_style(self):
#        '''
#        Test the capability of writing pygeostat plot styles into a config file under the user folder
#        '''
#        try:
#            self.plot_style.set_systemdefault()
#        except Exception as ex:
#            self.fail('Unable to set system defaults \n{}'.format(str(ex)))
#
#    def test_get_default_plot_style(self):
#        '''
#        test getting the default plot style
#        '''
#
#        try:
#            self.plot_style.get_systemdefault()
#        except Exception as ex:
#            self.fail('Unable to set system defaults \n{}'.format(str(ex)))
#
#if __name__ == '__main__':
#	subprocess.call([sys.executable, '-m', 'unittest', str(__file__), '-v'])
