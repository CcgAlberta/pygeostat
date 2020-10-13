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


class PygeostatDataFileTest(unittest.TestCase):
    '''
     Test suite for pygeostat data file class
    '''

    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)


    def tearDown(self):
        for file in os.listdir('.'):
            _, extention = os.path.splitext(file)
            if extention in ['.vtk', '.h5', '.vtu']:
                os.remove(file)

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            _ = gs.DataFile(flname='notfound')

    def test_load_gslib_attributes(self):
        dat = gs.ExampleData('oilsands')
        self.assertEqual(dat.x.lower(), 'east')
        self.assertEqual(dat.y.lower(), 'north')
        self.assertEqual(dat.z.lower(), 'elevation')

    def test_read_csv(self):
        try:
            dat= gs.read_file('https://people.sc.fsu.edu/~jburkardt/data/csv/addresses.csv')
        except Exception:
            self.fail('Unable to load the CSV file')

    def test_write_read_h5(self):
        data = {'East': [1000, 1200], 'North': [2000, 2500]}
        dat = gs.DataFile(data=data)

        try:
            dat.write_file('test.h5')
        except Exception:
            self.fail('Unable to write hdf5 output')

    def test_write_vtk(self):
        data = {'East': [1000, 1200], 'North': [2000, 2500], 'Elevation': [1,1]}
        dat = gs.DataFile(data=data, dftype='point')

        try:
            dat.write_file('test.vtk')
        except Exception:
            self.fail('Unable to write vtk output')

class PygeostatGridDefFromString(unittest.TestCase):
    '''
     Test suite for pygeostat GridDef class with grid string
    '''

    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.grid = gs.GridDef('''40  0.5  1
                     40  0.5  1
                     40  0.5  1''')
        self.testfl = gs.ExampleData('3d_grid', griddef=self.grid)

    def test_grid_count(self):
        assert self.testfl.griddef.count() == 64000
        self.assertEqual(self.grid.count(), 64000)

    def test_grid_origin(self):
        assert self.testfl.griddef.origin() == (0.0, 0.0, 0.0)
        assert self.grid.origin() == (0.0, 0.0, 0.0)

    def test_grid_blockvolume(self):
        assert self.testfl.griddef.block_volume == 1.0
        assert self.grid.block_volume == 1.0

    def test_grid_extents(self):
        assert self.testfl.griddef.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
        assert self.grid.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]

    def test_grid_indexcoords(self):
        assert self.testfl.griddef.get_index(10.0, 5.0, 15.0) == (22569, True)
        assert self.grid.get_index(10.0, 5.0, 15.0) == (22569, True)
        # Check for coordinates outside of grid
        assert self.testfl.griddef.get_index(10.0, 5.0, 50.0) == (-1, False)
        assert self.grid.get_index(10.0, 5.0, 50.0) == (-1, False)

    def test_grid_indicescoords(self):
        assert self.testfl.griddef.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
        assert self.grid.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
        # Check for coordinates outside of grid
        assert self.testfl.griddef.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
        assert self.grid.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)


class PygeostatGridDefFromArrayg(unittest.TestCase):
    '''
     Test suite for pygeostat GridDef class with grid array
    '''

    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.grid = gs.GridDef([40,0.5,1,40,0.5,1,40,0.5,1])
        self.testfl = gs.ExampleData('3d_grid', griddef=self.grid)

    def test_grid_count(self):
        assert self.testfl.griddef.count() == 64000
        self.assertEqual(self.grid.count(), 64000)

    def test_grid_origin(self):
        assert self.testfl.griddef.origin() == (0.0, 0.0, 0.0)
        assert self.grid.origin() == (0.0, 0.0, 0.0)

    def test_grid_blockvolume(self):
        assert self.testfl.griddef.block_volume == 1.0
        assert self.grid.block_volume == 1.0

    def test_grid_extents(self):
        assert self.testfl.griddef.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
        assert self.grid.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]

    def test_grid_indexcoords(self):
        assert self.testfl.griddef.get_index(10.0, 5.0, 15.0) == (22569, True)
        assert self.grid.get_index(10.0, 5.0, 15.0) == (22569, True)
        # Check for coordinates outside of grid
        assert self.testfl.griddef.get_index(10.0, 5.0, 50.0) == (-1, False)
        assert self.grid.get_index(10.0, 5.0, 50.0) == (-1, False)

    def test_grid_indicescoords(self):
        assert self.testfl.griddef.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
        assert self.grid.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
        # Check for coordinates outside of grid
        assert self.testfl.griddef.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
        assert self.grid.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)



if __name__ == '__main__':
	subprocess.call([sys.executable, '-m', 'unittest', str(__file__), '-v'])