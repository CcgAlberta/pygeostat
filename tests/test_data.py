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


# class PygeostatDataFileTest(unittest.TestCase):
#     '''
#      Test suite for pygeostat data file class
#     '''

#     def setUp(self):
#         warnings.simplefilter('ignore', category=ImportWarning)


#     def tearDown(self):
#         for file in os.listdir('.'):
#             _, extention = os.path.splitext(file)
#             if extention in ['.vtk', '.h5', '.vtu']:
#                 os.remove(file)

#     def test_file_not_found(self):
#         with self.assertRaises(FileNotFoundError):
#             _ = gs.DataFile(flname='notfound')

#     def test_load_gslib_attributes(self):
#         dat = gs.ExampleData('oilsands')
#         self.assertEqual(dat.x.lower(), 'east')
#         self.assertEqual(dat.y.lower(), 'north')
#         self.assertEqual(dat.z.lower(), 'elevation')

#     def test_read_csv(self):
#         try:
#             dat= gs.read_file('https://people.sc.fsu.edu/~jburkardt/data/csv/addresses.csv')
#         except Exception:
#             self.fail('Unable to load the CSV file')

#     def test_write_read_h5(self):
#         data = {'East': [1000, 1200], 'North': [2000, 2500]}
#         dat = gs.DataFile(data=data)

#         try:
#             dat.write_file('test.h5')
#         except Exception:
#             self.fail('Unable to write hdf5 output')

#     def test_write_vtk(self):
#         data = {'East': [1000, 1200], 'North': [2000, 2500], 'Elevation': [1,1]}
#         dat = gs.DataFile(data=data, dftype='point')

#         try:
#             dat.write_file('test.vtk')
#         except Exception:
#             self.fail('Unable to write vtk output')

#     def test_set_catdict(self):
#         data = {'East': [1000, 1200], 'North': [2000, 2500], 'Elevation': [1,1], 'Categorical': [0,1]}
#         dat = gs.DataFile(data=data, dftype='point', cat='Categorical')

#         try:
#             dat.setcatdict({0: 'Category 1', 1: 'Category 2'})
#         except Exception:
#             self.fail('Unable to set the categorical dictionary')

#     def test_set_catdict_missing_category(self):
#         '''
#         A test to check a missing categorical code in a dictionary raises the correct exception
#         '''

#         data = {'East': [1000, 1200, 1500], 'North': [2000, 2500, 2600], 'Elevation': [1,1, 2], 'Categorical': [0,1,-99]}
#         dat = gs.DataFile(data=data, dftype='point', cat='Categorical')

#         with self.assertRaises(ValueError):
#             dat.setcatdict({0: 'Category 1', 1: 'Category 2'})

#     def test_set_catdict_typeissue(self):
#         '''
#         A test to check that provided non numeric codes can be handled by a proper error message
#         '''
#         data = {'East': [1000, 1200, 1500], 'North': [2000, 2500, 2600], 'Elevation': [1,1, 2], 'Categorical': [0,1,'-99']}
#         dat = gs.DataFile(data=data, dftype='point', cat='Categorical')

#         with self.assertRaises(TypeError):
#             dat.setcatdict({0: 'Category 1', 1: 'Category 2', '-99': 'Unknown'})


class PygeostatDataFile_infergriddef_Test(unittest.TestCase):
    '''
     Test suite for pygeostat infergriddef of data file class
    '''
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)

    def compare(self, griddef1, griddef2):
        assert griddef1.nx   == griddef2.nx 
        assert griddef1.ny   == griddef2.ny 
        assert griddef1.nz   == griddef2.nz 
        assert griddef1.xmn  == griddef2.xmn 
        assert griddef1.ymn  == griddef2.ymn 
        assert griddef1.zmn  == griddef2.zmn
        assert griddef1.xsiz == griddef2.xsiz 
        assert griddef1.ysiz == griddef2.ysiz
        self.assertEqual(griddef1.zsiz,  griddef2.zsiz) 
           
    def test_infergriddef(self):
        df2d = gs.ExampleData("point2d_ind")
        df3d = gs.ExampleData("point3d_ind_mv")
        #where xsiz = ysiz = zsiz, a float can also be provided.    
        self.compare(df3d.infergriddef(blksize = 75),  
                     df3d.infergriddef(blksize = [75,75,75]))

        #where nx = ny = nz, an int can also be provided.
        self.compare(df3d.infergriddef(nblk = 60) , 
                     df3d.infergriddef(nblk = [60,60,60]))
          
        #for 3D data, infergriddef() must return a 3D grid definition even if zsiz is given as None or 0 or 1.
        self.compare(df3d.infergriddef(blksize = [50,60,1]),
                     gs.GridDef('''20  135.0  50.0 
                                   19  1230.0 60.0 
                                   82  310.5  1.0'''))

        self.compare(df3d.infergriddef(blksize = [50,60,1]),
                     df3d.infergriddef(blksize = [50,60,None]))

        self.compare(df3d.infergriddef(blksize = [50,60,1]),
                     df3d.infergriddef(blksize = [50,60,0]))             
                           
        # nz given as None or 0 or 1 returns a 2d grid that covers the vertical extent of a 3D dataset.
        self.compare(df3d.infergriddef(nblk = [50,60,1]),
                     gs.GridDef('''50  119.8   19.6 
                                   60  1209.1  18.2 
                                   1   350.85  81.7'''))

        self.compare(df3d.infergriddef(nblk = [50,60,1]),
                     df3d.infergriddef(nblk = [50,60,None]))

        self.compare(df3d.infergriddef(nblk = [50,60,1]),
                     df3d.infergriddef(nblk = [50,60,0]))


        #if nblk is not integer,infergriddef(), must convert it to int and provide correct grid definition 
        _nblk = [12.3, 5.7, 10.11]
        self.compare(df3d.infergriddef(nblk=_nblk),
                     df3d.infergriddef(nblk=[int(i)  for i in _nblk]))

        # If the grid is 2-D, zsiz must be provided as None. Otherwise it must raise exception.
        with self.assertRaises(ValueError):
            df2d.infergriddef(blksize = [50,60,4])

        with self.assertRaises(ValueError):
            df2d.infergriddef(blksize = 75)

        self.compare(df2d.infergriddef(blksize = [50,60,None]),  
                     gs.GridDef('''24   25.0    50.0  
                                   19   1230.0  60.0 
                                   1    0.5     1.0'''))   

        #If the grid is 2-D, nz must be provided as None. Otherwise it must raise exception.
        with self.assertRaises(ValueError):
            df2d.infergriddef(nblk = [50,60,1])       

        with self.assertRaises(ValueError):
            df2d.infergriddef(nblk = 60)    

        self.compare(df2d.infergriddef(nblk = [50,60,None]),  
                gs.GridDef('''50   11.9    23.8 
                              60   1209.1  18.2 
                              1    0.5     1.0'''))



# class PygeostatGridDefFromString(unittest.TestCase):
#     '''
#      Test suite for pygeostat GridDef class with grid string
#     '''

#     def setUp(self):
#         warnings.simplefilter('ignore', category=ImportWarning)
#         self.grid = gs.GridDef('''40  0.5  1
#                      40  0.5  1
#                      40  0.5  1''')
#         self.testfl = gs.ExampleData('3d_grid', griddef=self.grid)

#     def test_grid_count(self):
#         assert self.testfl.griddef.count() == 64000
#         self.assertEqual(self.grid.count(), 64000)

#     def test_grid_origin(self):
#         assert self.testfl.griddef.origin() == (0.0, 0.0, 0.0)
#         assert self.grid.origin() == (0.0, 0.0, 0.0)

#     def test_grid_blockvolume(self):
#         assert self.testfl.griddef.block_volume == 1.0
#         assert self.grid.block_volume == 1.0

#     def test_grid_extents(self):
#         assert self.testfl.griddef.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
#         assert self.grid.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]

#     def test_grid_indexcoords(self):
#         assert self.testfl.griddef.get_index(10.0, 5.0, 15.0) == (22569, True)
#         assert self.grid.get_index(10.0, 5.0, 15.0) == (22569, True)
#         # Check for coordinates outside of grid
#         assert self.testfl.griddef.get_index(10.0, 5.0, 50.0) == (-1, False)
#         assert self.grid.get_index(10.0, 5.0, 50.0) == (-1, False)

#     def test_grid_indicescoords(self):
#         assert self.testfl.griddef.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
#         assert self.grid.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
#         # Check for coordinates outside of grid
#         assert self.testfl.griddef.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
#         assert self.grid.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)


# class PygeostatGridDefFromArray(unittest.TestCase):
#     '''
#      Test suite for pygeostat GridDef class with grid array
#     '''

#     def setUp(self):
#         warnings.simplefilter('ignore', category=ImportWarning)
#         self.grid = gs.GridDef([40,0.5,1,40,0.5,1,40,0.5,1])
#         self.testfl = gs.ExampleData('3d_grid', griddef=self.grid)

#     def test_grid_count(self):
#         assert self.testfl.griddef.count() == 64000
#         self.assertEqual(self.grid.count(), 64000)

#     def test_grid_origin(self):
#         assert self.testfl.griddef.origin() == (0.0, 0.0, 0.0)
#         assert self.grid.origin() == (0.0, 0.0, 0.0)

#     def test_grid_blockvolume(self):
#         assert self.testfl.griddef.block_volume == 1.0
#         assert self.grid.block_volume == 1.0

#     def test_grid_extents(self):
#         assert self.testfl.griddef.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
#         assert self.grid.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]

#     def test_grid_indexcoords(self):
#         assert self.testfl.griddef.get_index(10.0, 5.0, 15.0) == (22569, True)
#         assert self.grid.get_index(10.0, 5.0, 15.0) == (22569, True)
#         # Check for coordinates outside of grid
#         assert self.testfl.griddef.get_index(10.0, 5.0, 50.0) == (-1, False)
#         assert self.grid.get_index(10.0, 5.0, 50.0) == (-1, False)

#     def test_grid_indicescoords(self):
#         assert self.testfl.griddef.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
#         assert self.grid.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
#         # Check for coordinates outside of grid
#         assert self.testfl.griddef.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
#         assert self.grid.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)



if __name__ == '__main__':
	subprocess.call([sys.executable, '-m', 'unittest', str(__file__), '-v'])