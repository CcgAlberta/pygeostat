#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest
import platform  
from pygeostat.data import DataFile, ExampleData
from pygeostat.data.grid_definition import GridDef

# Fixtures
@pytest.fixture
def cleanup_files():
    """Remove test files after each test."""
    yield
    for file in os.listdir('.'):
        _, extension = os.path.splitext(file)
        if extension in ['.vtk', '.h5', '.vtu', '.gsb']:
            os.remove(file)

@pytest.fixture
def grid_from_string():
    """Create GridDef from string for testing."""
    grid = GridDef('''40  0.5  1
                      40  0.5  1
                      40  0.5  1''')
    testfl = ExampleData('3d_grid', griddef=grid)
    return grid, testfl

@pytest.fixture
def grid_from_array():
    """Create GridDef from array for testing."""
    grid = GridDef([40, 0.5, 1, 40, 0.5, 1, 40, 0.5, 1])
    testfl = ExampleData('3d_grid', griddef=grid)
    return grid, testfl

# Helper functions
def compare_grids(griddef1, griddef2):
    """Compare two grid definitions for equality."""
    assert griddef1.nx == griddef2.nx
    assert griddef1.ny == griddef2.ny
    assert griddef1.nz == griddef2.nz
    assert griddef1.xmn == griddef2.xmn
    assert griddef1.ymn == griddef2.ymn
    assert griddef1.zmn == griddef2.zmn
    assert griddef1.xsiz == griddef2.xsiz
    assert griddef1.ysiz == griddef2.ysiz
    assert griddef1.zsiz == griddef2.zsiz


# Tests for DataFile functionality
def test_file_not_found():
    """Test error handling when file not found."""
    with pytest.raises(FileNotFoundError):
        DataFile(flname='notfound')

def test_load_gslib_attributes():
    """Test loading GSLIB file attributes."""
    dat = ExampleData('oilsands')
    assert dat.x.lower() == 'east'
    assert dat.y.lower() == 'north'
    assert dat.z.lower() == 'elevation'

# Skiping the test if not using windows.
## To remove if binary files are compiled for Linux
@pytest.mark.skipif("platform.system() != 'Windows'", 
                   reason="GSB format tests only run on Windows")
def test_write_read_gsb(cleanup_files):
    """Test writing and reading GSB format (Windows only)."""
    from pygeostat.data.iotools import write_gsb, read_gsb
    
    data = ExampleData('point3d_ind_mv')
    
    # Write GSB files
    write_gsb(data, 'test1.gsb')
    data.write_file('test2.gsb')
    
    # Verify files exist
    assert os.path.isfile('test1.gsb')
    assert os.path.isfile('test2.gsb')
    
    # Read GSB files
    data1 = read_gsb('test1.gsb')
    data2 = DataFile('test2.gsb')
    
    # Verify data
    import numpy as np
    assert np.all(data1 == data.data)
    assert np.all(data2.data == data.data)

def test_read_csv():
    """Test reading CSV file."""
    from pygeostat.data.iotools import read_file
    dat = read_file('https://people.sc.fsu.edu/~jburkardt/data/csv/addresses.csv')
    assert dat is not None

# Test to check parameters
## Need to improve this one to replace the test above
#def test_read_csv_with_trim():
#    """Test that read_csv applies the tmin parameter correctly."""
#    import tempfile
#    import pandas as pd
#    from pygeostat.data.iotools import read_csv
#    
#    # Create a temporary CSV file with values that should be trimmed
#    with tempfile.NamedTemporaryFile(suffix='.csv', mode='w+') as temp_file:
#        temp_file.write("A,B,C\n-1,0,1\n-2,3,4\n")
#        temp_file.flush()
#        
#        # Test with tmin=0 (values below 0 should become NaN)
#        result = read_csv(temp_file.name, tmin=0)
#        
#        # Verify the trimming was applied correctly
#        assert pd.isna(result.loc[0, 'A'])
#        assert pd.isna(result.loc[1, 'A'])
#        assert result.loc[0, 'B'] == 0
#        assert result.loc[0, 'C'] == 1

def test_write_read_h5(cleanup_files):
    """Test writing and reading HDF5 format."""
    data = {'East': [1000, 1200], 'North': [2000, 2500]}
    dat = DataFile(data=data)
    dat.write_file('test.h5')
    assert os.path.isfile('test.h5')


def test_write_vtk(cleanup_files):
    """Test writing VTK format."""
    data = {'East': [1000, 1200], 'North': [2000, 2500], 'Elevation': [1, 1]}
    dat = DataFile(data=data, dftype='point')
    dat.write_file('test.vtk')
    # TODO: Ref #100
    assert os.path.isfile('test.vtk.vtu') or os.path.isfile('test.vtk')


def test_set_catdict():
    """Test setting category dictionary."""
    data = {'East': [1000, 1200], 'North': [2000, 2500], 'Elevation': [1, 1], 'Categorical': [0, 1]}
    dat = DataFile(data=data, dftype='point', cat='Categorical')
    dat.setcatdict({0: 'Category 1', 1: 'Category 2'})
    assert dat.catdict == {0: 'Category 1', 1: 'Category 2'}


def test_set_catdict_missing_category():
    """Test error when category dictionary is missing values."""
    data = {
        'East': [1000, 1200, 1500], 
        'North': [2000, 2500, 2600], 
        'Elevation': [1, 1, 2], 
        'Categorical': [0, 1, -99]
    }
    dat = DataFile(data=data, dftype='point', cat='Categorical')

    with pytest.raises(ValueError):
        dat.setcatdict({0: 'Category 1', 1: 'Category 2'})


def test_set_catdict_typeissue():
    """Test error when category values have type issues."""
    data = {
        'East': [1000, 1200, 1500], 
        'North': [2000, 2500, 2600], 
        'Elevation': [1, 1, 2], 
        'Categorical': [0, 1, '-99']
    }
    dat = DataFile(data=data, dftype='point', cat='Categorical')

    with pytest.raises(TypeError):
        dat.setcatdict({0: 'Category 1', 1: 'Category 2', '-99': 'Unknown'})


# Tests for infergriddef functionality
def test_infergriddef_blksize():
    """Test infergriddef with block size."""
    df3d = ExampleData("point3d_ind_mv")
    
    # Where xsiz = ysiz = zsiz, a float can also be provided    
    compare_grids(
        df3d.infergriddef(blksize=75),  
        df3d.infergriddef(blksize=[75, 75, 75])
    )


def test_infergriddef_nblk():
    """Test infergriddef with number of blocks."""
    df3d = ExampleData("point3d_ind_mv")
    
    # Where nx = ny = nz, an int can also be provided
    compare_grids(
        df3d.infergriddef(nblk=60),
        df3d.infergriddef(nblk=[60, 60, 60])
    )


def test_infergriddef_3d_zsiz():
    """Test infergriddef with 3D zsiz variations."""
    df3d = ExampleData("point3d_ind_mv")
    
    # For 3D data, infergriddef() must return a 3D grid definition even if zsiz is given as None or 0 or 1
    compare_grids(
        df3d.infergriddef(blksize=[50, 60, 1]),
        GridDef('''20  135.0  50.0 
                19  1230.0 60.0 
                82  310.5  1.0''')
    )

    compare_grids(
        df3d.infergriddef(blksize=[50, 60, 1]),
        df3d.infergriddef(blksize=[50, 60, None])
    )

    compare_grids(
        df3d.infergriddef(blksize=[50, 60, 1]),
        df3d.infergriddef(blksize=[50, 60, 0])
    )


def test_infergriddef_3d_nz():
    """Test infergriddef with 3D nz variations."""
    df3d = ExampleData("point3d_ind_mv")
    
    # nz given as None or 0 or 1 returns a 2D grid that covers the vertical extent of a 3D dataset
    compare_grids(
        df3d.infergriddef(nblk=[50, 60, 1]),
        GridDef('''50  119.8   19.6 
                60  1209.1  18.2 
                1   350.85  81.7''')
    )

    compare_grids(
        df3d.infergriddef(nblk=[50, 60, 1]),
        df3d.infergriddef(nblk=[50, 60, None])
    )

    compare_grids(
        df3d.infergriddef(nblk=[50, 60, 1]),
        df3d.infergriddef(nblk=[50, 60, 0])
    )


def test_infergriddef_float_to_int():
    """Test infergriddef converting float to int."""
    df3d = ExampleData("point3d_ind_mv")
    
    # If nblk is not integer, infergriddef() must convert it to int and provide correct grid definition 
    _nblk = [12.3, 5.7, 10.11]
    compare_grids(
        df3d.infergriddef(nblk=_nblk),
        df3d.infergriddef(nblk=[int(i) for i in _nblk])
    )


def test_infergriddef_2d_zsiz():
    """Test infergriddef with 2D zsiz variations."""
    df2d = ExampleData("point2d_ind")
    
    # If data is 2-D, zsiz must be provided as None. Otherwise it must raise exception.
    with pytest.raises(ValueError):
        df2d.infergriddef(blksize=[50, 60, 4])

    with pytest.raises(ValueError):
        df2d.infergriddef(blksize=75)

    compare_grids(
        df2d.infergriddef(blksize=[50, 60, None]),  
        GridDef('''24   25.0    50.0  
                 19   1230.0  60.0 
                 1    0.5     1.0''')
    )   


def test_infergriddef_2d_nz():
    """Test infergriddef with 2D nz variations."""
    df2d = ExampleData("point2d_ind")
    
    # If data is 2-D, nz must be provided as None. Otherwise it must raise exception.
    with pytest.raises(ValueError):
        df2d.infergriddef(nblk=[50, 60, 1])       

    with pytest.raises(ValueError):
        df2d.infergriddef(nblk=60)    

    compare_grids(
        df2d.infergriddef(nblk=[50, 60, None]),
        GridDef('''50   11.9    23.8 
                 60   1209.1  18.2 
                 1    0.5     1.0''')
    )

# Tests for GridDef from string
def test_grid_count_string(grid_from_string):
    """Test grid count with string-defined grid."""
    grid, testfl = grid_from_string
    assert testfl.griddef.count() == 64000
    assert grid.count() == 64000


def test_grid_origin_string(grid_from_string):
    """Test grid origin with string-defined grid."""
    grid, testfl = grid_from_string
    assert testfl.griddef.origin() == (0.0, 0.0, 0.0)
    assert grid.origin() == (0.0, 0.0, 0.0)


def test_grid_blockvolume_string(grid_from_string):
    """Test grid block volume with string-defined grid."""
    grid, testfl = grid_from_string
    assert testfl.griddef.block_volume == 1.0
    assert grid.block_volume == 1.0


def test_grid_extents_string(grid_from_string):
    """Test grid extents with string-defined grid."""
    grid, testfl = grid_from_string
    assert testfl.griddef.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
    assert grid.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]


def test_grid_indexcoords_string(grid_from_string):
    """Test grid index coordinates with string-defined grid."""
    grid, testfl = grid_from_string
    assert testfl.griddef.get_index(10.0, 5.0, 15.0) == (22569, True)
    assert grid.get_index(10.0, 5.0, 15.0) == (22569, True)
    # Check for coordinates outside of grid
    assert testfl.griddef.get_index(10.0, 5.0, 50.0) == (-1, False)
    assert grid.get_index(10.0, 5.0, 50.0) == (-1, False)


def test_grid_indicescoords_string(grid_from_string):
    """Test grid indices coordinates with string-defined grid."""
    grid, testfl = grid_from_string
    assert testfl.griddef.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
    assert grid.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
    # Check for coordinates outside of grid
    assert testfl.griddef.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
    assert grid.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)


# Tests for GridDef from array
def test_grid_count_array(grid_from_array):
    """Test grid count with array-defined grid."""
    grid, testfl = grid_from_array
    assert testfl.griddef.count() == 64000
    assert grid.count() == 64000


def test_grid_origin_array(grid_from_array):
    """Test grid origin with array-defined grid."""
    grid, testfl = grid_from_array
    assert testfl.griddef.origin() == (0.0, 0.0, 0.0)
    assert grid.origin() == (0.0, 0.0, 0.0)


def test_grid_blockvolume_array(grid_from_array):
    """Test grid block volume with array-defined grid."""
    grid, testfl = grid_from_array
    assert testfl.griddef.block_volume == 1.0
    assert grid.block_volume == 1.0


def test_grid_extents_array(grid_from_array):
    """Test grid extents with array-defined grid."""
    grid, testfl = grid_from_array
    assert testfl.griddef.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
    assert grid.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]


def test_grid_indexcoords_array(grid_from_array):
    """Test grid index coordinates with array-defined grid."""
    grid, testfl = grid_from_array
    assert testfl.griddef.get_index(10.0, 5.0, 15.0) == (22569, True)
    assert grid.get_index(10.0, 5.0, 15.0) == (22569, True)
    # Check for coordinates outside of grid
    assert testfl.griddef.get_index(10.0, 5.0, 50.0) == (-1, False)
    assert grid.get_index(10.0, 5.0, 50.0) == (-1, False)


def test_grid_indicescoords_array(grid_from_array):
    """Test grid indices coordinates with array-defined grid."""
    grid, testfl = grid_from_array
    assert testfl.griddef.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
    assert grid.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
    # Check for coordinates outside of grid
    assert testfl.griddef.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
    assert grid.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)



#class PygeostatDataFileTest(unittest.TestCase):
#    '''
#     Test suite for pygeostat data file class
#    '''
#
#    def setUp(self):
#        warnings.simplefilter('ignore', category=ImportWarning)
#
#
#    def tearDown(self):
#        for file in os.listdir('.'):
#            _, extention = os.path.splitext(file)
#            if extention in ['.vtk', '.h5', '.vtu', '.gsb']:
#                os.remove(file)
#
#    def test_file_not_found(self):
#        with self.assertRaises(FileNotFoundError):
#            _ = gs.DataFile(flname='notfound')
#
#    def test_load_gslib_attributes(self):
#        dat = gs.ExampleData('oilsands')
#        self.assertEqual(dat.x.lower(), 'east')
#        self.assertEqual(dat.y.lower(), 'north')
#        self.assertEqual(dat.z.lower(), 'elevation')
#
#
#    def test_write_read_gsb(self):
#        data = gs.ExampleData('point3d_ind_mv')
#        try:
#            gs.write_gsb(data, 'test1.gsb')
#        except Exception:
#            self.fail('Unable to write gsb output')
#        
#        try:
#            data.write_file('test2.gsb')
#        except Exception:
#            self.fail('Unable to write gsb output')  
#        
#        import os
#        assert os.path.isfile('test1.gsb'), 'file does not exist'
#        assert os.path.isfile('test1.gsb'), 'file does not exist'
#        
#        try:
#            data1 = gs.read_gsb('test1.gsb')
#        except Exception:
#            self.fail('Unable to read the gsb file')    
#                    
#        try:
#            data2 = gs.DataFile('test2.gsb')
#        except Exception:
#            self.fail('Unable to read the gsb file')      
#        
#        assert np.all(data1 == data.data)
#        assert np.all(data2.data == data.data)
#
#
#    def test_read_csv(self):
#        try:
#            dat= gs.read_file('https://people.sc.fsu.edu/~jburkardt/data/csv/addresses.csv')
#        except Exception:
#            self.fail('Unable to load the CSV file')
#
#    def test_write_read_h5(self):
#        data = {'East': [1000, 1200], 'North': [2000, 2500]}
#        dat = gs.DataFile(data=data)
#
#        try:
#            dat.write_file('test.h5')
#        except Exception:
#            self.fail('Unable to write hdf5 output')
#
#    def test_write_vtk(self):
#        data = {'East': [1000, 1200], 'North': [2000, 2500], 'Elevation': [1,1]}
#        dat = gs.DataFile(data=data, dftype='point')
#
#        try:
#            dat.write_file('test.vtk')
#        except Exception:
#            self.fail('Unable to write vtk output')
#
#    def test_set_catdict(self):
#        data = {'East': [1000, 1200], 'North': [2000, 2500], 'Elevation': [1,1], 'Categorical': [0,1]}
#        dat = gs.DataFile(data=data, dftype='point', cat='Categorical')
#
#        try:
#            dat.setcatdict({0: 'Category 1', 1: 'Category 2'})
#        except Exception:
#            self.fail('Unable to set the categorical dictionary')
#
#    def test_set_catdict_missing_category(self):
#        '''
#        A test to check a missing categorical code in a dictionary raises the correct exception
#        '''
#
#        data = {'East': [1000, 1200, 1500], 'North': [2000, 2500, 2600], 'Elevation': [1,1, 2], 'Categorical': [0,1,-99]}
#        dat = gs.DataFile(data=data, dftype='point', cat='Categorical')
#
#        with self.assertRaises(ValueError):
#            dat.setcatdict({0: 'Category 1', 1: 'Category 2'})
#
#    def test_set_catdict_typeissue(self):
#        '''
#        A test to check that provided non numeric codes can be handled by a proper error message
#        '''
#        data = {'East': [1000, 1200, 1500], 'North': [2000, 2500, 2600], 'Elevation': [1,1, 2], 'Categorical': [0,1,'-99']}
#        dat = gs.DataFile(data=data, dftype='point', cat='Categorical')
#
#        with self.assertRaises(TypeError):
#            dat.setcatdict({0: 'Category 1', 1: 'Category 2', '-99': 'Unknown'})
#
#
#class PygeostatDataFileInfergriddefTest(unittest.TestCase):
#    '''
#     Test suite for pygeostat infergriddef of data file class
#    '''
#    def setUp(self):
#        warnings.simplefilter('ignore', category=ImportWarning)
#
#    def compare(self, griddef1, griddef2):
#        assert griddef1.nx   == griddef2.nx 
#        assert griddef1.ny   == griddef2.ny 
#        assert griddef1.nz   == griddef2.nz 
#        assert griddef1.xmn  == griddef2.xmn 
#        assert griddef1.ymn  == griddef2.ymn 
#        assert griddef1.zmn  == griddef2.zmn
#        assert griddef1.xsiz == griddef2.xsiz 
#        assert griddef1.ysiz == griddef2.ysiz
#        self.assertEqual(griddef1.zsiz,  griddef2.zsiz) 
#           
#    def test_infergriddef(self):
#        df2d = gs.ExampleData("point2d_ind")
#        df3d = gs.ExampleData("point3d_ind_mv")
#        #where xsiz = ysiz = zsiz, a float can also be provided.    
#        self.compare(df3d.infergriddef(blksize = 75),  
#                     df3d.infergriddef(blksize = [75,75,75]))
#
#        #where nx = ny = nz, an int can also be provided.
#        self.compare(df3d.infergriddef(nblk = 60) , 
#                     df3d.infergriddef(nblk = [60,60,60]))
#          
#        #for 3D data, infergriddef() must return a 3D grid definition even if zsiz is given as None or 0 or 1.
#        self.compare(df3d.infergriddef(blksize = [50,60,1]),
#                     gs.GridDef('''20  135.0  50.0 
#                                   19  1230.0 60.0 
#                                   82  310.5  1.0'''))
#
#        self.compare(df3d.infergriddef(blksize = [50,60,1]),
#                     df3d.infergriddef(blksize = [50,60,None]))
#
#        self.compare(df3d.infergriddef(blksize = [50,60,1]),
#                     df3d.infergriddef(blksize = [50,60,0]))             
#                           
#        # nz given as None or 0 or 1 returns a 2D grid that covers the vertical extent of a 3D dataset.
#        self.compare(df3d.infergriddef(nblk = [50,60,1]),
#                     gs.GridDef('''50  119.8   19.6 
#                                   60  1209.1  18.2 
#                                   1   350.85  81.7'''))
#
#        self.compare(df3d.infergriddef(nblk = [50,60,1]),
#                     df3d.infergriddef(nblk = [50,60,None]))
#
#        self.compare(df3d.infergriddef(nblk = [50,60,1]),
#                     df3d.infergriddef(nblk = [50,60,0]))
#
#
#        #if nblk is not integer,infergriddef(), must convert it to int and provide correct grid definition 
#        _nblk = [12.3, 5.7, 10.11]
#        self.compare(df3d.infergriddef(nblk=_nblk),
#                     df3d.infergriddef(nblk=[int(i)  for i in _nblk]))
#
#        # If data is 2-D, zsiz must be provided as None. Otherwise it must raise exception.
#        with self.assertRaises(ValueError):
#            df2d.infergriddef(blksize = [50,60,4])
#
#        with self.assertRaises(ValueError):
#            df2d.infergriddef(blksize = 75)
#
#        self.compare(df2d.infergriddef(blksize = [50,60,None]),  
#                     gs.GridDef('''24   25.0    50.0  
#                                   19   1230.0  60.0 
#                                   1    0.5     1.0'''))   
#
#        #If data is 2-D, nz must be provided as None. Otherwise it must raise exception.
#        with self.assertRaises(ValueError):
#            df2d.infergriddef(nblk = [50,60,1])       
#
#        with self.assertRaises(ValueError):
#            df2d.infergriddef(nblk = 60)    
#
#        self.compare(df2d.infergriddef(nblk = [50,60,None]),  
#                gs.GridDef('''50   11.9    23.8 
#                              60   1209.1  18.2 
#                              1    0.5     1.0'''))
#
#
#
#class PygeostatGridDefFromString(unittest.TestCase):
#    '''
#     Test suite for pygeostat GridDef class with grid string
#    '''
#
#    def setUp(self):
#        warnings.simplefilter('ignore', category=ImportWarning)
#        self.grid = gs.GridDef('''40  0.5  1
#                     40  0.5  1
#                     40  0.5  1''')
#        self.testfl = gs.ExampleData('3d_grid', griddef=self.grid)
#
#    def test_grid_count(self):
#        assert self.testfl.griddef.count() == 64000
#        self.assertEqual(self.grid.count(), 64000)
#
#    def test_grid_origin(self):
#        assert self.testfl.griddef.origin() == (0.0, 0.0, 0.0)
#        assert self.grid.origin() == (0.0, 0.0, 0.0)
#
#    def test_grid_blockvolume(self):
#        assert self.testfl.griddef.block_volume == 1.0
#        assert self.grid.block_volume == 1.0
#
#    def test_grid_extents(self):
#        assert self.testfl.griddef.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
#        assert self.grid.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
#
#    def test_grid_indexcoords(self):
#        assert self.testfl.griddef.get_index(10.0, 5.0, 15.0) == (22569, True)
#        assert self.grid.get_index(10.0, 5.0, 15.0) == (22569, True)
#        # Check for coordinates outside of grid
#        assert self.testfl.griddef.get_index(10.0, 5.0, 50.0) == (-1, False)
#        assert self.grid.get_index(10.0, 5.0, 50.0) == (-1, False)
#
#    def test_grid_indicescoords(self):
#        assert self.testfl.griddef.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
#        assert self.grid.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
#        # Check for coordinates outside of grid
#        assert self.testfl.griddef.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
#        assert self.grid.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
#
#
#class PygeostatGridDefFromArray(unittest.TestCase):
#    '''
#     Test suite for pygeostat GridDef class with grid array
#    '''
#
#    def setUp(self):
#        warnings.simplefilter('ignore', category=ImportWarning)
#        self.grid = gs.GridDef([40,0.5,1,40,0.5,1,40,0.5,1])
#        self.testfl = gs.ExampleData('3d_grid', griddef=self.grid)
#
#    def test_grid_count(self):
#        assert self.testfl.griddef.count() == 64000
#        self.assertEqual(self.grid.count(), 64000)
#
#    def test_grid_origin(self):
#        assert self.testfl.griddef.origin() == (0.0, 0.0, 0.0)
#        assert self.grid.origin() == (0.0, 0.0, 0.0)
#
#    def test_grid_blockvolume(self):
#        assert self.testfl.griddef.block_volume == 1.0
#        assert self.grid.block_volume == 1.0
#
#    def test_grid_extents(self):
#        assert self.testfl.griddef.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
#        assert self.grid.extents() == [(0.0, 40.0), (0.0, 40.0), (0.0, 40.0)]
#
#    def test_grid_indexcoords(self):
#        assert self.testfl.griddef.get_index(10.0, 5.0, 15.0) == (22569, True)
#        assert self.grid.get_index(10.0, 5.0, 15.0) == (22569, True)
#        # Check for coordinates outside of grid
#        assert self.testfl.griddef.get_index(10.0, 5.0, 50.0) == (-1, False)
#        assert self.grid.get_index(10.0, 5.0, 50.0) == (-1, False)
#
#    def test_grid_indicescoords(self):
#        assert self.testfl.griddef.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
#        assert self.grid.get_index3d(10.0, 5.0, 15.0) == (9, 4, 14, True)
#        # Check for coordinates outside of grid
#        assert self.testfl.griddef.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
#        assert self.grid.get_index3d(10.0, 5.0, 50.0) == (-1, -1, -1, False)
#
#
#
#if __name__ == '__main__':
#	subprocess.call([sys.executable, '-m', 'unittest', str(__file__), '-v'])
