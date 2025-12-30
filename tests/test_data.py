#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for DataFile and GridDef functionality in pygeostat.

This module contains comprehensive tests for:
- DataFile class (data loading, manipulation, I/O)
- GridDef class (grid definitions, transformations, random sampling)
"""

import os
import platform

import numpy as np
import pytest

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


# Tests for GridDef.random_indices method
def test_random_indices_1d_basic(grid_from_string):
    """Test random_indices returns valid 1D indices."""
    grid, _ = grid_from_string
    n = 100
    indices = grid.random_indices(dim=1, n=n)

    # Should return exactly n indices
    assert len(indices) == n

    # All indices should be within valid range (0 to grid.count())
    # Note: random_integers is inclusive on both ends, but should be <= count()
    assert all(0 <= idx <= grid.count() for idx in indices)


def test_random_indices_3d_basic(grid_from_string):
    """Test random_indices returns valid 3D indices."""
    grid, _ = grid_from_string
    n = 50
    xind, yind, zind = grid.random_indices(dim=3, n=n)

    # Should return exactly n indices for each dimension
    assert len(xind) == n
    assert len(yind) == n
    assert len(zind) == n

    # All x indices should be within valid range (0 to nx)
    assert all(0 <= x <= grid.nx for x in xind)

    # All y indices should be within valid range (0 to ny)
    assert all(0 <= y <= grid.ny for y in yind)

    # All z indices should be within valid range (0 to nz)
    assert all(0 <= z <= grid.nz for z in zind)


def test_random_indices_reproducible_with_seed(grid_from_string):
    """Test random_indices produces reproducible results with seed."""
    grid, _ = grid_from_string

    # Generate indices with same seed twice
    indices1 = grid.random_indices(dim=1, n=20, seed=42)
    indices2 = grid.random_indices(dim=1, n=20, seed=42)

    # Should get identical results
    np.testing.assert_array_equal(indices1, indices2)


def test_random_indices_3d_reproducible_with_seed(grid_from_string):
    """Test 3D random_indices produces reproducible results with seed."""
    grid, _ = grid_from_string

    # Generate indices with same seed twice
    x1, y1, z1 = grid.random_indices(dim=3, n=20, seed=42)
    x2, y2, z2 = grid.random_indices(dim=3, n=20, seed=42)

    # Should get identical results
    np.testing.assert_array_equal(x1, x2)
    np.testing.assert_array_equal(y1, y2)
    np.testing.assert_array_equal(z1, z2)


def test_random_indices_with_buffer_x(grid_from_string):
    """Test random_indices respects x buffer."""
    grid, _ = grid_from_string
    buffer = 5

    xind, yind, zind = grid.random_indices(dim=3, n=100, buffer_x=buffer, seed=42)

    # All x indices should respect buffer (buffer to nx - buffer)
    assert all(buffer <= x <= grid.nx - buffer for x in xind)


def test_random_indices_with_buffer_y(grid_from_string):
    """Test random_indices respects y buffer."""
    grid, _ = grid_from_string
    buffer = 5

    xind, yind, zind = grid.random_indices(dim=3, n=100, buffer_y=buffer, seed=42)

    # All y indices should respect buffer
    assert all(buffer <= y <= grid.ny - buffer for y in yind)


def test_random_indices_with_buffer_z(grid_from_string):
    """Test random_indices respects z buffer."""
    grid, _ = grid_from_string
    buffer = 5

    xind, yind, zind = grid.random_indices(dim=3, n=100, buffer_z=buffer, seed=42)

    # All z indices should respect buffer
    assert all(buffer <= z <= grid.nz - buffer for z in zind)


def test_random_indices_with_all_buffers(grid_from_string):
    """Test random_indices respects all buffers simultaneously."""
    grid, _ = grid_from_string
    bx, by, bz = 3, 4, 5

    xind, yind, zind = grid.random_indices(dim=3, n=100,
                                           buffer_x=bx,
                                           buffer_y=by,
                                           buffer_z=bz,
                                           seed=42)

    # All indices should respect their respective buffers
    assert all(bx <= x <= grid.nx - bx for x in xind)
    assert all(by <= y <= grid.ny - by for y in yind)
    assert all(bz <= z <= grid.nz - bz for z in zind)


def test_random_indices_buffer_as_tuple(grid_from_string):
    """Test random_indices with buffer as tuple (lower, upper)."""
    grid, _ = grid_from_string
    buffer = (2, 3)  # Different buffer on each side

    xind, yind, zind = grid.random_indices(dim=3, n=100, buffer_x=buffer, seed=42)

    # All x indices should respect asymmetric buffer
    assert all(buffer[0] <= x <= grid.nx - buffer[1] for x in xind)


def test_random_indices_edge_case_n_zero(grid_from_string):
    """Test random_indices with n=0 returns empty array."""
    grid, _ = grid_from_string

    indices = grid.random_indices(dim=1, n=0)
    assert len(indices) == 0


def test_random_indices_edge_case_n_one(grid_from_string):
    """Test random_indices with n=1 returns single index."""
    grid, _ = grid_from_string

    indices = grid.random_indices(dim=1, n=1, seed=42)
    assert len(indices) == 1
    assert 0 <= indices[0] <= grid.count()


def test_random_indices_3d_edge_case_n_one(grid_from_string):
    """Test 3D random_indices with n=1 returns single index per dimension."""
    grid, _ = grid_from_string

    xind, yind, zind = grid.random_indices(dim=3, n=1, seed=42)
    assert len(xind) == 1
    assert len(yind) == 1
    assert len(zind) == 1


def test_random_indices_invalid_seed_raises_error(grid_from_string):
    """Test random_indices raises error for invalid seed."""
    grid, _ = grid_from_string

    with pytest.raises(ValueError):
        grid.random_indices(dim=1, n=10, seed="invalid")


def test_random_indices_invalid_buffer_x_raises_error(grid_from_string):
    """Test random_indices raises error for invalid buffer_x."""
    grid, _ = grid_from_string

    with pytest.raises(ValueError):
        grid.random_indices(dim=3, n=10, buffer_x="invalid")


def test_random_indices_buffer_too_large_raises_error(grid_from_string):
    """Test random_indices raises error when buffer is too large."""
    grid, _ = grid_from_string

    # Buffer larger than grid should raise ValueError
    with pytest.raises(ValueError):
        grid.random_indices(dim=3, n=10, buffer_x=50)  # grid.nx = 40


