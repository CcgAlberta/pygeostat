#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""grid_definition.py: Contains the grid definition class used with pygeostat DataFiles"""

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np


class GridDef:

    """
    Class containing GSLIB grid definition.

    Given either a GSLIB style grid string or an array arranged as [nx, xmn, xsiz, ny, ymn, ysiz,
    nz, zmn, zsiz] initialize the GridDef Class

    """

    def __init__(self, griddef=None, grid_str=None, grid_arr=None, grid_file=None):

        if isinstance(griddef, str) and grid_str is None:
            grid_str = griddef
        elif isinstance(griddef, list) and grid_arr is None:
            grid_arr = griddef

        if grid_file is not None:
            if grid_str is not None:
                print('WARNING: Both a string and an file supplied to grid initialization')
                print('  Assuming the file')
            with open(grid_file) as f:
                grid_str = f.read()

        if grid_str is not None:
            temp_grid_arr = self._parse_grid_string(grid_str)
            if grid_arr is not None:
                if grid_arr != temp_grid_arr:
                    print('WARNING: Both a string and an array supplied to grid initialization')
                    print('  Assuming that the array should be used')
                    temp_grid_arr = grid_arr
        elif grid_arr is not None:
            temp_grid_arr = grid_arr
            if grid_file is not None:
                print('WARNING: Both a string and an file supplied to grid initialization')
                print('  Assuming the array')
        else:
            print('WARNING: No grid definition supplied')
            print('  Assuming zeroes')
            temp_grid_arr = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        self.nx = int(temp_grid_arr[0])
        self.xmn = float(temp_grid_arr[1])
        self.xsiz = float(temp_grid_arr[2])
        self.ny = int(temp_grid_arr[3])
        self.ymn = float(temp_grid_arr[4])
        self.ysiz = float(temp_grid_arr[5])
        self.nz = int(temp_grid_arr[6])
        self.zmn = float(temp_grid_arr[7])
        self.zsiz = float(temp_grid_arr[8])
        self.xlimits = self.extents()[0]
        self.ylimits = self.extents()[1]
        self.zlimits = self.extents()[2]
        self.xlength = (self.xlimits[1] - self.xlimits[0])
        self.ylength = (self.ylimits[1] - self.ylimits[0])
        self.zlength = (self.zlimits[1] - self.zlimits[0])

    def __str__(self):
        """Return the string representation of the current grid definition"""
        return "{} {} {} \n{} {} {} \n{} {} {}".format(self.nx, self.xmn, self.xsiz, self.ny,
                                                       self.ymn, self.ysiz, self.nz, self.zmn,
                                                       self.zsiz)

    def __repr__(self):
        " pretty printing for jupyter notebook workflows"
        return "Pygeostat GridDef:\n{}".format(self)

    def copy(self):
        """ return a copy of this object """
        return GridDef(grid_str=str(self))

    def _parse_grid_string(self, grid_str):
        """Parse the GSLIB style grid string into a list of floats"""
        arr = []
        for line in str.splitlines(grid_str):
            # Check that we have data
            if len(line) > 0:
                # Append the first 3 entries
                arr.append(line.split()[0])
                arr.append(line.split()[1])
                arr.append(line.split()[2])
        assert len(arr) == 9
        return arr

    def count(self):
        """Return the number of blocks in the current grid definition"""
        return self.nx * self.ny * self.nz

    def origin(self):
        """Return the origin of the current grid definition"""
        [(xmin, xmax), (ymin, ymax), (zmin, zmax)] = self.extents()
        return (xmin, ymin, zmin)


    def _get_grid_array(self):
        """ Return the array of grid parameters (nx, xmn, xsiz, ny, ymn, ysiz, nz, zmn, zsiz)  """
        return (self.nx, self.xmn, self.xsiz, self.ny, self.ymn, self.ysiz, self.nz, self.zmn,
                self.zsiz)
    grid_array = property(_get_grid_array)


    def _get_block_volume(self):
        """Return the volume of one block in the current grid definition"""
        return self.xsiz * self.ysiz * self.zsiz

    block_volume = property(_get_block_volume)

    def convert_to_2d(self, orient='xy'):
        """
        Flattens a grid to 2D by default on the xy plane, returning a new 2D GridDef object
        """
        if orient == 'xy':
            return GridDef('%i %f %f\n%i %f %f\n1 0.5 1' % (self.nx, self.xmn, self.xsiz,
                                                            self.ny, self.ymn, self.ysiz))
        elif orient == 'xz':
            return GridDef('%i %f %f\n%i %f %f\n1 0.5 1' % (self.nx, self.xmn, self.xsiz,
                                                            self.nz, self.zmn, self.zsiz))
        elif orient == 'yz':
            return GridDef('%i %f %f\n%i %f %f\n1 0.5 1' % (self.ny, self.ymn, self.ysiz,
                                                            self.nz, self.zmn, self.zsiz))

    def extents(self, orient=None):
        """
        Return the extents of the current grid definition.

        Parameters:
            orient(str): acceptable is 'x','y', or 'z' to give the dimensions along that direction

        Returns:
            various tuples based on what was passed

        """
        xmin = self.xmn - self.xsiz / 2.0
        xmax = xmin + (self.nx) * self.xsiz
        ymin = self.ymn - self.ysiz / 2.0
        ymax = ymin + (self.ny) * self.ysiz
        zmin = self.zmn - self.zsiz / 2.0
        zmax = zmin + (self.nz) * self.zsiz
        if orient is None:
            return [(xmin, xmax), (ymin, ymax), (zmin, zmax)]
        elif orient == 'x':
            return (xmin, xmax)
        elif orient == 'y':
            return (ymin, ymax)
        elif orient == 'z':
            return (zmin, zmax)
        else:
            raise ValueError('Unacceptable orient was provided. x, y and z are the accpetable options')

    def outline_points(self, orient='xy'):
        """
        return the xpts and ypts for plotline an outline of the current grid
        definition in the defined orientation

        Parameters:
            orient(str): The orientation to return the corner points of the grid
                 ``'xy'``, ``'xz'``, ``'yz'`` are the only accepted values

        Returns:
            xpts, ypts (float, float): the corner points of the grid in the specified orientation

        """

        if orient.lower() == 'xy':
            xpts = [self.xlimits[0], self.xlimits[0], self.xlimits[1],
                    self.xlimits[1], self.xlimits[0]]
            ypts = [self.ylimits[0], self.ylimits[1], self.ylimits[1],
                    self.ylimits[0], self.ylimits[0]]
        elif orient.lower() == 'xz':
            xpts = [self.xlimits[0], self.xlimits[0], self.xlimits[1],
                    self.xlimits[1], self.xlimits[0]]
            ypts = [self.zlimits[0], self.zlimits[1], self.zlimits[1],
                    self.zlimits[0], self.zlimits[0]]
        elif orient.lower() == 'yz':
            xpts = [self.ylimits[0], self.ylimits[0], self.ylimits[1],
                    self.ylimits[1], self.ylimits[0]]
            ypts = [self.zlimits[0], self.zlimits[1], self.zlimits[1],
                    self.zlimits[0], self.zlimits[0]]
        else:
            raise ValueError('''`xy`, `xz`, `yz` are the only accepted values for orient''')

        return xpts, ypts

    def index3d_to_index1d(self, ix, iy, iz):
        """
        Will return the 0-indexed 1-Dimensional grid index given the 3 dimensional indices
        as well as a True/False indicator for inside or outside the grid definition.

        Parameters:
            ix (int or numpy.ndarray): x index or n-length array of x indices
            iy (int or numpy.ndarray): y index or n-length array of y indices
            iz (int or numpy.ndarray): z index or n-length array of z indices

        Returns:
            idx (int or numpy.ndarray): 1-d index or n-length array of 1-d indices
            ingrid (bool or numpy.ndarray): in (True) or out (False) of grid, returned if
            ingrid=True

        Examples:

            Calculate a 1d grid index based on a 3d index (integers). Returns the index as an
            integer, as well as a boolean of whether the index is in the grid:

            >>> ix, iy, iz = 2, 4, 0
            >>> idx, ingrid = griddef.index3d_to_index1d(ix, iy, iz)

            Calculate 1d grid indices based on 3d indices (numpy arrays). Returns an array
            of indices, as well as a boolean array of whether the index is in the grid:

            >>> ix, iy, iz = np.arange(0, 5), np.arange(0, 5), np.zeros(5)
            >>> idx, ingrid = griddef.index3d_to_index1d(ix, iy, iz)

        """
        def isall(objtype, objects): return all(isinstance(ob, objtype) for ob in objects)
        if isall(int, [ix, iy, iz]):
            idx = int(ix + iy * self.nx + iz * self.nx * self.ny)
            if idx < self.count() and idx >= 0:
                ingrid = True
            else:
                ingrid = False
        elif isall(np.ndarray, [ix, iy, iz]):
            # Ensure the past indices are integers
            ix = ix.astype(int, casting='safe')
            iy = ix.astype(int, casting='safe')
            iz = ix.astype(int, casting='safe')
            idx = np.add(ix, iy * self.nx, iz * self.nx * self.ny, dtype=int)
            ingrid = np.logical_and(idx >= 0, idx < self.count())
        else:
            raise TypeError(('ix, iy and iz should be the same type, which is either'
                             'int or numpy.ndarray of matching length'))
        return idx, ingrid

    def index1d_to_index3d(self, idx):
        """Will return the 3-dimensional indices given a 1 dimensional index as well as
        a True/False indicator for inside or outside the grid definition.

        Parameters:
            idx (int or np.ndarray): 1 dimensional index or numpy array of indices

        Returns:
            ix (int or np.ndarray): x index or numpy array of x indices
            iy (int or np.ndarray): y index or numpy array of y indices
            iz (int or np.ndarray): z index or numpy array of z indices
            ingrid (bool or np.ndarray): In (True) or Out (False) of grid

        Examples:

            Calculate a 3d grid index based on a 1d index (integers). Returns the 3d index as an
            integer, as well as a boolean of whether the index is in the grid:

            >>> idx = 918
            >>> ix, iy, iz, ingrid = griddef.index1d_to_index3d(idx)

            Calculate 3d grid indices based on 1d indices (numpy array). Returns a arrays
            of 3d indices, as well as a boolean array of whether the indices are in the grid:

            >>> idx = np.array([0, 230, 460, 690, 920])
            >>> ix, iy, iz, ingrid = griddef.index1d_to_index3d(idx)

        """
        if isinstance(idx, (int, np.integer)):
            ix = int((idx % (self.nx * self.ny)) % self.nx)
            iy = int((idx % (self.nx * self.ny)) / self.nx)
            iz = int((idx / (self.nx * self.ny)))
            if idx < self.count() and idx >= 0:
                ingrid = True
            else:
                ingrid = False
        elif isinstance(idx, np.ndarray):
            ix = np.mod(np.mod(idx, self.nx * self.ny), self.nx).astype(int)
            iy = np.floor(np.mod(idx, self.nx * self.ny) / self.nx).astype(int)
            iz = np.floor(idx / (self.nx * self.ny)).astype(int)
            ingrid = np.logical_and(idx >= 0, idx < self.count())
        else:
            raise TypeError(('ix, iy and iz should be the same type, which is either'
                             'int or numpy.ndarray of matching length'))
        return ix, iy, iz, ingrid

    def get_index(self, x, y, z):
        """
        Will return the 0-indexed 1-Dimensional grid index given 3 coordinates as well
        as a True/False indicator for inside or outside the grid definition.

        Uses:
            pygeostat.get_index3d()

        Parameters:
            x (float or numpy.ndarray): x-coordinate value or numpy array of x-coordinate values
            y (float or numpy.ndarray): y-coordinate value or numpy array of y-coordinate values
            z (float or numpy.ndarray): z-coordinate value or numpy array of z-coordinate values

        Returns:
            idx (int or numpy.ndarray): 1-d index of the grid block containing the coordinates
            ingrid (bool or numpy.ndarray): in (True) or out (False) of the grid

        Examples:

            Calculate a 1d grid index based on an input coordinate (floats). Returns the index as an
            integer, as well as a boolean of whether the coordinate is in the grid:

            >>> x, y, z = 30.5, 12.5, 0.5
            >>> idx, ingrid = griddef.gridIndexCoords(x, y, z)

            Calculate 1d grid indices based on input coordinates (numpy arrays). Returns the indices
            as an array, as well as a boolean array of whether the coordinates are in the grid:

            >>> idx = np.array([0, 230, 460, 690, 920])
            >>> ix, iy, iz, ingrid = griddef.index1d_to_index3d(idx)
    """
        # Get the indices
        ix, iy, iz, ingrid = self.get_index3d(x, y, z)
        if isinstance(ix, np.ndarray):
            idx = ix + iy * self.nx + iz * self.nx * self.ny
            # Handle values outside the grid
            idx = np.where(ingrid, idx, -1)
            idx = idx.astype(int)
        else:
            if ingrid:
                idx = int(ix + iy * self.nx + iz * self.nx * self.ny)
            else:
                idx = -1
        return idx, ingrid

    def get_index3d(self, x, y, z):
        """Retruns the 0-indexed 3-Dimensional grid indices given 3 coordinates as well as a
        True/False indicator for inside or outside the grid definition.

        Parameters:
            x (float or numpy.ndarray): x-coordinate value or numpy array of x-coordinate values
            y (float or numpy.ndarray): y-coordinate value or numpy array of y-coordinate values
            z (float or numpy.ndarray): z-coordinate value or numpy array of z-coordinate values

        Returns:
            ix (int or numpy.ndarray): x index of the grid block
            iy (int or numpy.ndarray): y index of the grid block
            iz (int or numpy.ndarray): z index of the grid block
            ingrid (bool or numpy.ndarray): in (True) or out (False) of the grid

        Examples:

            Calculate a 3d grid index based on an input coordinate (floats). Returns the index as
            integers, as well as a boolean of whether the coordinate is in the grid:

            >>> x, y, z = 30.5, 12.5, 0.5
            >>> ix, iy, iz, ingrid = griddef.get_index3d(x, y, z)

            Calculate 3d grid indices based on input coordinates (numpy arrays). Returns the indices
            as arrays, as well as a boolean array of whether the coordinates are in the grid:

            >>> x, y = np.linspace(30.5, 100.5, 5), np.linspace(30.5, 100.5, 5)
            >>> ix, iy, iz, ingrid = griddef.get_index3d(x, y, z)

        """
        if hasattr(x, 'values'):
            x = x.values
        if hasattr(y, 'values'):
            y = y.values
        if hasattr(z, 'values'):
            z = z.values
        if all(isinstance(c, np.ndarray) for c in [x, y, z]):
            # Vectorized ix, iy, iz calcs
            ix = np.ceil((x - self.xmn) / self.xsiz + 0.5)
            iy = np.ceil((y - self.ymn) / self.ysiz + 0.5)
            iz = np.ceil((z - self.zmn) / self.zsiz + 0.5)
            # Check if on the minimum edge
            ix = np.where((ix == 0) & (x == self.xmn - (self.xsiz / 2.0)), 1, ix)
            iy = np.where((iy == 0) & (y == self.ymn - (self.ysiz / 2.0)), 1, iy)
            iz = np.where((iz == 0) & (z == self.zmn - (self.zsiz / 2.0)), 1, iz)
            # Check if within the grid
            ingrid = np.where((ix > 0) & (ix <= self.nx) & (iy > 0) & (
                iy <= self.ny) & (iz > 0) & (iz <= self.nz), True, False)
            # Get the final indices
            ix = ix - 1
            iy = iy - 1
            iz = iz - 1
            # Handle values outside the grid
            ix = np.where(ingrid, ix, -1)
            ix = ix.astype(int)
            iy = np.where(ingrid, iy, -1)
            iy = iy.astype(int)
            iz = np.where(ingrid, iz, -1)
            iz = iz.astype(int)
        else:  # default unvectorized version
            ix = np.ceil((x - self.xmn) / self.xsiz + 0.5)
            iy = np.ceil((y - self.ymn) / self.ysiz + 0.5)
            iz = np.ceil((z - self.zmn) / self.zsiz + 0.5)
            # Check if on minimum edge,
            if ix == 0 and x == self.xmn - (self.xsiz / 2.0):
                ix = 1
            if iy == 0 and y == self.ymn - (self.ysiz / 2.0):
                iy = 1
            if iz == 0 and z == self.zmn - (self.zsiz / 2.0):
                iz = 1
            # Check if within the grid
            if ix > 0 and ix <= self.nx and iy > 0 and iy <= self.ny and iz > 0 and iz <= self.nz:
                ix = int(ix - 1)
                iy = int(iy - 1)
                iz = int(iz - 1)
                ingrid = True
            else:
                ix = -1
                iy = -1
                iz = -1
                ingrid = False
        return ix, iy, iz, ingrid

    def get_coordinates(self, ix=None, iy=None, iz=None, idx=None):
        """Returns the 3 coordinate values based on the passed grid index. If no indices are
        passed, then the coordinates of all grid nodes are returned. Either all 3-D indices
        must all be passed (ix, iy, iz), the 1-D index is passed (idx), or all kwargs are None
        (returns all coordinates).

        Parameters:
            ix (int): x index
            iy (int): y index
            iz (int): z index
            idx (int) : 1-D index

        Returns:
            x (float or numpy.ndarray): x-coordinate value, or values if all grid nodes are returned
            y (float or numpy.ndarray): y-coordinate value, or values if all grid nodes are returned
            z (float or numpy.ndarray): z-coordinate value, or values if all grid nodes are returned

        Note:
            The option to return all grid node coordinates is  memory intensive for > 60 M cell grids.

        Usage:

            Generate a grid definition, and generate the coordinate arrays:

            >>> griddef = gs.GridDef(gridstr="50 0.5 1 \\n50 0.5 1 \\n50 0.5 1")
            >>> x, y, z = griddef.get_coordinates()

            Generate coordinates (floats) corresponding with a specific 3d index:

            >>> x, y, z = griddef.get_coordinates(1, 1, 1)

        """
        if all([test is not None for test in [ix, iy, iz]]):
            if idx is not None:
                raise ValueError('ix, iy and iz must be None if idx is not None (or vice versa)')
            x = ix * self.xsiz + self.xmn
            y = iy * self.ysiz + self.ymn
            z = iz * self.zsiz + self.zmn
        elif idx is not None:
            if any([test is not None for test in [ix, iy, iz]]):
                raise ValueError('ix, iy and iz must be None if idx is not None (or vice versa)')
            ix, iy, iz, _ = self.index1d_to_index3d(idx)
            x = ix * self.xsiz + self.xmn
            y = iy * self.ysiz + self.ymn
            z = iz * self.zsiz + self.zmn
        elif all([test is None for test in [ix, iy, iz, idx]]):
            xmn = self.xmn
            xmx = xmn + self.xsiz * (self.nx - 1)
            x = np.linspace(xmn, xmx, self.nx)
            ymn = self.ymn
            ymx = ymn + self.ysiz * (self.ny - 1)
            y = np.linspace(ymn, ymx, self.ny)
            zmn = self.zmn
            zmx = zmn + self.zsiz * (self.nz - 1)
            z = np.linspace(zmn, zmx, self.nz)
            # Things get a bit weird here. Python does y-fastest so the X and Y are switched here
            y, x, z = np.meshgrid(y, x, z, indexing='xy')
            nshp = self.count()
            x = np.reshape(x, nshp, order='F')
            y = np.reshape(y, nshp, order='F')
            z = np.reshape(z, nshp, order='F')
        else:
            raise ValueError('invalid kwarg combination - refer to the docstring!')
        return x, y, z

    def get_vertical_indices(self, x, y):
        """Returns grid indices corresponding with drilling a vertical
        drill hole intersecting all z blocks on the way down

        Parameters:
            x (float): x-coordinate value
            y (float): y-coordinate value

        Returns:
            indices (dict): grid indices of vertical 'drill hole'"""
        # Get the extents
        [(xmin, xmax), (ymin, ymax), (zmin, zmax)] = self.extents()
        # Check if we are within the bounds of the grid
        if x < xmin or x > xmax or y < ymin or y > ymax:
            return []
        z = np.linspace(zmin + self.zsiz * 0.5, zmax - self.zsiz * 0.5, num=self.nz)
        x, y = x * np.ones(z.shape), y * np.ones(z.shape),
        idx, _ = self.get_index(x, y, z)
        return idx

    def get_slice_coordinate(self, orient, index):
        """Returns the real coordinate for a slice given the index in the grid and orientation
        NOTE: assumes 0-indexed slice coordinates are passed.

        Parameters:
            orient (str) : orientation of the current slice ('xy', 'yz', 'xz')
            index (int) : index of the current slice

        Returns:
            cord (float) : coordinate of the current slice

        Usage:
        >>> griddef.get_slice_coordinate('xy',10)

        """
        if orient == 'xy':  # want the z coordinate
            cord = self.zmn + self.zsiz * index
        elif orient == 'yz':  # want the x coordinate
            cord = self.xmn + self.xsiz * index
        elif orient == 'xz':  # want the y coordinate
            cord = self.ymn + self.ysiz * index

        return cord

    def get_slice_index(self, orient, coordinate):
        """Returns the index in the grid along the correct dimension for the specificed slice
        i.e. the z index for an 'xy' slice

        >>> griddef.get_slice_index('xy',700.2)

        Parameters:
            orient (str) : orientation of the current slice ('xy', 'yz', 'xz')
            coordinate (float) : index of the current slice

        Returns:
            index (int) : index to the specified slice

        """
        if orient == 'xy':  # want the z coordinate
            ix, iy, iz, infl = self.get_index3d(self.xmn, self.ymn, coordinate)
            index = iz
        elif orient == 'yz':  # want the x coordinate
            ix, iy, iz, infl = self.get_index3d(coordinate, self.ymn, self.zmn)
            index = ix
        elif orient == 'xz':  # want the y coordinate
            ix, iy, iz, infl = self.get_index3d(self.xmn, coordinate, self.zmn)
            index = iy
        return index

    def change_blocknum(self, nx_new, ny_new, nz_new, return_copy=False):
        """
        Function to change the blocksize of the current grid while retaining the original bounding
        box.

        Useful if attempting to work at a coarse grid (for speed) prior to obtaining a final
        estimate at the original resolution.

        Parameters:
            nx_new (int): New number of blocks in X direction
            ny_new (int): New number of blocks in Y direction
            nz_new (int): New number of blocks in Z direction

        Example:
            Define a grid:

            >>> import pygeostat as gs
            >>> grid = gs.GridDef(grid_str="100 5 10  \\n100 5 10  \\n100 5 5")

            Change the dimensions of the grid:

            >>> grid.changedim(50,50,50)
            >>> print(grid.nx,grid.xmn,grid.xsiz)
            >>> print(grid.ny,grid.ymn,grid.ysiz)
            >>> print(grid.nz,grid.zmn,grid.zsiz)

            Use the changed resolution grid in a parameter file:

            >>> parstr = "TestParFile \\n{grd}"
            >>> prog = gs.Program(program='./text.exe',parstr=parstr.format(grd=str(grid)))

        """
        bounds = np.array(self.extents())
        span = np.array([bounds[0, 1] - bounds[0, 0],
                         bounds[1, 1] - bounds[1, 0],
                         bounds[2, 1] - bounds[2, 0]])
        newsize = span / np.array([nx_new, ny_new, nz_new])
        newmin = bounds[:, 0] + 0.5 * newsize
        # assign the new dimensions to the class
        nx = int(nx_new)
        ny = int(ny_new)
        nz = int(nz_new)
        xmn = newmin[0]
        ymn = newmin[1]
        zmn = newmin[2]
        xsiz = newsize[0]
        ysiz = newsize[1]
        zsiz = newsize[2]

        if return_copy is False:
            self.__init__(grid_arr=[nx, xmn, xsiz, ny, ymn, ysiz, nz, zmn, zsiz])
        else:
            from . import GridDef as GridDef
            return GridDef(grid_arr=[nx, xmn, xsiz, ny, ymn, ysiz, nz, zmn, zsiz])

    def change_blocksize(self, xsiz_new, ysiz_new, zsiz_new, return_copy=False):
        """
        Function to change the size of individual blocks in the grid.
        Finds the new number of blocks
        given the target sizes in each direction.

        Parameters:
            xsiz_new (float): New size of blocks in X
            ysiz_new (float): New size of blocks in Y
            zsiz_new (float): New size of blocks in Z
            return_copy (bool): if True will return a copy instead of modifying self

        """
        bounds = np.array(self.extents())
        span = np.array([bounds[0, 1] - bounds[0, 0],
                         bounds[1, 1] - bounds[1, 0],
                         bounds[2, 1] - bounds[2, 0]])
        newblocks = np.ceil(span / np.array([xsiz_new, ysiz_new, zsiz_new]))
        newmin = bounds[:, 0] + 0.5 * np.array([xsiz_new, ysiz_new, zsiz_new])
        # assign the new dimensions to the class
        nx = int(newblocks[0])
        ny = int(newblocks[1])
        nz = int(newblocks[2])
        xmn = newmin[0]
        ymn = newmin[1]
        zmn = newmin[2]
        xsiz = xsiz_new
        ysiz = ysiz_new
        zsiz = zsiz_new

        if return_copy is False:
            self.__init__(grid_arr=[nx, xmn, xsiz, ny, ymn, ysiz, nz, zmn, zsiz])
        else:
            from . import GridDef as GridDef
            return GridDef(grid_arr=[nx, xmn, xsiz, ny, ymn, ysiz, nz, zmn, zsiz])

    def generate_grid_points(self):
        """
        Returns a (ncell x 3) set of real location grid points that conform to the standard gslib
        orderings, starting at the bottom of the model and increasing upwards, x-fastest, then y,
        then z.

        Note:
            Memory intensive for > 60 M cell grids. This function likely rarely needs to be used,
            and does not explicitly account for a mask in construction. Can be used to get the set
            of GSLIB evaluation locations from a griddef, for example.

        Usage:
            Generate a grid definition, and generate the array of locations:

            >>> griddef = gs.GridDef(grid_str="50 0.5 1 \\n50 0.5 1 \\n50 0.5 1")
            >>> gridxyz = griddef.generate_grid_points()

        """
        xmin = self.xmn
        xmax = xmin + self.xsiz * (self.nx - 1)
        xrng = np.linspace(xmin, xmax, self.nx)

        ymin = self.ymn
        ymax = ymin + self.ysiz * (self.ny - 1)
        yrng = np.linspace(ymin, ymax, self.ny)

        zmin = self.zmn
        zmax = zmin + self.zsiz * (self.nz - 1)
        zrng = np.linspace(zmin, zmax, self.nz)

        # Python does y-fastest so the X and Y are switched here
        y1, x1, z1 = np.meshgrid(yrng, xrng, zrng, indexing='xy')
        nshp = self.count()
        x = np.reshape(x1, nshp, order='F')
        y = np.reshape(y1, nshp, order='F')
        z = np.reshape(z1, nshp, order='F')
        xyz = np.c_[x, y, z]
        return xyz

    def pad_grid(self, nx_pad, ny_pad, nz_pad, return_copy=False):
        """
        Pad the grid on either side in all directions by the number of cells specified on input

        Parameters:
            nx_pad (int or tuple): number of cells in x direction to add to the grid on each side
            ny_pad (int or tuple): number of cells in y direction to add to the grid on each side
            nz_pad (int or tuple): number of cells in z direction to add to the grid on each side
            return_copy (bool): return copy or reinitialize self

        Examples:

            Generate a grid definition:

            >>> griddef = gs.GridDef(gridstr="50 0.5 1 \\n50 0.5 1 \\n1 0.5 1")

            Symmetrically pad the grid cells in the x and y directions

            >>> griddef2 = griddef2.padgrid(10, 10, 0, return_copy=True)

            Asymmetrically pad the grid with cells in the x and y directions

            >>> griddef2.padgrid((6, -5), 5, 0)

        """
        if isinstance(nx_pad, int):
            xmn = self.xmn - nx_pad * self.xsiz
            nx = self.nx + 2 * nx_pad
        elif isinstance(nx_pad, tuple):
            xmn = self.xmn - int(nx_pad[0]) * self.xsiz
            nx = self.nx + int(nx_pad[0] + nx_pad[1])
        else:
            raise ValueError('nx_pad must be either an integer or tuple')

        if isinstance(ny_pad, int):
            ymn = self.ymn - ny_pad * self.ysiz
            ny = self.ny + 2 * ny_pad
        elif isinstance(nx_pad, tuple):
            ymn = self.ymn - int(ny_pad[0]) * self.ysiz
            ny = self.ny + int(ny_pad[0] + ny_pad[1])
        else:
            raise ValueError('nx_pad must be either an integer or tuple')

        if isinstance(nz_pad, int):
            zmn = self.zmn - nz_pad * self.zsiz
            nz = self.nz + 2 * nz_pad
        elif isinstance(nx_pad, tuple):
            zmn = self.zmn - int(nz_pad[0]) * self.zsiz
            nz = self.nz + int(nz_pad[0] + nz_pad[1])
        else:
            raise ValueError('nx_pad must be either an integer or tuple')

        # reinitialize this object with the new parameters

        if return_copy is False:
            self.__init__(grid_arr=[nx, xmn, self.xsiz, ny, ymn, self.ysiz, nz, zmn, self.zsiz])
        else:
            from . import GridDef as GridDef
            return GridDef(grid_arr=[nx, xmn, self.xsiz, ny, ymn, self.ysiz, nz, zmn, self.zsiz])

    def random_points(self, n=100):
        """
        Generate random points from within the grid
        """
        x = self.xmn + np.random.rand(n) * self.nx * self.xsiz
        y = self.ymn + np.random.rand(n) * self.ny * self.ysiz
        if self.nz > 1:
            z = self.zmn + np.random.rand(n) * self.nz * self.zsiz
            return np.c_[x, y, z]
        else:
            return np.c_[x, y]

    def random_indices(self, dim=3, n=20, seed=None, buffer_x=None, buffer_y=None, buffer_z=None):
        """
        Generate a list of random indices from within the grid

        Parameters:
            dim (int): 1 returns a 1D index, 3 returns x, y, and z indexes.
            n (int): The number of indices to return for each dimension
            seed (int): If provided will initialized the random number generator with the seed.
                If not provided then you will get different values every time you run this function
            buffer_x (int): you can set a buffer if you don't want any random indices near
                the edge of the x border of the grid, Note: only works for dim=3
            buffer_y (int): you can set a buffer if you don't want any random indices near
                the edge of the y border of the grid, Note: only works for dim=3
            buffer_z (int): you can set a buffer if you don't want any random indices near
                the edge of the z border of the grid, Note: only works for dim=3

        Returns:
            "dim = 1"
            indice (list): a 1D indice of size `n`
            "dim = 3"
            xind (list): a list of indices in x dimension of size `n`
            yind (list): a list of indices in y dimension of size `n`
            zind (list): a list of indices in z dimension of size `n`
        """
        from .. utility.logging import printerr

        if seed:
            if isinstance(seed, int):
                np.random.seed(seed)
            else:
                raise ValueError('seed must be an integer')

        if buffer_x:
            if isinstance(buffer_x, int):
                bx = (buffer_x, buffer_x)
            elif isinstance(buffer_x, tuple):
                if isinstance(buffer_x[0], int) and isinstance(buffer_x[1], int):
                    bx = buffer_x
            else:
                raise ValueError('buffer_x is not an integer or a tuple of integers', buffer_x)

        if buffer_y:
            if isinstance(buffer_y, int):
                by = (buffer_y, buffer_y)
            elif isinstance(buffer_y, tuple):
                if isinstance(buffer_y[0], int) and isinstance(buffer_y[1], int):
                    by = buffer_y
            else:
                raise ValueError('buffer_y is not an integer or a tuple of integers', buffer_y)

        if buffer_z:
            if isinstance(buffer_z, int):
                bz = (buffer_z, buffer_z)
            elif isinstance(buffer_z, tuple):
                if isinstance(buffer_z[0], int) and isinstance(buffer_z[1], int):
                    bz = buffer_z
            else:
                raise ValueError('buffer_z is not an integer or a tuple of integers', buffer_z)

        if dim != 1:
            if dim != 3:
                printerr('Only 1 or 3 Dimensions allowed. Assuming you wanted 3 dimension',
                         errtype='warning')
                print('dim entered = ', dim)
            if buffer_x:
                xrange = (0 + bx[0], self.nx - bx[1])
                if xrange[0] > xrange[1] or xrange[0] == xrange[1]:
                    raise ValueError('buffer_x is to large for the grid causing the upper range'
                                     ' to be less then or equal to the lower range', xrange)
            else:
                xrange = (0, self.nx)
            if buffer_y:
                yrange = (0 + by[0], self.ny - by[1])
                if yrange[0] > yrange[1] or yrange[0] == yrange[1]:
                    raise ValueError('buffer_y is to large for the grid causing the upper range'
                                     '  to be less then or equal to the lower range', yrange)
            else:
                yrange = (0, self.ny)
            if buffer_z:
                zrange = (0 + bz[0], self.nz - bz[1])
                if zrange[0] > zrange[1] or zrange[0] == zrange[1]:
                    raise ValueError('buffer_z is to large for the grid causing the upper range'
                                     ' to be less then or equal to the lower range', zrange)
            else:
                zrange = (0, self.nz)
            rand_xind = np.random.random_integers(xrange[0], xrange[1], size=n)
            rand_yind = np.random.random_integers(yrange[0], yrange[1], size=n)
            rand_zind = np.random.random_integers(zrange[0], zrange[1], size=n)
        else:
            rand_indices = np.random.random_integers(0, self.count(), size=n)

        if dim == 1:
            return rand_indices
        else:
            return rand_xind, rand_yind, rand_zind
