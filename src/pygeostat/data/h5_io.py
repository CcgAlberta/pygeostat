#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
h5_io.py: Contains input/output functions for using HDF5 files within pygeostat
'''
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import textwrap

from ..utility.logging import printerr


def _fixh5path(h5path):
    """Fix an h5path for use in h5py file objects"""
    if h5path in [None, '', ' ', '/']:
        h5path = '/'
    if h5path[0] != '/':
        h5path = '/' + h5path
    if h5path[-1] != '/':
        h5path = h5path + '/'

    return h5path


def h5_combine_data(flname, h5paths, datasets=None):
    """
    Combine data into one DataFrame from multiple paths in a HDF5 file.

    Parameters:
        flname (str): Path of the HDF5 you wish to read from
        h5paths (list): A list of h5paths to combine. Forward slash (/) delimited path through the
            group hierarchy you wish to place the dataset(s) specified by the argument ``datasets``
            into. The dataset name cannot be passed using this argument, it is interpreted as a
            group name. A value of ``None`` places the dataset into the root directory of the
            HDF5 file.
        datasets (list of lists): If only a specific set of datasets from each path are desired then
            pass a list of lists of equal length as the h5paths list. An empty list within the list
            will cause all datasets in the corresponding path to be readin.

    Returns:
        DataFrame

    Example:

    >>> flname = 'drilldata.h5'
    ... h5paths = ['/Orig_data/series4870/', 'NS/Declus/series4870/']
    ... datasets = [['LOCATIONX', 'LOCATIONY', 'LOCATIONZ'], []]
    ... data = gs.h5_combine_data(flname, h5paths, datasets=datasets)
    """
    import h5py
    # Handle some default parameters
    if isinstance(h5paths, list):
        if len(h5paths) > 1:
            for path in h5paths:
                if path in [None, True, '', ' ', '/']:
                    path = '/'
        else:
            printerr("Multiple paths were not passed. If you are not combining data across paths"
                     "use `gs.read_h5()`", errtype='error')
    else:
        printerr("Multiple paths were not passed. If you are not combining data across paths use"
                 " `gs.read_h5()`", errtype='error')
    # Sanity checks
    if not h5py.is_hdf5(flname):
        raise IOError("The passed file path is not a valid HDF5 file.")

    for idx, path in enumerate(h5paths):
        if len(datasets[idx]) > 0:
            if idx == 0:
                data = read_h5(flname=flname, h5path=path,
                               datasets=datasets[idx])
            else:
                data2 = read_h5(flname=flname, h5path=path,
                                datasets=datasets[idx])
                data = pd.concat([data, data2], axis=1)
        else:
            if idx == 0:
                data = read_h5(flname=flname, h5path=path)
            else:
                data2 = read_h5(flname=flname, h5path=path)
                data = pd.concat([data, data2], axis=1)

    return data


def write_h5(data, flname, h5path=None, datasets=None, dtype=None, gridstr=None,
               trim_variable=None, var_min=-998.0):
    """Write data to an HDF5 file using the python package H5PY. The file is appended to and in
    the case that a dataset already exists, it is overwritten.

    Parameters:
        data: A 1-D np.array/pd.Series or a ``pd.DataFrame`` containing different variables as
            columns
        flname (str): Path of the HDF5 you wish to write to or create
        h5path (str): Forward slash (/) delimited path through the group hierarchy you wish to
            place the dataset(s) specified by the argument ``datasets`` into. The dataset name
            cannot be passed using this argument, it is interpreted as a group name. A value of
            ``None`` places the dataset into the root directory of the HDF5 file.
        datasets (str or list): Name of the dataset(s) to write out. If a ``pd.DataFrame`` is
            passed, the values passed by the argument ``datasets`` must match the DataFrame's
            columns.
        dtype (str): The data type to write. Currently, only the following values are permitted:
            ``['int32', 'float32', 'float64']``. If a ``pd.DataFrame`` is passed and this argument
            is left to it's default value of ``None``, the DataFrame's dtypes must be of the types
            listed above.
        gridstr (str): Grid definition string that is saved to the HDF5 file as an attribute
            of the group defined by the parameter ``h5path``.
        trim_variable (str): Variable to use for trimming the data. An index will be written to the
            h5file and will be used to rebuild dataset while only nontrimmed data will be written
            out
        var_min (float): minimum trimming limit usedif trim_variable is passed

    Examples:
        Write a single pd.Series or np.array to an HDF5 file:

        >>> gs.write_h5(array, 'file.h5', h5path='Modeled/Var1', datasets='Realization_0001')

        Write a whole ``pd.DataFrame`` in group (folder) 'OriginalData' that contains a dataset
        for every column in the ``pd.DataFrame``:

        >>> gs.write_h5('file.h5', DataFrame, h5path='OriginalData')
    """
    import h5py
    # Sanity checks
    if not isinstance(data, (np.ndarray, pd.Series, pd.DataFrame)):
        raise TypeError("The data passed must be a np.ndarray, pd.Series, or pd.DataFrame")
    if isinstance(data, (np.ndarray, pd.Series)) and ((isinstance(datasets, list)) and
                                                      (len(datasets) > 1)):
        raise TypeError("The passed data is 1-D yet more than one variable is specified")

    # Sort out some defaults
    if h5path in [None, '', ' ', '/']:
        h5path = '/'

    if datasets is None:
        if isinstance(data, pd.DataFrame):
            datasets = list(data.columns)
        elif isinstance(data, pd.Series):
            datasets = [data.name]
        else:
            datasets = ['data']
    else:
        if not isinstance(datasets, list):
            datasets = [datasets]

    # More sanity checks
    if isinstance(data, pd.DataFrame) and datasets is not None:
        errcols = set(datasets) - set(data.columns)
        if len(errcols) > 0:
            raise KeyError("The following datasets cannot be found in the pd.DataFrame: %s" %
                           errcols)

    # Open up the h5 file
    with h5py.File(flname, 'w') as h5file:
        # Create the groups
        if h5path == '/':
            group = h5file
        else:
            group = h5file.require_group(h5path)

        trim = False
        if trim_variable:
            if isinstance(data, pd.DataFrame):
                array = data[trim_variable]
                trim = True
            elif isinstance(data, pd.Series):
                array = data
                trim = True
            else:
                printerr('Data needs to be either pd.DataFrame, or pd.Series to trim'
                         'Data will be written out untrimmed', errtype='warning')

        if trim:
            h5_index = array.index[array > var_min].values.astype('int32')
            # Make sure the last cell is always included
            last_index = len(array) - 1
            if h5_index.max() < last_index:
                h5_index = np.append(h5_index, [last_index])
            array = np.atleast_2d(h5_index).T
            group.require_dataset(name='h5_index', data=array, shape=array.shape,
                                  dtype='int32')
        # Iterate over the data if required and write it to the HDF5 data
        for dataset in datasets:
            if isinstance(data, pd.DataFrame):
                array = data[[dataset]]
                dtype = str(array.dtypes[0])
            else:
                array = data
            if trim:
                array = array.loc[h5_index]
            # Determine default data type if required
            if dtype is None:
                dtype = str(array.dtype)
            if dtype == 'int64' or dtype == 'int32' or dtype == 'float32' or dtype == 'float64':
                if dtype == 'int64':
                    printerr("The data type 'int64' is not supported by the Fortran module"
                             " hdf5_io. The data will be written as 'int32'.",
                             errtype='warning')
                    array = array.astype('int32')
                array = np.atleast_2d(array).T
                group.require_dataset(name=dataset, data=array, shape=array.shape, dtype=dtype)
            else:
                raise NotImplementedError("The data type `%s` is not supported by the Fortran"
                                          " module hdf5_io. Supported data types: ['int32',"
                                          "'float32', 'float64']" % dtype)

        # Write the griddef if required
        if isinstance(gridstr, bytes):
            h5file[h5path].attrs['griddef'] = '%s' % gridstr


def read_h5(flname, h5path=None, datasets=None, fill_value=-999):
    """
    Return a 1-D array from an HDF5 file or build a ``pd.DataFrame()`` from a list of datasets in a
    single group.

    The argument ``h5path`` must be a path to a group. If 1 or more specific variables are desired
    to be loaded, pass a list to ``datasets`` to specify which to read.

    Parameters:
        flname (str): Path of the HDF5 you wish to write to or create
        h5path (str): Forward slash (/) delimited path through the group hierarchy you wish to
            read the dataset(s) specified by the argument ``datasets`` from. The dataset name
            cannot be passed using this argument, it is interpreted as a group name only. A
            value of ``None`` places the dataset into the root directory of the HDF5 file. A value
            of ``False`` loads a blank pd.DataFrame().
        datasets (str or list): Name of the dataset(s) to read from the group specified by
            ``h5path``. Does nothing if ``h5path`` points to a dataset.
        fill_value (float or np.NaN): value to fill in grid with if trimmed data was written out.
            default is -999

    Returns:
        data (pd.DataFrame): DataFrame containing one or more columns, each containing a single
        1-D array of a variable.
    """
    import h5py
    # Handle some default parameters
    trim = False
    if h5path in [None, True, '', ' ', '/']:
        h5path = '/'
    elif h5path is False:
        return pd.DataFrame()
    else:
        h5path = '/' + h5path
    if isinstance(datasets, str):
        datasets = [datasets]

    # Sanity checks
    if not h5py.is_hdf5(flname):
        raise IOError("The passed file path is not a valid HDF5 file.")

    # Open the file and read the data
    with h5py.File(flname, 'r') as h5store:
        # Sort out some defaults and make sure all the datasets required are present
        if isinstance(h5store[h5path], h5py.Group):
            group = True
            h5dsets = [x for x in list(h5store[h5path].keys()) if
                       isinstance(h5store['%s/%s' % (h5path, x)], h5py.Dataset)]
            if 'h5_index' in h5dsets:
                trim = True
                h5dsets.remove('h5_index')
            if datasets is None:
                datasets = h5dsets
            else:
                error = set(datasets) - set(h5dsets)
                if len(error) > 0:
                    raise KeyError("The following datasets cannot be found in h5path specified: %s"
                                   % error)
        elif isinstance(h5store[h5path], h5py.Dataset):
            group = False
            if h5path.endswith('/'):
                h5path = h5path[:-1]
            datasets = [h5path.rsplit('/', 1)[1]]
            if 'h5_index' in datasets:
                trim = True
                datasets.remove('h5_index')
        else:
            raise ValueError("The `h5path` does not point to a dataset or a group.")

        # Read the data
        data = []
        if trim:
            if group:
                path = h5path + '/h5_index'
            else:
                path = '/h5_index'
            # python
            h5_index = h5store[path][:].flatten()

        for dataset in datasets:
            if group:
                path = h5path + '/' + dataset
            else:
                path = h5path
            # python
            data.append(h5store[path][:].flatten())
        if len(data) > 0:
            data = np.column_stack(data)
            if trim:
                data = pd.DataFrame(data, columns=datasets, index=h5_index)
                data = data.reindex(index=range(h5_index.max() + 1), fill_value=fill_value)
            else:
                data = pd.DataFrame(data, columns=datasets)
        else:
            printerr("No data was found within the root directory. A value of `None` has been"
                     " returned", errtype='warning')
            data = None
    return data


def ish5dataset(h5fl, dataset, h5path=None):
    """
    Check to see if a dataset exits within an HDF5 file

    The argument ``h5path`` must be a path to a group and cannot contain the dataset name.
    Can only check for one dataset at a time.

    Parameters:
        flname (str): Path of the HDF5 you wish to check
        h5path (str): Forward slash (/) delimited path through the group hierarchy you wish to
            check for the specified dataset. The dataset name cannot be passed using this
            argument, it is interpreted as a group name only. A value of ``None`` places the
            dataset into the root directory of the HDF5 file.
        dataset (str): Name of the dataset to check for in the group specified by ``h5path``.

    Returns:
        exists (bool): Indicator if the specified dataset exists
    """
    import h5py
    # Handle some default parameters
    if h5path in [None, True, '', ' ', '/']:
        h5path = '/'
    else:
        h5path = '/' + h5path + '/'
    h5path = h5path + dataset
    # Sanity checks
    if os.path.isfile(h5fl):
        if not h5py.is_hdf5(h5fl):
            raise IOError("The passed file exists but s not a valid HDF5 file.")
        # Check for the dataset
        with h5py.File(h5fl, 'r') as h5store:
            if dataset in [x for x in h5store.keys()]:
                exists = isinstance(h5store[h5path], h5py.Dataset)
            else:
                exists = False
    else:
        exists = False

    return exists


class H5Store:
    """
    A simple class within pygeostat to manage and use HDF5 files.

    :ivar str flname:  Path to a HDF5 file to create or use
    :ivar ``h5py.File`` h5data: h5py File object
    :ivar dict paths: Dictionary containing all of the groups found in the HDF5 file that
        contain datasets

    Parameters:
        flname (str): Path to a HDF5 file to create or use

    Usage:

        Write a np.array or pd.Series to the HDF5 file:

        >>> H5Store['Group1/Group2/Var1'] = np.array()

        Write all the columns in a ``pd.DataFrame`` to the HDF5 file:

        >>> H5Store['Group1/Group2'] = pd.DataFrame()

        Retrieve a single 1-D array:

        >>> array = H5Store['Group1/Group2/Var1']

        Retrieve a single 1-D array within the root directory of the HDF5 file:

        >>> array = H5Store['Var1']

        Retrieve the first value from the array:

        >>> value = H5Store['Var1', 0]

        Retrieve a slice of values from the array:

        >>> values = H5Store['Var1', 10:15]
    """
    import h5py

    def __init__(self, flname, replace=False):

        self.flname = flname
        if self.h5py.is_hdf5(flname):
            if replace:
                os.remove(flname)
                self.h5data = None
                self.paths = None
            else:
                self._loadfile()
        else:
            self.h5data = None
            self.paths = None

    def __str__(self):
        """
        Print a nice list of groups and the datasets found within them using the variable
        ``self.paths``.

        Example:
            Print any groups found within the HDF5 file and the datasets within:

            >>> print(H5Store)
        """
        printout = (' ' * 21) + 'Groups and datasets found in HDF5 file\n' + ('-' * 80) + '\n'
        for group in sorted(self.paths):
            printout = printout + ('Group Path: \'%s\'' % group) + '\n'
            printout = printout + textwrap.fill('Datasets:   %s' % self.paths[group], width=80,
                                                subsequent_indent='             ') + '\n'
            printout = printout + ('-' * 80) + '\n'

        return printout

    def __getitem__(self, key):
        """
        Retrieve an array using the self[key] notation. The passed key is the path used to access
        the array desired and included direction through groups if required and the dataset name.
        The array may be selectively queried allowing a specific value or range of values to be
        loaded into the systems memory and not the whole array.

        Example:

            Retrieve a single 1-D array:

            >>> array = H5Store['Group1/Group2/Var1']

            Retrieve a single 1-D array within the root directory of the HDF5 file:

            >>> array = H5Store['Var1']

            Retrieve the first value from the array:

            >>> value = H5Store['Var1', 0]

            Retrieve a slice of values from the array:

            >>> values = H5Store['Var1', 10:15]
        """

        # Manage the input
        rng = None
        if isinstance(key, tuple):
            key, rng = key
        if key in ['', ' ']:
            key = '/'
        if not key.endswith("/"):
            key += "/"

        # Sanity checks
        if self.h5data is None:
            raise FileNotFoundError("No HDF5 file exists yet")
        if rng is not None:
            if not isinstance(self.h5data[key], self.h5py.Dataset):
                raise KeyError("The passed key must lead to a H5 dataset if a slice is desired")

        # Load the full group or just a single array
        if isinstance(self.h5data[key], self.h5py.Group):
            datasets = [x for x in list(self.h5data[key].keys()) if
                        isinstance(self.h5data['%s%s' % (key, x)], self.h5py.Dataset)]
            data = []
            for dataset in datasets:
                data.append(self.h5data['%s/%s' % (key, dataset)][:].flatten())
            return pd.DataFrame(np.column_stack(data), columns=datasets)
        elif rng is not None:
            return self.h5data[key][rng]
        else:
            return self.h5data[key][:]

    def __setitem__(self, key, value):
        """
        Write the the HDF5 file using the self[key] notation.

        If a pd.Series or np.array is passed, the last entry in the path is used as the dataset
        name. If a ``pd.DataFrame`` is passed, all columns are written to the path specified to
        datasets with their names retrieved from the ``pd.DataFrame``'s columns. If more flexible
        usage is required, please use :func:`gs.write_h5()<pygeostat.data.h5_io.write_h5>`.

        Example:
            Write a np.array or pd.Series to the HDF5 file:

            >>> H5Store['Group1/Group2/Var1'] = np.array()

            Write all the columns in a ``pd.DataFrame`` to the HDF5 file:

            >>> H5Store['Group1/Group2'] = pd.DataFrame()
        """

        # Handle different types of input objects
        if isinstance(value, pd.DataFrame):
            h5path = key
            variables = None
        elif isinstance(value, pd.Series):
            h5path = key
            variables = [value.name]
        else:
            key = key.rsplit('/', 1)
            if len(key) > 1:
                h5path, variables = key
            else:
                h5path = None
                variables = key[0]

        # Write the data to the HDF5 file
        write_h5(value, self.flname, h5path, variables)

        # Reload the file
        if isinstance(self.h5data, self.h5py.File):
            self.h5data.close()
        self._loadfile()

    def __enter__(self):
        " required for `with gs.H5Store(file) as h5store:` pattern "
        return self

    def __exit__(self, type, value, traceback):
        " required for `with gs.H5Store(file) as h5store:` pattern "
        self.close()

    def _loadfile(self):
        """
        Load the HDF5 file as h5py object and get the dataset paths
        """
        self.h5data = self.h5py.File(self.flname, 'r+')
        self.paths = self._get_paths()

    def _get_paths(self):
        """
        Build lookup dictionary of all the paths to the datasets found in a HDF5 file.
        """
        done = False
        paths = {}
        groups = ['/']
        while not done:
            for group in groups:
                children = ['%s%s/' % (group, x)for x in list(self.h5data[group].keys()) if
                            isinstance(self.h5data['%s%s' % (group, x)], self.h5py.Group)]
                groups.extend(children)
                datasets = [x for x in list(self.h5data[group].keys()) if
                            isinstance(self.h5data['%s%s' % (group, x)], self.h5py.Dataset)]
                if len(datasets) > 0:
                    paths[group] = datasets
                groups.remove(group)
            if len(groups) == 0:
                done = True

        return paths

    def close(self):
        """
        Release the open HDF5 file from python.
        """
        try:
            self.h5data.close()
        except:
            print("The H5 file wasn't open")

    def datasets(self, h5path=None):
        """
        Return the datasets found in the specified group.

        Keyword Arguments:
            h5path (str): Forward slash (/) delimited path through the group hierarchy you wish to
                retrieve the lists of datasets from. A dataset name cannot be passed using this
                argument, it is interpreted as a group name. A value of ``None`` places the dataset
                into the root directory of the HDF5 file.

        Returns:
            datasets (list): List of the datasets found within the specified `h5path`
        """
        h5path = _fixh5path(h5path)
        return self.paths[h5path]

    def iteritems(self, h5path=None, datasets=None, wildcard=None):
        """
        Produces an iterator that can be used to iterate over HDF5 datasets.

        Can use the parameter ``h5path`` to indicate which group to retrieve the datasets from.
        If a set of specific datasets are required, the parameter ``datasets`` will restrict the
        iterator to those. The parameter ``wildcard`` allows a string wild-card value to restrict
        which datasets are iterated over.

        Keyword Arguments:
            h5path (str): Forward slash (/) delimited path through the group hierarchy you wish to
                retrieve datasets from. A dataset name cannot be passed using this argument, it is
                interpreted as a group name. A value of ``None`` places the dataset into the root
                directory of the HDF5 file.
            datasets (list): List of specific dataset names found within the specified group to
                iterator over
            wildcard (str): String to search for within the names of the datasets found within the
                specified group to iterate over

        Examples:

            Load a HDF5 file to pygeostat:

            >>> data = gs.H5Store('data.h5')

            Iterate over all datasets within the root directory of a HDF5 file:

            >>> for dataset in data.iteritems():
            >>>     gs.histplt(dataset)

            Iterate over the datasets within a specific group that are realizations:

            >>> for dataset in data.iteritems(h5path='Simulation/NS_AU', wildcard='Realization'):
            >>>     gs.histplt(dataset)
        """
        # Handle some parameters
        h5path = _fixh5path(h5path)
        if isinstance(datasets, str):
            datasets = [datasets]
        # Sanity checks
        if datasets and wildcard:
            raise ValueError("Only `datasets` or `wildcard` can be passed")
        if datasets:
            error = set(datasets) - set(self.paths[h5path])
            if len(error) > 0:
                raise ValueError("The following datasets were specified but not found in the H5"
                                 " group: %s" % error)
        elif wildcard:
            datasets = []
            for dataset in self.paths[h5path]:
                if wildcard in dataset:
                    datasets.append(dataset)
            if len(datasets) == 0:
                raise ValueError("No datasets with the specified `wildcard` were found")
        else:
            datasets = self.paths[h5path]
            if len(datasets) == 0:
                raise ValueError("No datasets were found in the specified `h5path`")
        # Produce the generator
        for dataset in datasets:
            yield self.__getitem__(dataset)
