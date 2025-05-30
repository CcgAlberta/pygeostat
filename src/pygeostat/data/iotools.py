#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
iotools.py: Contains input/output utilities/functions for pygeostat. Many of which
are based off of Pandas builtin functions.
'''
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import warnings

import pandas as pd
import numpy as np
from .. pygeostat_parameters import Parameters

def read_file(flname, fltype=None, headeronly=False, delimiter=r'\s*', h5path=None, h5datasets=None,
             columns=None, ireal=1, griddef=None, tmin=None):
    '''
    Reads in a GSLIB-style Geo-EAS data file, CSV, GSB or HDF5 data files.

    Parameters:
        flname (str): Path (or name) of file to read.

    Keyword Args:
        fltype (str): Type of file to read: either ``csv``, ``gslib``, or ``hdf5``.
        headeronly (bool): If True, only reads in the 1st line from the data file
                     which is useful for just getting column numbers or testing. OR
                     it allows you to open a hdf5 object with Pandas HDFStore functionality
        delimiter (str): Delimiter specified instead of sniffing
        h5path (str): Forward slash (/) delimited path through the group hierarchy you wish to
            read the dataset(s) specified by the argument ``datasets`` from. The dataset name
            cannot be passed using this argument, it is interpreted as a group name only. A
            value of ``None`` places the dataset into the root directory of the HDF5 file. A value
            of ``False`` loads a blank pd.DataFrame().
        h5datasets (str or list): Name of the dataset(s) to read from the group specified by
            ``h5path``. Does nothing if ``h5path`` points to a dataset.
        column (list): List of column labels to use for resulting frame
        ireal (int): Number of realizaitons in the file
        griddef (GridDef): griddef for the realization
        tmin (float): values less than this number are convernted to NaN, since NaN's are
            natural handled within matplotlib, pandas, numpy, etc. If None, set to
            pygeostat.Parameters['data.tmin'].

    Returns:
        data (pandas.DataFrame): Pandas DataFrame object with input data.

    Note:
        Functions can also be called seperately with the following code

        >>> data.data = pygeostat.read_gslib(flname)
        >>> data.data = pygeostat.read_csv(flname)
        >>> data.data = pygeostat.read_h5(flname, h5path='')
        >>> data.data = pygeostat.read_gsb(flname)
        >>> data.data = pygeostat.open_hdf5(flname)

    Examples:
		>>> data.data = gs.read_gsb('testgsb.gsb')
		>>> data = gs.DataFile('testgsb.gsb')

    '''
    # Infer filetype if none specified based on the file extention
    if fltype is None:
        import os
        _, extention = os.path.splitext(flname)
        extention = extention.lower()
        if extention in ['.out', '.data']:
            fltype = 'gslib'
        elif extention in ['.csv']:
            fltype = 'csv'
        elif extention in ['.h5', '.hdf5', '.hd5']:
            fltype = 'hdf5'
        elif flname.lower().endswith('gsb') or flname.endswith('GSB'):
            fltype = 'gsb'
        else:
            # Otherwise just assume GSLIB for now but this may change
            fltype = 'gslib'

    # Call specific read functions based on data type
    if fltype == 'gslib':
        data = read_gslib(flname=flname, headeronly=headeronly, delimiter=delimiter,
                          tmin=tmin)
    elif fltype == 'csv':
        data = read_csv(flname=flname, headeronly=headeronly, tmin=tmin)
    elif fltype == 'hdf5' or fltype == 'hd5' or fltype == 'h5':
        from .h5_io import read_h5
        if headeronly:
            data = pd.DataFrame()
        else:
            # Note that tmin isn't implemented for h5 yet, due to my inexperience
            # with this underlying module/format
            data = read_h5(flname=flname, h5path=h5path, datasets=h5datasets)
    elif fltype == 'gsb':
        data = read_gsb(flname=flname, ireal=ireal, tmin=tmin)
    else:
        print('File type unsupported! Try "gslib", "csv" or "hdf5"')

    if columns is not None:
        data.columns = columns

    return data


def _test_file_open(filename):
    """ test and raise if a file cannot be opened """
    try:
        with open(filename, 'r'):
            pass
    except FileNotFoundError:
        raise FileNotFoundError('{} does not exist!'.format(filename))
    except IOError:
        raise IOError('Could not open {}'.format(filename))


def read_gslib(flname, headeronly=False, delimiter=r'\s*', tmin=None):
    '''Reads in a GSLIB-style Geo-EAS data file

    Parameters:
        flname (str): Path (or name) of file to read.

    Keyword Args:
        headeronly (bool): If True, only reads in the 1st line from the data file
                     which is useful for just getting column numbers or testing
        delimiter (str): Delimiter specified instead of sniffing
        tmin (float): values less than this number are convernted to NaN, since NaN's are
            natural handled within matplotlib, pandas, numpy, etc. If None, set to
            pygeostat.Parameters['data.tmin'].

    Returns:
        data (pandas.DataFrame): Pandas DataFrame object with input data.

    '''
    # Can the file be opened?
    _test_file_open(flname)
    # Only read in the header + 1 line of the data file?
    if headeronly:
        nrows = 1
    else:
        nrows = None
    # Reading the header for a GSLIB "Geo-EAS" format file
    with open(flname, 'r') as datafl:
        # Geo-EAS Header
        _ = datafl.readline().strip()  # Title
        nvar = int(datafl.readline().split()[0])
        varnames = []
        for _ in range(nvar):
            varnames.append(datafl.readline().strip())
        nrowtoskip = 2 + nvar

    engine = Parameters['data.io.pandas_engine']
    # Use the new pandas chunked reading method
    # with the "python" engine since the C engine sometimes does not work with
    # realizations with trailing spaces
    try:
        tpdf = pd.read_csv(flname, skiprows=nrowtoskip, header=None,
                           delimiter=delimiter, skipinitialspace=True,
                           nrows=nrows, engine=engine, chunksize=100000)
        data = pd.concat(tpdf, ignore_index=True)
    # Fallback to old method using the python engine
    except (ValueError, NotImplementedError):
        data = pd.read_csv(flname, skiprows=nrowtoskip, header=None, delimiter=r'\s*',
                           skipinitialspace=True, nrows=nrows, engine='python')

    # Replace values below tmin with nan
    data = _data_trim(data, tmin=tmin)

    # Assign variable names
    data.columns = varnames

    # return only the Pandas DataFrame
    return data


def read_csv(flname, headeronly=False, tmin=None):
    '''Reads in a GSLIB-style CSV data file.

    Parameters:
        flname (str): Path (or name) of file to read.

    Keyword Args:
        headeronly (bool): If True, only reads in the 1st line from the data file
                     which is useful for just getting column numbers or testing
        delimiter (str): Delimiter specified instead of sniffing
        tmin (float): values less than this number are convernted to NaN, since NaN's are
            natural handled within matplotlib, pandas, numpy, etc. If None, set to
            pygeostat.Parameters['data.tmin'].

    Returns:
        data (pandas.DataFrame): Pandas DataFrame object with input data.

    '''
    if headeronly:
        nrows = 1
    else:
        nrows = None
    # Reading a GSLIB "Geo-EAS" format file
    tpdf = pd.read_csv(flname, header=0, nrows=nrows, engine='c', chunksize=100000)
    data = pd.concat(tpdf, ignore_index=True)
    # Replace values below tmin with nan
    data = _data_trim(data, tmin=tmin)
    # Only return the DataFrame
    return data


def compile_pygsb():
    '''
    Compiles 'pygeostat/fortran/src/pygsb.f90' using 'pygeostat/fortran/compile.py'
    and tries to import pygsb.pyd

    Note:
        How to install a gfortran compiler:

        - Install chocolatey from:

            chocolatey.org/install

        (chocolatey is a package manager that let you install software using command prompt and PowerShell)
        
        - After installing chocolatey, then install the ‘gnu Fortran compiler’ by writing the below in a PowerShell:

            choco install mingw --version 8.1

            choco install visualstudio2019community

            choco install visualstudio2019-workload-vctools
        - When  installing "mingw" through "chocolatey", ensure that the path of the "mingw" 's "bin" folder is added to the environment variables path.

    '''
    import os
    import subprocess
    cwd1 = os.path.abspath(os.path.join(os.path.dirname(__file__), '../fortran'))
    
    if not os.path.isfile(os.path.join(cwd1,'pygsb.pyd')):
        compiler = 'gnu'
        subprocess.call('python compile.py -clean pygsb', cwd=cwd1)
        subprocess.call('python compile.py -compiler={} pygsb'.format(compiler), cwd=cwd1)
    else:
    # Ensure the right version of .pyd exists.     
        try:
            from ..fortran import pygsb as pygsb
        except ImportError:
            compiler = 'gnu'
            subprocess.call('python compile.py -clean pygsb', cwd=cwd1)
            subprocess.call('python compile.py -compiler={} pygsb'.format(compiler), cwd=cwd1)
    try:
        from ..fortran import pygsb as pygsb
    except ImportError:
        raise ImportError("Could not import 'pygsb' from 'pygeostat.fortran'('pygsb.f90' did not compile to create 'pygsb.pyd'. Consider installing a gfortran compiler with tools such as mingw (see documnetation to help you how install the Fortran compiler).")
    return pygsb


def isbinary(file):
    """
    From http://stackoverflow.com/a/7392391/5545005
    Its hard to understand what's going on here.. but it seems to work for gsb files ....
    H5 has a handy check but this fills the gap for gsb files when trying to read ascii
    """
    textchars = bytearray({7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)) - {0x7f})

    def is_binary_string(bytes):
        return bool(bytes.translate(None, textchars))
    return is_binary_string(open(file, 'rb').read(1024))


def read_gsb(flname, ireal=-1, tmin=None, null=None):
    '''Reads in a CCG GSB (GSLIB-Binary) file.

    Parameters:
        flname (str): Path (or name) of file to read.

    Keyword Args:
        ireal (int): 1-indexed realization number to read (reads 1 at a time), -1 to read all
        tmin (float): values less than this number are convernted to NaN, since NaN's are
            natural handled within matplotlib, pandas, numpy, etc. If None, set to
            pygeostat.Parameters['data.tmin'].
        null (float): when the gsb array has a keyout, on reconstruction this value fills the array
            in keyed out locations. If `None` taken from Parameters['data.null']

    Returns:
        data (pandas.DataFrame): Pandas DataFrame object with input data.

    .. codeauthor:: Jared Deutsch 2016-02-19
    '''
    pygsb = compile_pygsb()
    # Can the file be opened?
    _test_file_open(flname)

    if not isbinary(flname):
        raise ValueError('The file %s appears to be non-binary formatted ' % flname)

    # Load required header information
    nvar, nx, ny, nz, nreal, errorvalue = pygsb.pyreadgsbheader(flname)
    if errorvalue != 0:
        raise AssertionError("Error reading GSB header!, Error #{}".format(errorvalue))

    # RMB - fixing an observed bug, which occurs if any of the dimensions are zero
    if nx == 0:
        nx = 1
    if ny == 0:
        ny = 1
    if nz == 0:
        nz = 1
    if ireal == 0:
        raise ValueError(('ireal shoud be greater than 1 (read a specified 1-index realization)'
                          'or -1 (read all realizations)'))
    if null is None:
        null = Parameters.get('data.null')
    if null is None:
        null = -999.99  # matching the GSB fortran defaults
    # Get the data
    nxyz = nx * ny * nz
    if ireal <= -1:
        errorvalue, reals, vnames = pygsb.pyreadgsbdata(flname, nvar, nxyz, 1, null)
        for ireal in range(1, nreal):
            errorvalue, values, vnames = pygsb.pyreadgsbdata(
                flname, nvar, nxyz, ireal + 1, null)
            reals = np.append(reals, values, axis=1)
    else:
        errorvalue, reals, vnames = pygsb.pyreadgsbdata(flname, nvar, nxyz, ireal, null)
    if errorvalue != 0:
        raise AssertionError("Error reading GSB data!, Error #{}".format(errorvalue))

    # Clean up column names and trimmed data
    vnames = [''.join([v.decode("utf-8") for v in vname]).strip() for vname in vnames]

    # Convert to pandas dataframe
    data = pd.DataFrame(data=reals.transpose(), columns=vnames)

    # Replace values below tmin with nan
    data = _data_trim(data, tmin)

    # Only return the DataFrame
    return data


def _data_trim(data, tmin):
    """Replace values less than tmin with NaN. This routine is private
       to the iotools module, since it is only intended for use on the
       import of data.

    Parameters:
        data (pandas.DataFrame): data to trim
        tmin (int or float): values less than tmin are assigned np.nan

    """
    from pandas.api.types import is_numeric_dtype
    if tmin is None:
        # Use the Parameters tmin, which may also be None (leading to no trimming)
        tmin = Parameters['data.tmin']
    if tmin is None or not tmin:
        return data
    if not isinstance(tmin, int) and not isinstance(tmin, float):
        raise ValueError('tmin must be an integer or float!')
    if tmin is not None:
        if Parameters['data.fix_legacy_null'] and Parameters['data.null'] is not None:
            warned = False
            for legacy_null in Parameters['data.legacy_null']:
                if not warned and legacy_null != Parameters['data.null']:
                    for col in data.columns:
                        if is_numeric_dtype(data[col]):
                            if (np.isclose(data[col], legacy_null)).any():
                                warnings.warn(
                                    f'found {legacy_null} in {col} during read!')
                                warned = True
                                break
                data.replace(legacy_null, Parameters['data.null'], inplace=True)
        for col in data.columns:
            if is_numeric_dtype(data[col]):
                data.loc[data[col] < tmin, col] = np.nan
    return data


def write_gslib(data, flname, title=None, variables=None,
                fmt=None, sep=' ', null=None):
    """Writes out a GSLIB-style data file.

    Parameters:
        data (pygeostat.DataFile or pandas.DataFrame): data to write out
        flname (str): Path (or name) of file to write out.

    Keyword Args:
        title (str): Title for output file.
        variables (List(str)): List of variables to write out if only a subset is desired.
        fmt (str): Format to use for floating point numbers.
        sep (str): Delimiter to use for file output, generally don't need to change.
        null (float): NaN numbers are converted to this value prior to writing. If None, set
            to data.null. If data.Null is None, set to pygeostat.Parameters['data.null'].

    """
    from .data import DataFile as DataFile
    data = _data_fillnan(data, null)
    # If a DataFile is used, check the arguments
    if isinstance(data, DataFile):
        flname, variables, sep, fltype = data.check_datafile(flname, variables,
                                                             sep, 'gslib')
    else:
        # If variables is none, then get the columns
        if variables is None:
            variables = data.columns.tolist()
    # If title is STILL none after introspection, then at least include something
    if title is None:
        title = 'pygeostat_saved_data'
    if fmt is None:
        fmt = Parameters.get('data.write.python_floatfmt', '%.5f')
    # Force a write
    with open(flname, 'w') as outfl:
        # GSLIB Header
        outfl.write(title.strip() + '\n' + str(len(variables)) + '\n')
        for variable in variables:
            outfl.write(str(variable) + '\n')
        # GSLIB separated data
        if isinstance(data, DataFile):
            data.data[variables].to_csv(outfl, header=False, index=False, sep=sep,
                                        float_format=fmt, lineterminator='\n')
        else:
            data[variables].to_csv(outfl, header=False, index=False, sep=sep,
                                   float_format=fmt, lineterminator='\n')


def write_csv(data, flname, variables=None,
              fmt='%.5f', sep=',', fltype='csv', null=None):
    """Writes out a CSV or Excel (XLSX) data file.

    Parameters:
        data (pygeostat.DataFile or pandas.DataFrame): data to write out
        flname (str): Path (or name) of file to write out.

    Keyword Args:
        variables (List(str)): List of variables to write out if only a subset is desired.
        fmt (str): Format to use for floating point numbers.
        sep (str): Delimiter to use for file output, generally don't need to change.
        fltype (str): Type of file to write either ``csv`` or ``xlsx``.
        null (float): NaN numbers are converted to this value prior to writing. If None, set
            to data.null. If data.Null is None, set to pygeostat.Parameters['data.null'].

    """
    from .data import DataFile as DataFile
    data = _data_fillnan(data, null)
    # If a DataFile is used, check the arguments
    if isinstance(data, DataFile):
        flname, variables, sep, fltype = data.check_datafile(flname, variables,
                                                             sep, fltype)
    else:
        # If variables is none, then get the columns
        if variables is None:
            variables = data.columns.tolist()
    # Force a write
    with open(flname, 'w') as outfl:
        if fltype.lower() == 'csv':
            # CSV using pandas native CSV writer
            if isinstance(data, DataFile):
                data.data[variables].to_csv(outfl, header=True, index=False, sep=sep,
                                            float_format=fmt, lineterminator='\n')
            else:
                data[variables].to_csv(outfl, header=True, index=False, sep=sep,
                                       float_format=fmt, lineterminator='\n')
        elif (fltype.lower() == 'xlsx') or (fltype.lower() == 'excel'):
            # Excel file writer - only xlsx is supported naturally in pygeostat,
            # but XLWT could be used if XLS files needed to be created
            # NOTE: this requires the Python package openpyxl
            if isinstance(data, DataFile):
                data.data[variables].to_excel(flname, header=True, index=False,
                                              float_format=fmt)
            else:
                data[variables].to_excel(flname, header=True, index=False,
                                         float_format=fmt)


def write_gsb(data, flname, tvar=None, nreals=1, variables=None, griddef=None, fmt=0):
    """
    Writes out a GSB (GSLIB-Binary) style data file. NaN values of tvar are compressed
    in the output with no tmin now provided.

    Parameters:
        data (pygeostat.DataFile or pandas.DataFrame): data to write out
        flname (str): Path (or name) of file to write out.
        tvar (str): Variable to trim by or None for no trimming. Note that all variables are
            trimmed in the data file (for compression) when this variable is trimmed.
        nreals (int): number of realizations in data

    Keyword Args:
        griddef (pygeostat.griddef.GridDef): This is required if the data is gridded and you
            want other gsb programs to read it
        fmt (int): if 0 then will write out all variables as float 64. Otherwise should be
            an list with a length equal to number of variables and with the following format codes
            1=int32, 2=float32, 3=float64
        variables (List(str)): List of variables to write out if only a subset is desired.

    .. codeauthor:: Jared Deutsch 2016-02-19, modified by Ryan Barnett 2018-04-12
    """
    from .data import DataFile as DataFile
    pygsb = compile_pygsb()
    null = Parameters.get('data.null', None)
    data = _data_fillnan(data, null)
    # If variables is none, then get the columns
    # Also configure the data for output
    if not isinstance(data, DataFile):
        if variables is None:
            variables = data.columns.tolist()
        datamat = np.array(data[variables].values).transpose()
    else:
        if griddef is None:
            if data.griddef is not None:
                griddef = data.griddef
        if variables is None:
            variables = data.data.columns.tolist()
        else:
            if isinstance(variables, str):
                variables = [variables]
        datamat = np.array(data.data[variables].values).transpose()
    # format
    if fmt == 0 or fmt is None:
        vkinds = [3 for x in range(len(variables))]
    elif not isinstance(fmt, list):
        raise ValueError("fmt needs to be a list. You passed a %s" % type(fmt))
    else:
        vkinds = []
        for f in fmt:
            vkinds.append(f)
    # Dimensioning
    nvar = len(variables)
    if griddef is not None:
        nx = griddef.nx
        ny = griddef.ny
        nz = griddef.nz
        if np.size(datamat, 1) != griddef.count() * nreals:
            raise ValueError("the passed data has the wrong  number of elements for nx, ny, nz, "
                             "nreal: %i %i %i %i" % (nx, ny, nz, nreals))
    else:
        if nreals == 1:
            nx = len(data[variables[0]])
            ny = 1
            nz = 1
        else:
            nx = np.size(datamat, 1) / nreals
            ny = 1
            nz = 1

    # Trimming
    if isinstance(tvar, str):
        tvar = variables.index(tvar) + 1
        if tvar <= 0:
            raise ValueError("Could not trimming variable {} in variable list {}".format(tvar,
                                                                                         variables))
    elif tvar is None:
        tvar = 1
    elif isinstance(tvar, bool):
        if tvar:
            tvar = 1
        else:
            tvar = 0
    elif tvar == 0:
        tvar = 1
    else:
        raise ValueError("Invalid trimming variable {}".format(tvar))

    # Character conversion of strings
    cvariables = []
    for varname in [var.ljust(64) for var in variables]:
        cvariables.append([v for v in varname])

    tmin = Parameters.get('data.tmin', None)
    if tmin is None:
        tmin = -1e21
    # Can the file be opened?
    try:
        _test_file_open(flname)
    except FileNotFoundError:
        pass  # this case is okay .. I guess we're just checking if we can open it?
    # Save out the data
    errorvalue = pygsb.pywritegsbdata(gsbfl=flname, datamat=datamat, vnames=cvariables, nvar=nvar,
                                      nreal=nreals, nx=nx, ny=ny, nz=nz, tmin=tmin,
                                      tmax=1.0e21, tvar=tvar, vkinds=vkinds)
    if errorvalue != 0:
        raise AssertionError("Error writing GSB data!, Error #{}".format(errorvalue))


def write_vtk(data, flname, dftype=None, x=None, y=None, z=None, variables=None, griddef=None,
              null=None, vdtype=None, cdtype=None):
    """
    Writes out an XML VTK data file. A required dependency is pyevtk, which may be installed using
    the following command:

        >>> pip install pyevtk

    Users are also recommended to install the latest Paraview, as versions from 2017 were
    observed to have odd precision bugs with the XML format.

    Parameters:
        data (pygeostat.DataFile): data to write out
        flname (str): Path (or name) of file to write out (without extension)

    Keyword Args:
        dftype (str): type of datafile options ``grid`` or ``point``, which if None,
            is drawn from data.dftype
        x (str): name of the x-coordinate, which is used if ``point``. Drawn
            from data.x if the kwarg=None. If not provided by these means for ```sgrid```,
            calculated via ``sim.griddef.get_coordinates()``.
        y (str): name of the y-coordinate, which is used if ``point``. Drawn
            from data.y if the kwarg=None. If not provided by these means for ```sgrid```,
            calculated via sim.griddef.get_coordinates().
        z (str): name of the z-coordinate, which is used if ``point``. Drawn
            from data.z if the kwarg=None. If not provided by these means for ```sgrid```,
            calculated via sim.griddef.get_coordinates().
        griddef (pygeostat.GridDef): grid definition, which is required if ``grid``.
            Drawn from data.griddef if the kwarg=None.
        variables (list or str): List or string of variables to write out. If None, then all
            columns aside from coordinates are written out by default.
        null (float): NaNs are converted to this value prior to writing. If None, set to
            pygeostat.Parameters['data.null_vtk'].
        vdtype (dict(str)) : Dictionary of the format {'varname': dtype}, where dtype is
            a numpy data format. May be used for reducing file size, by
            specifying ``int``, ``float32``, etc. If a format string is provided instead
            of a dictionary, that format is applied to all variables. This is not applied
            to coordinate variables (if applicable). If None, the value is drawn from
            Parameters['data.write_vtk.vdtype'].
        cdtype (str) : Numpy format to use for the output of coordinates, where valid formats
            are ``float64`` (default) and ``float32``. The later is recommended for reducing
            file sizes, but may not provide the requisite precision for UTM coordinates. If None,
            the value is drawn from Parameters['data.write_vtk.cdtype'].

    ``dftype`` should be one of:

    1. 'point' (irregular points) where ``data.x``, ``data.y`` and ``data.z`` are columns in ``data.data``
    2. 'grid' (regular or rectilinear grid) where ``data.griddef`` must be initialized
    3. 'sgrid' (structured grid) where ``data.x``, ``data.y`` and ``data.z`` are columns in ``data.data``. ``data.griddef`` should also be initialized, although only ``griddef.nx``, ``griddef.ny`` and ``griddef.nz`` are utilized (since the grid is assumed to not be regular)

    """
    from copy import deepcopy
    # Ensure dependencies are satisfied
    try:
        from pyevtk.hl import pointsToVTK, gridToVTK
    except:
        raise Exception('Looks like pyevtk is not installed, use:\n' +
                        '>>> pip install pyevtk\n')
    # Place in a new variable so potential changes aren't returned
    dat1 = deepcopy(data)
    # Ensure valid VTK option
    if dftype is None:
        if dat1.dftype == 'point':
            dftype = 'point'
        elif dat1.dftype == 'grid':
            dftype = 'grid'
        elif dat1.dftype == 'sgrid':
            dftype = 'sgrid'
        else:
            raise ValueError(" specified data.dftype is not yet available or invalid;\n" +
                             " see the docstring!")
    else:
        if not any([dftype == test for test in ['point', 'grid', 'sgrid']]):
            raise ValueError(" specified dftype is not yet available or invalid;\n" +
                             " see the docstring!")
    # Ensure a grid is present if necessary
    if dftype in ['grid', 'sgrid']:
        if griddef is None:
            griddef = dat1.griddef
            if griddef is None:
                raise ValueError(('griddef must be a kwarg or exist as data.griddef if '
                                  'dftype is grid or sgrid!'))
    # Ensure coordinates are present if necessary
    if dftype in ['point', 'sgrid']:
        if x is None:
            x = dat1.x
        if y is None:
            y = dat1.y
        if z is None:
            z = dat1.z
        test = dict(x=x, y=y, z=z)
        if dftype == 'point':
            for coord in test.keys():
                if test[coord] is None:
                    raise ValueError(('{} must be a kwarg or exist as data.{} if'
                                      'dftype is point'.format(coord, coord)))
        else:
            t = [t is None for t in test.values()]
            if all(t):
                raise ValueError(('at least one grid coordinate should be irregular (specified '
                                  'as a column in data) if using dftype="sgrid". Specify one '
                                  'or multiple irregular coordinate columns, or use dftype="grid"'))
            else:
                # A coordinate is not provided, so assume it follows the regular grid
                tx, ty, tz = griddef.get_coordinates()
                if x is None:
                    x, dat1['X'] = 'X', tx
                if y is None:
                    y, dat1['Y'] = 'Y', ty
                if z is None:
                    z, dat1['Z'] = 'Z', tz
    # Replace nan with null values
    if null is None:
        null = Parameters['data.null_vtk']
    dat1 = _data_fillnan(dat1, null)
    # Check the grid length
    if dftype != 'point':
        nx = griddef.nx
        ny = griddef.ny
        nz = griddef.nz
        nxyz = nx * ny * nz
        if dat1.shape[0] != nxyz:
            raise ValueError('nx*ny*nz should be equal to dat1.shape[0]!')
    if vdtype is None:
        vdtype = Parameters['data.write_vtk.vdtype']
    if cdtype is None:
        cdtype = Parameters['data.write_vtk.cdtype']
    # Generate a dictionary of the variables, including their speficied precision
    columns = dat1.data.columns.tolist()
    if variables is None:
        variables = []
        for col in columns:
            if col == x or col == y or col == z:
                continue
            variables.append(col)
        if len(variables) == 0:
            # Paraview requires some dat1...
            variables = [x]
    elif isinstance(variables, str):
        variables = [variables]
    # Generate a dictionary of variables precisions
    vdtyped = {col: 'float64' for col in columns}
    if isinstance(vdtype, dict):
        for var in vdtype:
            vdtyped[var] = vdtype[var]
    elif isinstance(vdtype, str):
        for var in vdtyped:
            if var == x or var == y or var == var == z:
                vdtyped[var] = cdtype
            else:
                vdtyped[var] = vdtype
    vardict = dict()
    for var in variables:
        if any(var in i for i in columns):
            t = np.ascontiguousarray(dat1.data[var].values)
            vardict[var] = t.astype(vdtyped[var])
            if dftype != 'point':
                vardict[var] = vardict[var].reshape((nx, ny, nz), order='F')
        else:
            raise KeyError(var + " is not a column name in dat1.data!")
    # General checks and variable assembly completed - write out the
    # data type
    if dftype == 'point':
        # Require specified coodinates for each point
        if any(x in i for i in columns):
            x = np.ascontiguousarray(dat1.data[x].values)
            x = x.astype(cdtype)
        else:
            raise KeyError(" x is not a column name in data.data!")
        if any(y in i for i in columns):
            y = np.ascontiguousarray(dat1.data[y].values)
            y = y.astype(cdtype)
        else:
            raise KeyError(" y is not a column name in data.data!")
        if any(z in i for i in columns):
            z = np.ascontiguousarray(dat1.data[z].values)
            z = z.astype(cdtype)
        else:
            raise KeyError(" z is not a column name in data.data!")
    if dftype == 'point':
        pointsToVTK(flname, x, y, z, data=vardict)
    elif dftype == 'grid':
        # Rectilinear VTK type
        grid_arr = data.griddef.grid_array
        xsiz, ysiz, zsiz = grid_arr[2], grid_arr[5], grid_arr[8]
        xmin, ymin, zmin = grid_arr[1] - xsiz * .5, grid_arr[4] - \
            ysiz * .5, grid_arr[7] - zsiz * .5
        xmax, ymax, zmax = xmin + nx * xsiz, ymin + ny * ysiz, zmin + nz * zsiz
        x = np.arange(xmin, xmax + xsiz, xsiz, dtype=cdtype)
        y = np.arange(ymin, ymax + ysiz, ysiz, dtype=cdtype)
        z = np.arange(zmin, zmax + zsiz, zsiz, dtype=cdtype)
        gridToVTK(flname, x, y, z, cellData=vardict)
    else:
        pass


def _data_fillnan(data, null):
    """Replace nan values prior to writing out. Private since only intended
       for use within the iotools

    Parameters:
        data (pandas.DataFrame or pygeostat.DataFrame): data that may have NaN values
        null (int or float): NaN occurence are replaced with this value

    """
    import pygeostat as gs
    import copy
    if null is None:
        if isinstance(data, gs.DataFile):
            if data.null is None:
                # Use the Parameters null, which may also be None (leading to no nan assignment)
                null = Parameters['data.null']
            else:
                # Use the DataFile.null, which may also be None (leading to no nan assignment)
                null = data.null
        else:
            null = Parameters['data.null']
    if null is None or not null:
        return data
    if not isinstance(null, int) and not isinstance(null, float):
        raise ValueError('null must be an integer or float!')
    if null is not None:
        data = copy.deepcopy(data)
        if isinstance(data, gs.DataFile):
            data.data.fillna(value=null, inplace=True)
        elif isinstance(data, pd.DataFrame):
            data.fillna(value=null, inplace=True)
        else:
            # Shouldn't get here... but just in case
            raise ValueError('only DataFrame and DataFile is implemented!')
    return data


def write_hvtk(data, flname, griddef, variables=None):
    """
    Writes out an H5 file and corresponding xdmf file that Paraview can read. Currently only
    supports 3D gridded datasets. This function will fail if the length of the DataFile or
    DataFrame does not equal ``griddef.count()``.

    The extension xdmf is silently enforced. Any other extension passed is replaced.

    Parameters:
        data (pd.DataFrame): The DataFrame to writeout
        flname (str): Path (or name) of file to write out.
        griddef (GridDef): Grid definitions for the realizations to be written out
        variables (str or list): optional set of variables to write out from the DataFrame

    """
    import os
    import subprocess
    from .h5_io import write_h5
    # setup the temporary strings to write things too:
    temp_attr = """           <Attribute Name="{varname}"
             AttributeType="Scalar"
             Center="Cell">
           <DataItem Dimensions="{nx} {ny} {nz}" NumberType="Float" Precision="8" Format="HDF">
            {h5flname}:{attrvarpath}
           </DataItem>
         </Attribute>"""
    xmlstr = """<?xml version="1.0" ?>
    <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
    <Xdmf Version="2.0">
     <Domain>
       <Grid Name="mesh" GridType="Uniform">
         <Topology TopologyType="3DCoRectMesh" Dimensions="{nz} {ny} {nx}"></Topology>
          <Geometry Type="ORIGIN_DXDYDZ">
            <DataItem Format="XML" Dimensions="{dim}">{zmn} {ymn} {xmn}</DataItem>
            <DataItem Format="XML" Dimensions="{dim}">{zsiz}  {ysiz}  {xsiz}</DataItem>
          </Geometry>
{attributes}
      </Grid>
     </Domain>
    </Xdmf>
    """
    # Checks
    if not isinstance(data, pd.DataFrame):
        print('ERROR: Pass a pd.DataFrame for writeout!')
        return
    if 1 in (griddef.nx, griddef.ny, griddef.nz) or len(data) != griddef.count():
        print('ERROR: This function only supports 3D grids! Check data or griddefs!')
        return

    attributes = ''
    ext = flname[flname.rfind('.'):]
    if ext.lower() != '.xdmf':
        flname = flname.replace(ext, '.xdmf')
        print('WARNING: Paraview h5-reader likes xdmf extensions - \nnew file'
              ' is: %s' % flname)
    outh5name = flname.replace('.xdmf', '.hvtk')

    # Always overwrites the old hvtk file, but has to clear it first
    if os.path.isfile(outh5name):
        try:
            os.remove(outh5name)
        except:
            try:
                subprocess.call('rm ' + outh5name)
            except:
                print('ERROR: os.remove failed, and `rm` not found on this system. Please '
                      'remove %s manually' % outh5name)
                return

    # do some string enchantments
    if '/' in outh5name or '\\' in outh5name:
        h5path = outh5name
        # take `\` and `/` out of the path
        while '/' in h5path:
            h5path = h5path[h5path.rfind('/') + 1:]
        while '\\' in h5path:
            h5path = h5path[h5path.rfind('\\') + 1:]
    else:
        h5path = outh5name
    h5path = h5path.replace(' ', '')

    # write out all the data in the dataframe if no variables are passed
    if variables is not None:
        if isinstance(variables, str):
            variables = [variables]
        data = data[variables]

    # write out each variable to the hvtk data storage
    for var in data.columns:
        variable = var.lower().replace(' ', '')
        varpath = 'data/%s' % variable
        write_h5(data[var].values, outh5name, h5path=varpath)
        varpath += '/data'
        attributes += temp_attr.format(nx=griddef.nx, ny=griddef.ny, nz=griddef.nz,
                                       h5flname=h5path, varname=variable,
                                       attrvarpath=varpath) + '\n'

    # write out the XDMF pointing to the correct attributes in the hvtk
    xmlstr_nx = griddef.nx + 1
    xmlstr_ny = griddef.ny + 1
    xmlstr_nz = griddef.nz + 1
    xmlstr_xmn = griddef.xmn - (griddef.xsiz / 2)
    xmlstr_ymn = griddef.ymn - (griddef.ysiz / 2)
    xmlstr_zmn = griddef.zmn - (griddef.zsiz / 2)
    xmlstr_xsiz = griddef.xsiz
    xmlstr_ysiz = griddef.ysiz
    xmlstr_zsiz = griddef.zsiz

    xmlstr = xmlstr.format(attributes=attributes, nx=xmlstr_nx, ny=xmlstr_ny, nz=xmlstr_nz,
                           xmn=xmlstr_xmn, ymn=xmlstr_ymn, zmn=xmlstr_zmn, xsiz=xmlstr_xsiz,
                           ysiz=xmlstr_ysiz, zsiz=xmlstr_zsiz, dim=3)

    with open(flname, 'w') as fh:
        fh.write(xmlstr)


def file_nlines(flname):
    '''Open a file and get the total number of lines. Seems pretty fast. Copied from stackoverflow
    http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python

    Parameters:
        flname (str): Name of the file to read
    '''
    f = open(flname)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    return lines


def writeout_gslib_gmm(gmm, outfile):
    """
    Writeout a fitted Gaussian mixture to the format consistent with ``gmmfit`` from the CCG
    Knowledge Base. Assume ``gmm`` is a an ``sklearn.mixture.GaussianMixture`` class fitted to data

    Note:
        Recently GMM was replaced with GaussianMixture, and there are subtle differences in
        attributes between the different versions..

    Parameters:
        gmm (GaussianMixture): a fitted mixture model
        outfile (str): the output file

    """
    if not gmm.converged_:
        print('GMM has not converged! ')
        return
    gmmstr = ''
    nvar = np.size(gmm.means_, 1)
    # write the header:
    gmmstr += 'sklearn.mixture.GMM Model for CCG gmm programs \n'
    gmmstr += '%i  %i \n' % (gmm.n_components, nvar)
    for i in range(gmm.n_components):
        # weights on one line
        gmmstr += '%.16g  %.16g \n' % (0.0, gmm.weights_[i])
        # means on one line
        gmmstr += ' '.join(['%.16g' % val for val in gmm.means_[i, :]])
        gmmstr += '\n'
        # covariance on another line
        for j in range(nvar):
            for k in range(j, nvar):
                gmmstr += '%.16g ' % gmm.covariances_[i, j, k]
        gmmstr += '\n'
    with open(outfile, 'w') as file:
        file.write(gmmstr)
    file.close()


def readvarg(vargflname, vargnum):
    """
    Reads in a variogram file and returns a pygeostat data file or list of pygeostat DataFiles with
    the variograms

    The parameter `vargnum` can be one of the following:

        * a single variogram number to return one pygeostat DataFile
        * a list of variogram numbers to return a list of DataFiles
        * the string 'all' to return a list with all variograms

    The returned `gs.DataFrame` objects contain the following parameters:

        * data['Number'] - lag number
        * data['Distance'] - lag distance
        * data['Value'] - variogram value
        * data['Points'] - number of points
        * data['Head'] - head value of variogram
        * data['Tail'] - tail value of variogram
        * tailname - tail variable name
        * headname - head variable name
        * direction - direction number

    or for gamv2004 variograms:

        * azm - azimuth
        * azmtol - azimuth tolerance
        * azmbandw - azimuth bandwidth
        * dip - dip
        * diptol - dip tolerance
        * dipbandw - dip bandwidth
        * lags - number of lags
        * lagdist - lag distance
        * lagtol - lag tolerance
        * vartype - variogram type
        * tailvar - tail variable number
        * headvar - head variable number
        * indcatcut - indicator category or cutoff (None if missing)

    Parameters:
        vargflname: name of variogram file in either gamv or gamv2004 style
        vargnum: permissible vargnum parameter as specified above

    Returns:
        vargs: A gs.DataFile (or list of DataFiles) with data set to the variogram data as
            as specified above
    """
    import io
    from . data import DataFile

    # What variogram number/numbers are we reading:
    if isinstance(vargnum, str):
        if vargnum.lower() == 'all':
            # Count the number of variograms
            vargnum = 1
            with open(vargflname, 'r') as vargfl:
                for line in vargfl:
                    if line[:1].isalpha():
                        vargnum += 1
            vargnum = list(range(1, vargnum))
        else:
            raise Exception('Invalid vargnum', vargnum)
    # Convert vargnum to list if just reading a single variogram
    if isinstance(vargnum, int):
        vargnum = [vargnum]
    # Read in each variogram
    vargidx = 0
    vargs = []
    with open(vargflname, 'r') as vargfl:
        savedline = None
        for line in vargfl:
            if not (line[:1].isalpha()) and savedline is None:
                continue
            # Increment the variogram index number
            vargidx += 1
            # Are we reading this variogram?
            if vargidx not in vargnum:
                savedline = None
                continue
            varg = DataFile(flname=vargflname, readfl=False, fltype='gslib')
            varg.vartype = '-1'
            # Get the head variable, tail variable and direction
            if savedline is None:
                headerline = line
            else:
                headerline = savedline
            # Try and get header data, might not work for vmodel outputs
            try:
                varg.tailname = headerline.partition('tail:')[2].partition('head:')[0].strip()
                if 'direction:' in headerline:
                    varg.headname = headerline.partition('head:')[2].partition('direction:')[0].strip()
                else:
                    varg.headname = headerline.partition('head:')[2].partition('direction')[0].strip()
                varg.direction = int(headerline.partition('direction')[2])
            except ValueError:
                pass
            if savedline is None:
                line = next(vargfl)
            # Check for gamv2004 output type
            if line.startswith('-HDIR'):
                # Read in gamv2004 specific information
                varg.azm = line.split()[1]
                varg.azmtol = line.split()[2]
                varg.azmbandw = line.split()[3]
                line = next(vargfl)
                varg.dip = line.split()[1]
                varg.diptol = line.split()[2]
                varg.dipbandw = line.split()[3]
                line = next(vargfl)
                varg.lags = line.split()[1]
                varg.lagdist = line.split()[2]
                varg.lagtol = line.split()[3]
                line = next(vargfl)
                varg.vartype = line.split()[1]
                varg.tailvar = line.split()[2]
                varg.headvar = line.split()[3]
                try:
                    varg.indcatcut = line.split()[3]
                except IndexError:
                    varg.indcatcut = None
                vargstr = next(vargfl)
            else:
                vargstr = line
            while True:
                line = None
                # Try and read the next line
                try:
                    line = next(vargfl)
                except StopIteration:
                    break
                # Add this line to the variogram string
                if line is not None:
                    if not line[:1].isalpha():
                        vargstr = vargstr + line
                    else:
                        savedline = line
                        break
            iovarg = io.StringIO(vargstr)
            # Read in the data
            if varg.vartype.strip() != '4':
                varg.data = pd.read_csv(iovarg, header=None, delim_whitespace=True,
                                        skipinitialspace=True,
                                        names=['Number', 'Distance', 'Value', 'Points',
                                               'Head', 'Tail'],
                                        engine='python')
            else:
                varg.data = pd.read_csv(iovarg, header=None, delim_whitespace=True,
                                        skipinitialspace=True,
                                        names=['Number', 'Distance', 'Value', 'Points',
                                               'Head', 'Tail', 'HeadVar', 'TailVar'],
                                        engine='python')
            vargs.append(varg)
    # Return with variogram(s)
    if len(vargs) == 1:
        return vargs[0]
    else:
        return vargs
