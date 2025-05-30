#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''h5_pytables.py: Old HDF5 funcionality'''

#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import pandas as pd


class DataFile():
    """
    === HDF5 ===

    Using the HDF5 file format has its own possitive features. For one it reads
    and writes much faster then using the ASCII format. Also attributes
    (like the grid definition) can be saved with the file and pygeostat will
    autoload those attributes. All files for a single project can also be saved
    in the same file as separate tables and each table or subset of a table can
    be read in as needed.
    HDF5 file simple read example

    >>> datafl = gs.DataFile(flname='../data/oilsands_out.hdf5', dftype='hdf5',
    >>>                      hdf_key='gslib_data')

    To view the HDF5 header information (tables stored in the file)

    >>> datafl.store

    If you have a HDF5 file with multiple tables and you just want to read in the file information
    to view what tables are in the file and any attributes saved to the file you can do a header
    style only read

    >>> datafl = gs.DataFile(flname='../data/oilsands_out.hdf5', dftype='hdf5',
    >>>                      headeronly=True, hdf5_header=True)

    Then to see what tables are written in the hdf5 file

    >>> datafl.store
    """



    def retrieve_hdf5(self, hdf_key):
        '''Retrieve a new group from self.store and load it into the self.data location
        of the gs.DataFile()

        Note:
            This currently won't work if you don't have a self.store in the DataFile.
            You need to open a .h5 file with the header option in order to use this.

        Parameters:
            hdf_key (str): name of the group to load into self.data

        Examples:

        load the initial data file (say the 1st realization)

        >>> datafl = gs.DataFile(flname='sgsim.out', hdf5_header=True, hdf_key='Realization_1')
        >>> # Look at some summary statistics of the first realization
        >>> datafl.describe()

        Later you decide you want to look at the 4th realization so now you retrieve that data

        >>> datafl.retrieve_hdf5('Realization_4')
        >>> # Look at the new summary statistics for the 4th realization
        >>> datafl.describe()

        '''

        if self.store is None:
            raise Exception('There is no self.store - Please read in an hdf5 file with the header',
                            ' option')
        else:
            from ..data import readfile as readfile
            self.data = readfile(self.flname, fltype='hdf5', hdf_key=hdf_key)

        return self

        from .data import open_h5 as open_h5
        # read the file and import data and/or header
        if headeronly is True:
            if hdf5_header is True:
                self.store = open_h5(flname)
                self.data = pd.DataFrame()
            elif hdf5_header is False:
                self.data = readfile(self.flname, fltype=fltype, headeronly=headeronly,
                                     delimiter=delimiter, hdf_key=hdf_key, columns=columns)
        elif headeronly is False:
            if hdf5_header is True:
                self.data = readfile(self.flname, fltype=fltype, headeronly=headeronly,
                                     delimiter=delimiter, hdf_key=hdf_key, columns=columns)
                self.store = open_h5(flname)
            elif hdf5_header is False:
                self.data = readfile(self.flname, fltype=fltype, headeronly=headeronly,
                                     delimiter=delimiter, hdf_key=hdf_key, columns=columns)


        # If a hdf5 file was read in check to see if any file attributes were loaded
        if hasattr(self, 'store'):
            if hasattr(self.store.root._v_attrs, 'gridstr'):
                gridstr = self.store.root._v_attrs.gridstr
                from .data import GridDef as GridDef
                if self.griddef is not None:
                    print('Warning: Grid definition is being rewritten by griddef saved in h5 file')
                if gridstr == 'None':
                    print('Warning: unable to set grid definition due to empty grid string saved',
                          ' in h5 file')
                else:
                    self.griddef = GridDef(gridstr)
            # set attributes call

            self.attrs = self.store.root._v_attrs

def write_h5_pytables(data, flname, variables=None, hdf_key='gslib_data', mode='w', hdf_append=False, complevel=None, complib=None, griddef=None):
    """Writes out a HDF5 binary data file. This is useful for greatly increasing
    the read and write speads to memory which is a common bottleneck. There are
    some usage notes and warnings that should be reviewed.

    For more information on HDF5 - search youtube for some of the python hdf5 review videos.

    For a HDF5 file viewing software (outside of python) with a GUI interface check out
    the HDFView software at http://www.hdfgroup.org

    Notes:
        For out of the box best speed - use the default of no compression.

        If you want faster speed then reading and writing text (ie. gslib or csv) but
        also have some file size concerns then the best current out of the box option
        is using the "zlib" compression library.

        The "blosc" compression library is supposed to have both faster read/write speeds
        than no compression and it provides some compression. It produces faster read/write
        speeds by reorganizing the structure of the binary file in the same way the computer
        memory stores information. *Warning* Currently the "blosc" compression library used
        is the one packaged in pytables. The current pytables build (3.2.0) causes python to
        crash when you try to use the "blosc" compression. Until a fix comes out downgrade
        pytables to version 3.1.1
        This can be done with conda using the following code in your command line of choice

        >>> conda install pytables=3.1.1

        Then saying yes to the downgrade warnings.

    Parameters:
        data (pygeostat.DataFile or pandas.DataFrame): data to write out
        flname (str): Path (or name) of file to write out.

    Keyword Args:
        variables (List(str)): List of variables to write out if only a subset is desired.
        hdf_key (str): Name of the group to write the data to in hdf5 file. If appending then
            multiple "groups" can be written to the same file.
        mode (str): write mode. `r` read only, `w` write a new file (existing file will be
            overwritten), `a` append (if no file exists a new file will be created), `r+ similar
            to `a` but file must already exist
        hdf_append (bool): For Table formats, append the input data to the existing
        complevel (int): Compression level from 0 (least) to 9 (greatest). Default: 9
        complib (str): Compression Library. None, "zlib", or "blosc". See Note for warning.
            Default: "blosc"
        chuncksize (int): Write out by chunks using the chunksize. This can significantly
            lower your memory usage in writing.

    Things To Look Into Adding:

        1. In the future we might add the  keyword argument "expectedrows"
        which will allow pandas/pytables to optimize the read/write performance.

        2. Adding the attributes (like grid definitions!) to the HDF file

    To write out a new file

    >>> gs.write_hdf5(data, flname, hdf_key='group name')

    To append data to an already existing hdf5 file

    >>> gs.write_hdf5(data, flname, hdf_key='group name', mode='a')

    To append rows to an existing table in a hdf5 file

    >>> gs.write_hdf5(data, flname, hdf_key'group name in file to append to', mode='a',
    >>>               hdf_append=True)

    """
    from .data import DataFile as DataFile
    # If a DataFile is used, check the arguments
    if complevel is None:
        complevel = 9
    if complib is None:
        complib = "blosc"
    if isinstance(data, DataFile):
        if variables is None:
            variables = data.data.columns.tolist()
    else:
        if variables is None:
            variables = data.columns.tolist()
    # Force a write
    tempstore = pd.HDFStore(flname, mode=mode)
    if isinstance(data, DataFile):
        tempstore.append(hdf_key, data.data[variables], index=True, DataColumns=True,
                         complib=complib, complevel=complevel, append=hdf_append)
    else:
        tempstore.append(hdf_key, data[variables], index=True, DataColumns=True,
                         complib=complib, complevel=complevel, append=hdf_append)
    # if a griddef is available or passed then append the info to the file

    if griddef is not None:
        tempstore.root._v_attrs.gridstr = str(griddef)
    elif hasattr(data, 'griddef'):
        if hasattr(data, 'store'):
            try:
                data.store.root._v_attrs.gridstr
            except AttributeError:
                print('Warning: overwritting griddef in h5 file')
                tempstore.root._v_attrs.gridstr = str(data.griddef)
        else:
            tempstore.root._v_attrs.gridstr = str(data.griddef)

    tempstore.close()

def open_h5(flname):
    '''Open a hdf5 without reading any of the data into memory.

    Parameters:
        flname (str): Path (or name) of file to read.

    Returns:
        store (pandas.pytables.HDFStore): Pytables HDFStore object

    Quickly read in a HDF5 file and look at the available tables

    >>> store = gs.open_hdf5(flname)
    >>> store

    To create a DataFile object with no data and then review what tables are stored in the file

    >>> datafl = gs.DataFile(flname=filename, headeronly=True, hdf5_header=True)
    >>> datafl.store

    To load a table into the DataFile as well as the store info

    >>> datafl = gs.DataFile(flname=filename, hdf5_header=True)
    '''
    store = pd.HDFStore(flname)

    return store


def read_h5_pytables(flname, hdf_key='gslib_data'):
    '''Reads in a GSLIB-style HDF5 binary data file. This is useful for greatly increasing
    the read and write speads to memory which is a common bottleneck. There are
    some usage notes and warnings that should be reviewed.

    For more information on HDF5 - search youtube for some of the python hdf5 review videos.

    For a HDF5 file viewing software (outside of python) with a GUI interface check out
    the HDFView software at http://www.hdfgroup.org

    Parameters:
        flname (str): Path (or name) of file to read.

    Keyword Args:
        hdf_key (str): The name of the group in the file which the data is stored under

    Returns:
        data (pandas.DataFrame): Pandas DataFrame object with input data.

    Note:
        To open the file without reading any data into memory use ``gs.open_hdf5``
        This function uses the builtin functionality of pandas which in turn utilizes Pytables.
        Pytables is row oriented so it seems like you can't append more column to a current table
        once in the h5 file. You can append more rows though. You would need to overwrite the
        and create a new table if you want more columns.
        For some great discussion on large data workflows check out this stack overflow question.
        http://stackoverflow.com/questions/14262433/large-data-work-flows-using-pandas
    '''
    # Can the file be opened?
    try:
        with open(flname, 'r'):
            pass
    except FileNotFoundError:
        raise Exception(' '.join([__name__, 'File could not be opened:', flname]))
    # Read in the data
    data = pd.read_hdf(flname, key=hdf_key)
    # varnames = data.columns.tolist()

    # Only return the DataFrame
    return data


def convert_rlzns_to_h5(flname, writeflname, griddef=None, num_real=0, real_size=0, fltype=None, complevel=9, complib='blosc', label_padding=3, subfolder='BT_Units/', label_prefix='Realization_', quiet=False):
    """
    A function to convert a text based file of multiple realizations to the HDF5 format. Each
    realization will be saved under the model directory in the specified subfolder as a seperate
    table in the file.

    Parameters:
        flname (str): Name of the realization file to convert to hdf5
        writeflname (str): File name for output hdf5 file
        griddef (obj): pygeostat.GridDef grid definition
        num_real (int): Number of realizations in the file
        real_size (int): Length of each realization
        fltype (str): Type of file being read `gslib` or `csv`. default is `gslib`
        complevel (0-9): Level of compression in the hdf5 file. 0 = None, 9 = Most compressed
        comlib (str): compression library to use `blosc`, `zlib` etc. Can use any pytables
            compression library
        label_padding (int): Zero padding for creating the labels. By default it is set to three
            so the first realization will be labeled as `Realization_001`
        label_prefix (str): The group name and file structure for storing in the H5 file. By default
            this is set up as `Models/Realization_xxx`

    Note:
        Needs either a grid definiton or number of realizations.

        If number of realizations and realization size are not passed then it will perform a quick
        read of the file to determine how many rows of data are in the file and use this to figure
        out missing information

    Examples:
        A simple fast call:

        >>> gs.convert_rlzns_to_h5('sgsim.out', 'project.h5', griddef=grid, num_real=100)
        >>> '100 Realizations written to project.h5'
        >>> 'Realizations written to the following h5 directory: Models/BT_Units/'
    """
    import warnings
    from .datautils.labels.make_labels import make_labels as make_labels
    # ------------------------------------------------
    # Sanity checks
    # ------------------------------------------------
    # Check to make sure either a griddef is give or the number of realizations
    if griddef is not None:
        real_size = griddef.count()
        gridstr = str(griddef)
    if griddef is None:
        if num_real is 0:
            raise Exception('ERROR: Needs either a grid definition or the number of realizations')

    # If no filetype is given guess at what type of file
    if fltype is None:
        # Guess a file type based on the file name if obvious
        if flname.endswith('.out') or flname.endswith('.dat'):
            fltype = 'GSLIB'
        elif flname.endswith('.csv'):
            fltype = 'CSV'
        elif flname.endswith('.h5') or flname.endswith('.hdf5') or flname.endswith('.hd5'):
            fltype = 'HDF5'
        else:
            # Otherwise just assume GSLIB for now but this may change
            fltype = 'GSLIB'

    # Raise an exception if the file is already a hdf5 file
    if fltype.lower() == 'hdf5':
        raise Exception('Not currently set up to handle hdf5 as INPUT file. Only as OUTPUT file')

    # ------------------------------------------------
    # Process header information from the files
    # ------------------------------------------------
    # If gslib style then skip the header info lines
    if fltype.lower() == 'gslib':
        with open(flname, 'r') as datafl:
            # Geo-EAS Header
            _ = datafl.readline().strip()  # Title
            nvar = int(datafl.readline().split()[0])
            varnames = []
            for _ in range(nvar):
                varnames.append(datafl.readline().strip())
            nrowtoskip = 2 + nvar
            header_lines = nrowtoskip
        delimiter = r'\s*'
        header = None

    # Set read options if csv type file
    if fltype.lower() == 'csv':
        temp = pd.read_csv(flname, header=0, nrows=1)
        varnames = temp.columns.tolist()
        del temp
        nrowtoskip = 0
        header_lines = 1
        delimiter = ','
        header = 0

    # ------------------------------------------------
    # Set initial options and find any missing info
    # ------------------------------------------------

    # If either real_size or num_real is missing figure out what they are
    if real_size is 0:
        lines = iotools.file_nlines(flname)
        data_size = lines - header_lines
        real_size = int(data_size / num_real)
    if num_real is 0:
        lines = iotools.file_nlines(flname)
        data_size = lines - header_lines
        num_real = int(data_size / real_size)
    # a couple of random variables including a list of group names
    len_var = len(varnames)
    # This allows us to be dumb and will still work if we forget how many realizations we have
    if num_real < 1000:
        num_labels = 1000
    else:
        num_labels = num_real
        if label_padding < 4:  # increase the padding if we have more than 1000 realizations
            label_padding = 4
    table_name = make_labels(label_prefix, num_labels, padding=label_padding)
    idx = 0
    # ------------------------------------------------
    # Iterratively write to the hdf5 file
    # ------------------------------------------------

    # Open a HDFfile to write to. Note will overwrite any existing file
    tempstore = pd.HDFStore(writeflname, 'w')

    # Because we can.... Write some grid information to the hdf5 file
    tempstore.root._v_attrs.realsize = real_size
    tempstore.root._v_attrs.numreal = num_real
    if gridstr:
        tempstore.root._v_attrs.gridstr = gridstr

    # Iterate through the read file chunk by chunk and write it to a hdf5 file
    for chunk in pd.read_csv(flname, skiprows=nrowtoskip, skipinitialspace=True,
                             chunksize=real_size, header=header, engine='python',
                             delimiter=delimiter):
        # Check to see if any trailing spaces got read in as a new column
        len_ch_col = len(chunk.columns)
        if len_ch_col > len_var:
            didx = len(varnames)
            del chunk[didx]
            warnings.warn('There is probably trailing spaces in the GSLIB file. Found ', len_ch_col,
                          'columns, but there is header information for only ', len_var,
                          'columns. Last column was assumed to be trailing spaces and '
                          'therefore was deleted')
        # Set column names and write to h5 file
        chunk.columns = varnames
        tempstore.append(table_name[idx], chunk, header=False, index=True, DataColumns=True,
                         complevel=complevel, complib=complib)
        idx += 1

    # Check to make sure we wrote the correct number of realizations to file
    tables_appended = idx - 1
    if tables_appended != num_real:
        print('WARNING: Number of realizations =', num_real,
              ' But number of realizations appended to h5 file =', tables_appended)

    # ------------------------------------------------
    # Clean up after ourselves
    # ------------------------------------------------
    tempstore.close()

    if quiet is False:
        print((idx - 1), ' Realizations written to ', flname)
        print('Realizations written to the following h5 directory: Models/' + subfolder)
