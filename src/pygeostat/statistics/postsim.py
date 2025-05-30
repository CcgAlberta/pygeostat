#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A postsim style utility
"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import pandas as pd
import numpy as np

from ..data import DataFile
from ..data import iotools as iotools

def postsim_multfiles(file_base_or_list, output_name, Nr=None, file_ending=None, fltype=None,
                      output_fltype=None, zero_padding=0, variables=None, var_min=None):
    '''The multiple file postsim function uses recursive statistics for memory management and
    coolness factor. See http://people.revoledu.com/kardi/tutorial/RecursiveStatistic/
    This function will take multiple realizations and post process the results into mean and
    variance for each variable. You can either pass it a list of files to iterate through or a
    filebase name and the number of realizations.

    Parameters:
        file_base_or_list (list) or (str): List of files or path + base name of sequentially named
            files
        output_name (str): ath (or name) of file to write output to.
        Nr (int): Number of realizations. Needed if file base name is passed.
        file_ending (str): file ending (ex. `"out"`). Used if file base name is passed. Period is
            not included.
        fltype (str): Type of data file: either ``csv``, ``gslib``, ``hdf5``, or ``gsb``. Used if
            file base name is passed and `file_ending` is not used.
        output_fltype (str): Type of output data file: either ``csv``, ``gslib``, ``hdf5``, or
            ``gsb``.
        zero_padding (int): Number of zeros to padd number in sequentially named files with. Default
            is 0.
        variables (str): List of variables to process.
        var_min (list) or (float): Minimum trimming limit to use. If one value is passed it will
            apply the trimming limit to all variables. Or a list of trimming limit for each variable
            can be passed.
    '''

    if isinstance(file_base_or_list, list):
        N = 0
        for filename in file_base_or_list:
            N += 1
            dt = DataFile(flname=filename)

            # Create the post sim file in the first call
            if N == 1:
                # setup index length and variables
                blk_count = len(dt.data)
                if variables:
                    if variables not in dt.data.columns.tolist():
                        raise KeyError('Variables passed do not match columns in datafile')
                else:
                    variables = dt.data.columns.tolist()
                columns = list()
                for var in variables:
                    columns.append(var + '_mean')
                    columns.append(var + '_variance')
                columns.append('Nr')
                # create pandas file
                postsim = pd.DataFrame(index=np.arange(blk_count), columns=columns)
                # Set nan values
                if var_min:
                    if isinstance(var_min, list):
                        if len(var_min) != len(variables):
                            raise KeyError('length of var_min list does not equal number of'
                                           ' variables being processed')
                        else:
                            trim_dict = dict(zip(variables, var_min))
                            for var in variables:
                                dt.data.setnan(variables=var, tmin=trim_dict[var])
                    elif isinstance(var_min, (int, float)):
                        dt.data.setnan(variables=variables, tmin=var_min)
                    else:
                        raise KeyError('var_min must be either a list or a number')
                # Initialize the first loop
                postsim.Nr = N
                for var in variables:
                    col_m = var + '_mean'
                    col_v = var + '_variance'
                    postsim[col_m] = dt.data[var]
                    postsim[col_v] = 0
                continue
            # ----Continue with next files----
            # Set nan values
            if var_min:
                if isinstance(var_min, list):
                    if len(var_min) != len(variables):
                        raise KeyError('length of var_min list does not equal number of'
                                       ' variables being processed')
                    else:
                        trim_dict = dict(zip(variables, var_min))
                        for var in variables:
                            dt.data.setnan(variables=var, tmin=trim_dict[var])
                elif isinstance(var_min, (int, float)):
                    dt.data.setnan(variables=variables, tmin=var_min)
                else:
                    raise KeyError('var_min must be either a list or a number')
            # calculate the stats
            postsim.Nr = N
            for var in variables:
                col_m = var + '_mean'
                col_v = var + '_variance'
                # Calculate the arithmetic average
                postsim.left_arg = ((N - 1) / N) * postsim[col_m]
                postsim.right_arg = (1 / N) * dt.data[var]
                postsim[col_m] = postsim.left_arg + postsim.right_arg
                # Calculate the variance
                postsim.left_arg = ((N - 1) / N) * postsim[col_v]
                postsim.temp = dt.data[var]
                postsim.right_arg = postsim.temp - postsim[col_m]
                postsim[col_v] = (postsim.left_arg + (1 / (N - 1)) *
                                  postsim.right_arg * postsim.right_arg)
    elif isinstance(file_base_or_list, str):
        for N in range(1, Nr + 1):
            # Create the post sim file in the first call
            if N == 1:
                if file_ending is None:
                    if fltype.lower() == 'gslib':
                        file_ending = 'out'
                    elif fltype.lower() == 'csv':
                        file_ending = 'csv'
                    elif fltype.lower() == 'gsb':
                        file_ending = 'gsb'
                    elif fltype.lower() == 'h5' or fltype.lower() == 'hdf5':
                        file_ending = 'h5'
                    else:
                        raise KeyError('Either file_ending or fltype is needed')
                filename = '{base}{number}.{ending}'.format(base=file_base_or_list,
                                                            number=str(N).zfill(zero_padding),
                                                            ending=file_ending)
                dt = DataFile(flname=filename)
                # Figure out the index (variable names... and total number of blocks)
                blk_count = len(dt.data)
                if variables:
                    if variables not in dt.data.columns.tolist():
                        raise KeyError('Variables passed do not match columns in datafile')
                else:
                    variables = dt.data.columns.tolist()
                columns = list()
                for var in variables:
                    columns.append(var + '_mean')
                    columns.append(var + '_variance')
                columns.append('Nr')
                columns.append('left_arg')
                columns.append('right_arg')
                columns.append('temp')
                # create pandas dataframe
                postsim = pd.DataFrame(index=np.arange(blk_count), columns=columns)
                # Set nan values
                if var_min:
                    if isinstance(var_min, list):
                        if len(var_min) != len(variables):
                            raise KeyError('length of var_min list does not equal number of'
                                           ' variables being processed')
                        else:
                            trim_dict = dict(zip(variables, var_min))
                            for var in variables:
                                dt.data.setnan(variables=var, tmin=trim_dict[var])
                    elif isinstance(var_min, (int, float)):
                        dt.data.setnan(variables=variables, tmin=var_min)
                    else:
                        raise KeyError('var_min must be either a list or a number')
                # Initialize the first loop
                postsim.Nr = N
                for var in variables:
                    col_m = var + '_mean'
                    col_v = var + '_variance'
                    postsim[col_m] = dt.data[var]
                    postsim[col_v] = 0
                continue
            # ----Continue with next files----
            # open the next file
            filename = '{base}{number}.{ending}'.format(base=file_base_or_list,
                                                        number=str(N).zfill(zero_padding),
                                                        ending=file_ending)
            dt = DataFile(flname=filename)
            # Set nan values
            if var_min:
                if isinstance(var_min, list):
                    if len(var_min) != len(variables):
                        raise KeyError('length of var_min list does not equal number of'
                                       ' variables being processed')
                    else:
                        trim_dict = dict(zip(variables, var_min))
                        for var in variables:
                            dt.data.setnan(variables=var, tmin=trim_dict[var])
                elif isinstance(var_min, (int, float)):
                    dt.data.setnan(variables=variables, tmin=var_min)
                else:
                    raise KeyError('var_min must be either a list or a number')
            # calculate the stats
            postsim.Nr = N
            for var in variables:
                col_m = var + '_mean'
                col_v = var + '_variance'
                # Calculate the arithmetic average
                postsim.left_arg = ((N - 1) / N) * postsim[col_m]
                postsim.right_arg = (1 / N) * dt.data[var]
                postsim[col_m] = postsim.left_arg + postsim.right_arg
                # Calculate the variance
                postsim.left_arg = ((N - 1) / N) * postsim[col_v]
                postsim.temp = dt.data[var]
                postsim.right_arg = postsim.temp - postsim[col_m]
                postsim[col_v] = (postsim.left_arg + (1 / (N - 1)) *
                                  postsim.right_arg * postsim.right_arg)
    else:

        raise TypeError('file_base_or_list must be either a list or a string')
    # Write out the results
    columns.remove('left_arg')
    columns.remove('right_arg')
    columns.remove('temp')
    if output_fltype is None:
        if fltype:
            output_fltype = fltype
        elif file_ending == '.out':
            output_fltype = 'gslib'
        elif file_ending == 'csv' or file_ending == 'gsb' or file_ending == 'h5':
            output_fltype = file_ending
        else:
            raise KeyError('Unable to figure out what output_fltype needs to be. Please pass'
                           ' a value')

    if output_fltype.lower() == 'gslib':
        iotools.write_gslib(postsim, output_name, variables=columns)
    elif output_fltype.lower() == 'csv':
        iotools.write_csv(postsim, output_name, variables=columns)
    elif output_fltype.lower() == 'gsb':
        iotools.write_gsb(postsim, output_name, tvar='Nr', variables=columns)
    elif output_fltype.lower() == 'h5' or fltype.lower() == 'hdf5':
        iotools.write_h5(postsim, output_name, variables=columns)
    else:
        raise NotImplementedError('output_fltype did not match any of the implemented filetypes')
