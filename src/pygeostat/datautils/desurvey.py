#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''desurvey.py: Contains utilities for desurveying drill hole data'''

__author__ = 'pygeostat development team'
__date__ = '2014-04-03'
__version__ = '3.000'

from math import sin, cos, radians, ceil
import numpy as np
import pandas as pd
from ..data import DataFile
import warnings
from collections import OrderedDict
import xml.etree.cElementTree as ET
import copy
from .. pygeostat_parameters import Parameters


class Drillhole(object):
    '''This class contains specific drill hole data and metadata.

    Drill hole classes may be created or generated using pygeostat functions.
    This is primarily used for desurveying drill hole data.

    Parameters:
        holeid (str): Hole Id
        collarx (numeric): x coordinate
        collary (numeric): y coordinate
        collarz (numeric): z coordinate

    .. codeauthor: pygeostat development team 2014-04-03'''

    def __init__(self, holeid=None, collarx=None, collary=None, collarz=None):
        self.holeid = holeid
        self.collarx = float(collarx)
        self.collary = float(collary)
        self.collarz = float(collarz)
        self.data = pd.DataFrame()

    def __str__(self):
        return('Drillhole {} collared at ({},{},{})'.format(self.holeid, self.collarx,
                                                            self.collary, self.collarz))

    def delx(self, dist, inclination, bearing):
        '''Calculate change in X

        Parameters:
            self
            dist (float): distance along direction
            inclination (float)
            bearing (float)

        Returns:
            Value (float): A float value of the change in distance

        .. codeauthor: pygeostat development team 2014-04-03'''
        val = dist * cos(radians(inclination)) * sin(radians(bearing))
        return val

    def dely(self, dist, inclination, bearing):
        '''Calculate change in y

        Parameters:
            self:
            dist (float): distance along direction
            inclination (float)
            bearing (float)

        Returns:
            Value (float): A float value of the change in distance

        .. codeauthor: pygeostat development team 2014-04-03'''
        val = dist * cos(radians(inclination)) * cos(radians(bearing))
        return val

    def delz(self, dist, inclination, bearing):
        '''Calculate change in z

        Parameters:
            self
            dist (float): distance along direction
            inclination (float)
            bearing (float): Not currently used in the calculation

        Returns:
            Value (float): A float value of the change in distance

        .. codeauthor: pygeostat development team 2014-04-03'''
        val = dist * sin(radians(inclination))
        return val

    def getxyz(self, along, null):
        '''Return X, Y, Z values given the length along the drill hole

        Parameters:
            self:
            along (numeric): Distance along the drill hole

        Returns:
            x, y, z (float)

        .. codeauthor: pygeostat development team 2014-04-03'''
        x, y, z = (self.collarx, self.collary, self.collarz)
        # Check some base cases
        if (along < 0.) or (along > float(self.data['Along'][-1:])):
            x, y, z = (null, null, null)
            print('ERROR: This along is invalid ', along, 'for', self.holeid)
        elif along == 0.:
            pass
        else:
            loc = 0.
            idx = 0
            while loc < along:
                try:
                    alongcheck = self.data['Along'][idx + 1]
                except KeyError:
                    alongcheck = 1e21
                    print('WARNING: Hole {} only has 1 survey station, duplicating it'.format(self.holeid))
                    idx = idx - 1
                if along >= alongcheck:
                    dist = self.data['Along'][idx + 1] - loc
                    loc = self.data['Along'][idx + 1]
                    idx += 1
                    x += self.delx(dist, self.data['Inclination'][idx],
                                   self.data['Azimuth'][idx])
                    y += self.dely(dist, self.data['Inclination'][idx],
                                   self.data['Azimuth'][idx])
                    z += self.delz(dist, self.data['Inclination'][idx],
                                   self.data['Azimuth'][idx])
                else:
                    dist = along - loc
                    loc = along
                    idx += 1
                    x += self.delx(dist, self.data['Inclination'][idx],
                                   self.data['Azimuth'][idx])
                    y += self.dely(dist, self.data['Inclination'][idx],
                                   self.data['Azimuth'][idx])
                    z += self.delz(dist, self.data['Inclination'][idx],
                                   self.data['Azimuth'][idx])
        return(x, y, z)


def set_desurvey(collarfl, surveyfl, along_name, azimuth_name, inclination_name):
    '''Sets up for desurveying returning a dictionary of DrillHole objects

    Parameters:
        collarfl (DataFile): pygeostat DataFile with dh, x, y, z all set
        surveyfl (DataFile): pygeostat DataFile with dh set
        along_name (str): variable name of "along" in survey file
        azimuth_name (str): variable name of "azimuth/bearing" in survey file
        inclination_name (str): variable name of "inclination" in survey file

    Returns:
        drillholes (dict): Returns a dictionary of DrillHole objects

    .. codeauthor:: pygeostat development team - 2014-04-03'''
    print(("\nWARNING: drillhole routines are untested with recent Pygeostat changes!\n"
           "           Proceed with caution!"))
    # Check for missing components
    for item in [collarfl.dh, collarfl.x, collarfl.y, collarfl.z,
                 surveyfl.dh, surveyfl.data[along_name],
                 surveyfl.data[azimuth_name], surveyfl.data[inclination_name]]:
        pass
    # Create drill hole objects with collar information
    drillholes = {}
    for rowidx, row in collarfl.data.iterrows():
        drillholes[row[collarfl.dh]] = Drillhole(holeid=row[collarfl.dh], collarx=row[collarfl.x],
                                                 collary=row[collarfl.y], collarz=row[collarfl.z])
    # Get along, azimuth and inclination data for this drill hole
    for dhid in drillholes:
        # Find rows with this drill hole ID
        rows = surveyfl.data[surveyfl.data[surveyfl.dh] == dhid].index.values
        # Check that rows were actually found
        if len(rows) == 0:
            raise Exception('No survey information found for' + drillholes[dhid])
        # Attach the data to the Drillhole
        drillholes[dhid].data = pd.DataFrame({'Along': surveyfl.data[along_name].ix[rows],
                                              'Azimuth': surveyfl.data[azimuth_name].ix[rows],
                                              'Inclination': surveyfl.data[inclination_name].ix[rows]})
        # Sort by along to ensure ordered values and reset the index
        drillholes[dhid].data = drillholes[dhid].data.sort_values(by='Along')\
                                                .reset_index()[['Along', 'Azimuth', 'Inclination']]
    return drillholes


def get_desurvey(datafl, drillholes, inplace=True, x='X', y='Y', z='Z', null=None):
    '''Gets the desurvey of a DataFile given froms and tos.

    Parameters:
        datafl (DataFile): pygeostat DataFile to be desurveyed with one of the options:

                        1) just ifrom or just ito set to desurvey at exactly that distance
                        2) both ifrom and ito set to desurvey halfway in between (ie: midpoint)

        drillholes (dict): dictionary of Drillhole objects, obtained from set_desurvey

    Keyword Args:
        inplace (bool): modifies the datafl.data to have X, Y and Z values
            otherwise returns a pandas dataframe with X, Y and Z

    .. codeauthor: pygeostat development team 2014-04-03'''

    # Create a temporary data file - this allows variable creation
    # without leaving behind a mess or overwriting anything accidentally!
    # There are probably cleaner ways to do this, but desurveying is normally a 1-off.
    print(("\nWARNING: drillhole routines are untested with recent Pygeostat changes!\n"
           "           Proceed with caution!"))
    if null is None:
        null = Parameters['data.null']
    tempfl = pd.DataFrame(datafl.data[datafl.dh])
    # Get the midpoint, or 'along' value to use
    if datafl.ifrom is not None and datafl.ito is not None:
        tempfl['mid'] = (datafl.data[datafl.ito] + datafl.data[datafl.ifrom]) * 0.5
    elif datafl.ifrom is not None and datafl.ito is None:
        tempfl['mid'] = datafl.data[datafl.ifrom]
    elif datafl.ifrom is None and datafl.ito is not None:
        tempfl['mid'] = datafl.data[datafl.ito]
    else:
        raise Exception('Both ifrom and ito are missing from datafl!')

    # Define a local desurvey function
    def desurvey(row):
        'Define a local desurvey function'
        # Get the corresponding drill hole
        drillhole = drillholes[row[datafl.dh]]
        # Get the X,Y,Z coordinates
        xloc, yloc, zloc = drillhole.getxyz(row['mid'], null)
        return pd.Series({x: xloc, y: yloc, z: zloc})
    # Desurvey
    xyzpoints = tempfl.apply(desurvey, axis=1)
    if not inplace:
        return xyzpoints
    else:
        # Assign data
        datafl.data[x] = xyzpoints[x]
        datafl.data[y] = xyzpoints[y]
        datafl.data[z] = xyzpoints[z]
        # Set datafl.x, datafl.y, datafl.z
        if datafl.x is None:
            datafl.x = x
        if datafl.y is None:
            datafl.y = y
        if datafl.z is None:
            datafl.z = z
        return None


def set_comps(datafl, complength):
    '''Returns pandas dataframe with drillhole, compfrom, compto

    Parameters:
        datafl (DataFile): pygeostat DataFile with dh, ifrom and ito set
        complength (numeric): length of composites

    Returns:
        comps (DataFrame): pandas DataFrame with dh, ifrom and ito set for each composite

    .. codeauthor: pygeostat development team 2014-04-03'''
    print(("\nWARNING: drillhole routines are untested with recent Pygeostat changes!\n"
           "           Proceed with caution!"))
    # Get list of drill holes
    dhids = datafl.unique_cats(datafl.data[datafl.dh])
    # Determine start (smallest 'from' value) and end (largest 'to' value) of each dh
    starts = {}
    ends = {}
    nsamples = {}
    totalsamples = 0
    for dhid in dhids:
        starts[dhid] = min(datafl.data[datafl.data[datafl.dh] == dhid][datafl.ifrom])
        ends[dhid] = max(datafl.data[datafl.data[datafl.dh] == dhid][datafl.ito])
        # Determine number of composites in each
        nsamples[dhid] = ceil((ends[dhid] - starts[dhid]) / complength)
        ends[dhid] = nsamples[dhid] * complength + starts[dhid]
        totalsamples += nsamples[dhid]
    # Initialize pandas dataframe
    datadict = OrderedDict()
    datadict[datafl.dh] = [dhid for dhid in range(int(totalsamples))]
    datadict[datafl.ifrom] = np.ones(int(totalsamples))
    datadict[datafl.ito] = np.ones(int(totalsamples))
    comps = pd.DataFrame(datadict)
    # Composites
    startidx = 0
    for dhid in dhids:
        # Pandas slices using INCLUSIVE INDEXING therefore subtract 1 and add 1 later
        endidx = startidx + nsamples[dhid] - 1
        # Set drill hole number, from and to
        comps.loc[startidx:endidx, datafl.dh] = dhid
        comps.loc[startidx:endidx, datafl.ifrom] = np.linspace(starts[dhid],
                                                               ends[dhid] - complength,
                                                               num=nsamples[dhid], endpoint=True)
        comps.loc[startidx:endidx, datafl.ito] = np.linspace(starts[dhid] + complength, ends[dhid],
                                                             num=nsamples[dhid], endpoint=True)
        startidx = endidx + 1
    return comps


def get_comps(comps, datafl, vartypes='continuous', null=None, nprocess=None):
    '''Returns a pandas DataFrame with upscaled composites parallelizing across
       drillholes.

    Parameters:
        datafl (DataFile): pygeostat DataFile with dh, ifrom, ito and at least 1 variable to upscale
        comps (DataFrame): pandas DataFrame with datafl.dh, datafl.ifrom, datafl.ito
        corresponding to the composite locations

    Keyword Args:
        vartypes (str OR dict): 'continuous', 'categorical' or a dictionary of variables like:
        {'Cu':'continuous','Facies','categorial'}
        null: the null value (values less than or equal to this will not be used)

    Returns:
        upscaled (DataFrame): Pandas DataFrame with values of upscaled variable

    .. codeauthor: pygeostat development team 2014-04-03'''
    import multiprocessing as mp
    print(("\nWARNING: drillhole routines are untested with recent Pygeostat changes!\n"
           "           Proceed with caution!"))
    if null is None:
        null = Parameters['data.null']
    # Get a list of drill holes for composites
    dhids = sorted(list(set(comps[datafl.dh])))
    # Determine the number of composites
    ncomps = len(comps[datafl.dh])
    # Create a dictionary of variable upscaling methods if not supplied
    if isinstance(vartypes, str):
        vartypedict = {}
        for var in datafl.data.columns:
            if var not in [datafl.dh, datafl.ifrom, datafl.ito]:
                vartypedict[var] = vartypes.lower()
    else:
        vartypedict = vartypes
    # Determine number of variables
    nvar = len(vartypedict)
    # Initialize a pandas DataFrame for upscaled values
    upscaled = pd.DataFrame(np.ones([ncomps, nvar]) * null, columns=list(vartypedict))
    upscaled[datafl.dh] = comps[datafl.dh]
    upscaled[datafl.ifrom] = comps[datafl.ifrom]
    upscaled[datafl.ito] = comps[datafl.ito]

    # Group data by drill hole ID
    grouped_comps = comps.groupby([datafl.dh], squeeze=True, as_index=False)
    grouped_data = datafl.data.groupby([datafl.dh], squeeze=True, as_index=False)

    # Upscale all drill holes in parallel asynchronously
    if nprocess is None:
        nprocess = Parameters['config.nprocess']
    # Create a pool for processing, calling with None to use all CPUs
    for var, vartype in vartypedict.items():
        # print('Assembling pool for {}...'.format(var))
        pool = mp.Pool(processes=nprocess)
        results = {}
        number_of_valid_dhs = 0
        for dhid, compdf in grouped_comps:
            # Get the data which corresponds to this composite
            try:
                datadf = grouped_data.get_group(dhid)
                number_of_valid_dhs += 1
            except KeyError:
                # No data for this composite drill hole so skip it
                continue
            # Add on the GSLIB call to execute asynchronously
            results[dhid] = pool.apply_async(upscale_worker, (compdf, datadf, var, vartype,
                                                              datafl.ifrom, datafl.ito, null))
        # Close the list of processes for the pool
        pool.close()
        # "Join" the list of processes to execute and wait for completion
        pool.join()
        # Tabulate the results that are not null (ie: a float)
        try:
            all_results = (pd.concat([group.get() for name, group in results.items()
                                      if type(group.get()) is not float]))
        except ValueError:
            # No data were found
            warnings.warn('No valid samples found for variable {}!!!'.format(var))
            return(upscaled)
        # Join this to our database
        upscaled[var] = upscaled[var].replace({null: np.nan})
        upscaled[var].fillna(all_results, inplace=True)
        upscaled[var].fillna(value=null, inplace=True)

    return(upscaled)


def upscale_worker(compdata, dhdata, var, vartype, ifrom, ito, null=None):
    '''Returns the upscaled value contained in [compfrom,compto]

    Parameters:
        compdata (DataFrame): pandas DataFrame with columns for ifrom, ito for the DH of interest
        dhdata (DataFrame): pandas DataFrame with columns for ifrom, ito and var for the DH
            of interest
        var (str): variable name
        vartype (str): 'continuous' or 'categorical'

    Keyword Args:
        null: value to return if no data in interval

    Returns:
        compdata[var]: pandas series of upscaled variable using either:
            continuous - length weighted linear average or
            categorical - most frequent category

    .. codeauthor: pygeostat development team 2014-04-03'''
    if null is None:
        null = Parameters['data.null']
    # Initialize the upscaled variable
    compdata[var] = np.ones(len(compdata[ifrom])) * null
    # Iterate over composite rows
    for compidx, comprow in compdata.iterrows():
        compfrom = comprow[ifrom]
        compto = comprow[ito]
        # Define local function to get length in interval

        def length_in_interval(datarow):
            'finds the length in the interval along the drill hole - tda'
            # Is the variable value valid?
            if datarow[var] <= null:
                return 0.0
            # Case 1 - entire sample interval in composite
            if (datarow[ifrom] >= compfrom) and (datarow[ito] <= compto):
                return datarow[ito] - datarow[ifrom]
            # Case 2 - partially bottom of interval
            elif ((datarow[ifrom] >= compfrom) and (datarow[ifrom] < compto) and
                  (datarow[ito] > compto)):

                return compto - datarow[ifrom]
            # Case 3 - partially in top of interval
            elif ((datarow[ifrom] < compfrom) and (datarow[ito] > compfrom) and
                  (datarow[ito] <= compto)):

                return datarow[ito] - compfrom
            # Case 4 - composite bracketed by sample interval
            elif (datarow[ifrom] <= compfrom) and (datarow[ito] >= compto):
                return compto - compfrom
            # Case 5 - no length in interval
            else:
                return 0.0
        # Get length of each sample in composite
        try:
            dhdata['Length'] = dhdata.apply(length_in_interval, axis=1)
        except ValueError as e:  # No valid data in drill hole found - null assignment!
            return null
        # If no data was found, return null value
        if sum(dhdata['Length']) <= 0.0:
            upscaled = null
        else:
            # Upscale the data
            # - Continuous
            if vartype.startswith('con'):
                upscaled = np.dot(dhdata[var], dhdata['Length']) / sum(dhdata['Length'])
            # - Categorical
            elif vartype.startswith('cat'):
                unique_cats = list(set(dhdata[dhdata['Length'] > 0.0][var]))
                unique_cat_wts = [sum(dhdata[dhdata[var] == cat]['Length']) for cat in unique_cats]
                upscaled = unique_cats[unique_cat_wts.index(max(unique_cat_wts))]
            # - Invalid variable type
            else:
                raise Exception('Invalid vartype', vartype)
        # Save this value
        compdata.loc[compidx, var] = upscaled

    # Return all values for this drill hole
    return(compdata[var])


def serial_upscale(compdh, compfrom, compto, datafl, var, vartype, null=None):
    '''Returns the upscaled value contained in [compfrom,compto]

    Parameters:
        compdh: composite drill hole
        compfrom: composite 'from'
        compto: composite 'to'
        datafl (DataFile): pygeostat DataFile with dh, ifrom, ito and var
        var (str): variable name
        vartype (str): 'continuous' or 'categorical'

    Keyword Args:
        null: value to return if no data in interval

    Returns:
        upscaled: value of upscaled variable using either:
            continuous - length weighted linear average or
            categorical - most frequent category

    .. codeauthor: pygeostat development team 2014-04-03'''
    if null is None:
        null = Parameters['data.null']
    # Get the variable, froms and tos which are in that drill hole
    dhdata = datafl.data[datafl.data[datafl.dh] == compdh][[datafl.ifrom, datafl.ito, var]]
    # Define local function to get length in interval

    def length_in_interval(datarow):
        'finds the length in the interval along the drill hole - tda'
        # Is the variable value valid?
        if datarow[var] <= null:
            return 0.0
        # Case 1 - entire sample interval in composite
        if (datarow[datafl.ifrom] >= compfrom) and (datarow[datafl.ito] <= compto):
            return datarow[datafl.ito] - datarow[datafl.ifrom]
        # Case 2 - partially bottom of interval
        elif ((datarow[datafl.ifrom] >= compfrom) and (datarow[datafl.ifrom] < compto) and
              (datarow[datafl.ito] > compto)):

            return compto - datarow[datafl.ifrom]
        # Case 3 - partially in top of interval
        elif ((datarow[datafl.ifrom] < compfrom) and (datarow[datafl.ito] > compfrom) and
              (datarow[datafl.ito] <= compto)):

            return datarow[datafl.ito] - compfrom
        # Case 4 - composite bracketed by sample interval
        elif (datarow[datafl.ifrom] <= compfrom) and (datarow[datafl.ito] >= compto):
            return compto - compfrom
        # Case 5 - no length in interval
        else:
            return 0.0
    # Get length of each sample in composite
    try:
        dhdata['Length'] = dhdata.apply(length_in_interval, axis=1)
    except ValueError:  # No valid data in drill hole found - null assignment!
        return null
    # If no data was found, return null value
    if sum(dhdata['Length']) <= 0.0:
        return null
    # Upscale the data
    # - Continuous
    if vartype.startswith('con'):
        upscaled = np.dot(dhdata[var], dhdata['Length']) / sum(dhdata['Length'])
    # - Categorical
    elif vartype.startswith('cat'):
        unique_cats = list(set(dhdata[dhdata['Length'] > 0.0][var]))
        unique_cat_wts = [sum(dhdata[dhdata[var] == cat]['Length']) for cat in unique_cats]
        upscaled = unique_cats[unique_cat_wts.index(max(unique_cat_wts))]
    # - Invalid variable type
    else:
        raise Exception('Invalid vartype', vartype)
    return upscaled


def serial_get_comps(comps, datafl, vartypes='continuous', null=None):
    '''Returns a pandas DataFrame with upscaled composites

    Parameters:
        datafl (DataFile): pygeostat DataFile with dh, ifrom, ito and at least 1 variable to upscale
        comps (DataFrame): pandas DataFrame with datafl.dh, datafl.ifrom, datafl.ito
        corresponding to the composite locations

    Keyword Args:
        vartypes (str OR dict): 'continuous', 'categorical' or a dictionary of variables like:
        {'Cu':'continuous','Facies','categorial'}
        null: the null value (values less than or equal to this will not be used)

    Returns:
        upscaled (DataFrame): Pandas DataFrame with values of upscaled variable

    .. codeauthor: pygeostat development team 2014-04-03'''
    if null is None:
        null = Parameters['data.null']
    # Get a list of drill holes for composites
    dhids = sorted(list(set(comps[datafl.dh])))
    # Determine the number of composites
    ncomps = len(comps[datafl.dh])
    # Create a dictionary of variable upscaling methods if not supplied
    if isinstance(vartypes, str):
        vartypedict = {}
        for var in datafl.data.columns:
            if var not in [datafl.dh, datafl.ifrom, datafl.ito]:
                vartypedict[var] = vartypes.lower()
    else:
        vartypedict = vartypes
    # Determine number of variables
    nvar = len(vartypedict)
    # Initialize a pandas DataFrame for upscaled values
    upscaled = pd.DataFrame(np.ones([ncomps, nvar]) * null, columns=list(vartypedict))
    upscaled[datafl.dh] = comps[datafl.dh]
    upscaled[datafl.ifrom] = comps[datafl.ifrom]
    upscaled[datafl.ito] = comps[datafl.ito]
    # Upscale and return
    for var, vartype in vartypedict.items():
        # Define the variable specific upscaling function to apply
        def varupscale(compseries):
            'Define the variable specific upscaling function to apply'
            return(upscale(compseries[datafl.dh], compseries[datafl.ifrom],
                           compseries[datafl.ito], datafl, var, vartype, null))
        # Apply this to the variable
        upscaled[var] = comps.apply(varupscale, axis=1)
    return upscaled


def fast_comps(comps, datafl, null=None):
    '''Returns a pandas DataFrame with upscaled composites
    and ASSUMES NO MISSING VALUES

    Parameters:
        datafl (DataFile): pygeostat DataFile with dh, ifrom, ito and at least 1 variable to upscale
        comps (DataFrame): pandas DataFrame with datafl.dh, datafl.ifrom, datafl.ito
        corresponding to the composite locations

    Keyword Args:
        null: the null value (values less than or equal to this will not be used)

    Returns:
        upscaled (DataFrame): Pandas DataFrame with value of upscaled variable

    .. codeauthor: pygeostat development team 2014-04-03'''

    if null is None:
        null = Parameters['data.null']
    # Get a list of drill holes for composites
    dhids = sorted(list(set(comps[datafl.dh])))
    # Determine the number of composites
    ncomps = len(comps[datafl.dh])
    # Create a dictionary of variable upscaling methods if not supplied
    vartypedict = {}
    for var in datafl.data.columns:
        if var not in [datafl.dh, datafl.ifrom, datafl.ito]:
            vartypedict[var] = 'continuous'
    # Determine number of variables
    nvar = len(vartypedict)
    # Initialize a pandas DataFrame for upscaled values
    upscaled = pd.DataFrame(np.ones([ncomps, nvar]) * null, columns=list(vartypedict))
    upscaled[datafl.dh] = comps[datafl.dh]
    upscaled[datafl.ifrom] = comps[datafl.ifrom]
    upscaled[datafl.ito] = comps[datafl.ito]
    # Now upscale
    for rowidx, row in comps.iterrows():
        compfrom = row[datafl.ifrom]
        compto = row[datafl.ito]
        # Get the variable, froms and tos which are in that drill hole
        dhdata = datafl.data[datafl.data[datafl.dh] == row[datafl.dh]]
        # Define local function to get length in interval

        def length_in_interval(datarow):
            'Define local function to get length in interval'
            # Is the variable value valid?
            if datarow[var] <= null:
                return 0.0

            # Case 1 - entire sample interval in composite
            if (datarow[datafl.ifrom] >= compfrom) and (datarow[datafl.ito] <= compto):
                return datarow[datafl.ito] - datarow[datafl.ifrom]

            # Case 2 - partially bottom of interval
            elif ((datarow[datafl.ifrom] >= compfrom) and (datarow[datafl.ifrom] < compto) and
                  (datarow[datafl.ito] > compto)):

                return compto - datarow[datafl.ifrom]

            # Case 3 - partially in top of interval
            elif ((datarow[datafl.ifrom] < compfrom) and (datarow[datafl.ito] > compfrom) and
                  (datarow[datafl.ito] <= compto)):

                return datarow[datafl.ito] - compfrom

            # Case 4 - composite bracketed by sample interval
            elif (datarow[datafl.ifrom] <= compfrom) and (datarow[datafl.ito] >= compto):
                return compto - compfrom

            # Case 5 - no length in interval
            else:
                return 0.0
        # Get length of each sample in composite
        dhdata['Length'] = dhdata.apply(length_in_interval, axis=1)
        # If no data was found, return null value
        if sum(dhdata['Length']) <= 0.0:
            for var, vartype in vartypedict.items():
                upscaled.loc[rowidx, var] = null
        # Upscale the data
        for var, vartype in vartypedict.items():
            upscaled.loc[rowidx, var] = np.dot(dhdata[var],
                                               dhdata['Length']) / sum(dhdata['Length'])
    return upscaled


def write_vtp(datafl, vartypes, drillholes, outflname, complength=5.0, null=None):
    '''Generates a VTP file compatible with ParaView.

    Parameters:
        datafl (DataFile): pygeostat DataFile with dh, ifrom, ito and at least 1 variable to upscale
        vartypes (dict): dictionary of variable types like {'Grade':'continuous', 'Category':'categorical'}
        drillholes (dict): dictionary of drillhole IDs to pygeostat DrillHole objects
                           likely obtained by running:
                           drillholes = gs.set_desurvey(collarfl, surveyfl, 'Depth', 'Azimuth', 'Inclination')
        outflname (str): output VTP file to generate

    Keyword Args:
        null: the null value (values equal to this will not be used)
        complength: this is the composite length which everything is regularized to which helps
                    prevent visual artefacts and long triangles with a tube filter

    .. codeauthor: pygeostat development team 2017'''
    from ..datautils import fastcomps as gsfastcomps
    if null is None:
        null = Parameters['data.null']

    # "Composite" it down to reduce visual artefacts
    tmpfl = copy.deepcopy(datafl)
    comps = set_comps(datafl, complength)
    tmpfl.data = gsfastcomps.get_comps(comps, datafl, vartypes=vartypes, null=null)

    # Slice out any assays exceeding maximum hole depth
    tmpdfs = []
    for dhid, drillhole in drillholes.iteritems():
        tmpdfs.append(tmpfl.data[(tmpfl.data[tmpfl.dh] == dhid) &
                                 (tmpfl.data[tmpfl.ifrom] <= max(drillhole.data['Along'])) &
                                 (tmpfl.data[tmpfl.ito] <= max(drillhole.data['Along']))])
    tmpfl.data = pd.concat(tmpdfs, ignore_index=False)

    # Desurvey in all permutations...
    # Mid desurveying
    get_desurvey(tmpfl, drillholes, inplace=True, x='MIDX', y='MIDY', z='MIDZ')
    # Top desurveying
    tmp_name = tmpfl.ito
    tmpfl.ito = None
    get_desurvey(tmpfl, drillholes, inplace=True, x='TOPX', y='TOPY', z='TOPZ')
    # Bottom desurveying
    tmpfl.ito = tmp_name
    tmp_name = tmpfl.ifrom
    tmpfl.ifrom = None
    get_desurvey(tmpfl, drillholes, inplace=True, x='BOTTOMX', y='BOTTOMY', z='BOTTOMZ')
    tmpfl.ifrom = tmp_name

    # Now generate the VTK file
    root = ET.Element('VTKFile', type="PolyData", version="1.0", byte_order="LittleEndian", header_type="UInt64")
    pdata = ET.SubElement(root, 'PolyData')

    # Each drillhole is a separate "Piece" of "PolyData" in the VTK File
    for dhid, drillhole in drillholes.iteritems():
        dhdf = tmpfl.data[tmpfl.data[tmpfl.dh] == dhid]
        nverts = len(dhdf) + 2 # Add 2 since we have the top and bottom extra vertices

        # Header
        piece = ET.SubElement(pdata, 'Piece', NumberOfPoints="{}".format(nverts), NumberOfVerts="0",
                              NumberOfLines="1", NumberOfStrips="0", NumberOfPolys="0")

        # Vertex coordinates
        points = ET.SubElement(piece, 'Points')
        darray = ET.SubElement(points, 'DataArray', type="Float32", Name="Points",
                               NumberOfComponents="3", format="ascii")

        # Top and bottom vertices are added to make the drill hole match the true length
        ptarray = ' '.join(map(str, dhdf[['TOPX', 'TOPY', 'TOPZ']].values[0]))
        ptarray = ptarray + ' ' + ' '.join([' '.join(map(str, row[['MIDX', 'MIDY', 'MIDZ']].values)) for \
                                      rowidx, row in dhdf.iterrows()])
        ptarray = ptarray + ' ' + ' '.join(map(str, dhdf[['BOTTOMX', 'BOTTOMY', 'BOTTOMZ']].values[-1]))
        darray.text = ptarray

        # Variable values
        pointdata = ET.SubElement(piece, 'PointData', Scalars="DATA")
        for varname, vartype in vartypes.iteritems():
            if vartype.lower().startswith('cat'):
                vtktype = 'Int32'
                ptarray = str(int(dhdf[varname].values[0])) + ' ' + \
                          ' ' .join(map(str, map(int, dhdf[varname].values))) + ' ' + \
                          str(int(dhdf[varname].values[-1]))
            else:
                vtktype = 'Float32'
                ptarray = str(dhdf[varname].values[0]) + ' ' + \
                          ' ' .join(map(str, dhdf[varname].values)) + ' ' + \
                          str(dhdf[varname].values[-1])

            darray = ET.SubElement(pointdata, 'DataArray', type=vtktype, Name=varname, format="ascii")
            darray.text = ptarray

        # Line values
        lines = ET.SubElement(piece, 'Lines')
        darray = ET.SubElement(lines, 'DataArray', type="Int64", Name="connectivity", format="ascii")
        darray.text = ' '.join(map(str, range(0,nverts)))
        darray = ET.SubElement(lines, 'DataArray', type="Int64", Name="offsets", format="ascii")
        darray.text = '{}'.format(nverts)

    tree = ET.ElementTree(root)
    tree.write(outflname)
