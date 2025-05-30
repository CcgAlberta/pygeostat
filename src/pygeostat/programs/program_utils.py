#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''programs.py: Contains methods for running GSLIB and CCG programs'''
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np

from ..utility.logging import log_progress, printerr
from .. pygeostat_parameters import Parameters


def runparallel(gslibprogram, kwargslist, nprocess=None, mute=False, progressbar=False):
    '''
    Run a set of gslib program calls in parallel

    Parameters:
        gslibprogram (Program): name of GSLIB/CCG program to run
        kwargslist (list of dictionaries): list of keyword arguments which will be used to call
            gslibprogram.run(kwarg)
        nprocess (int): number of threads to spawn.  Drawn from Parameters['config.nprocess'] if None.

    Examples:

        Setting up the calling parameters. This example is based off the example used in
        :func:`gs.Program() <pygeostat.programs.programs.Program>`

        >>> callpars = []
        >>> # For each variable we want to run in parallel, assemble a dictionary of parameters and
        ... # append to callpars
        >>> for variable in ['Bitumen','Fines','Chlorides']:
        >>>     # Establish the parameter file for this variable
        >>>     mypars = {'datafl':datafl.flname,
        ...               'varcol':datafl.gscol(variable),
        ...               'varname':variable,
        ...               'outfl':'histplt_'+variable+'.ps'}
        >>>    # Assemble the arguments for the GSLIB call and add the arguments to the list of calls
        >>>    callpars.append({'parstr':histpltpar.format(**mypars),
        ...                     'parfile':'histplt_'+variable+'.par',
        ...                     'testfilename':datafl.flname})

        Now run in parallel

        >>> histplt = gs.Program(program='histplt', parfile='histplt.par')
        >>> gs.runparallel(histplt, callpars)
        '''
    import multiprocessing as mp
    from os import remove

    progstr = gslibprogram.program
    if '/' in progstr:
        progstr = progstr[progstr.rfind('/') + 1:]
    if '.' in progstr:
        progstr = progstr[:progstr.rfind('.')]
    if not mute:
        print('Creating parallel processes')
    # Create a pool for processing, calling with None to use all CPUs
    if nprocess is None:
        nprocess = Parameters['config.nprocess']
    pool = mp.Pool(processes=nprocess)
    # Before executing, check to see if a parfile is present in each
    # kwarg. If not present then append a parfile that is unique for
    # each execution
    temppars = []
    for i, kwargs in enumerate(kwargslist):
        if 'parfile' not in kwargs.keys():
            temppar = progstr + "_" + temp_gslib_file("par")
            kwargslist[i]['parfile'] = temppar
            temppars.append(temppar)
    if progressbar:
        pbar = log_progress(range(len(kwargslist)))

        def pbar_update(*args, **kwargs):
            pbar.update()
    else:
        pbar_update = None
    # Now call in parallel for each kwargs
    for kwargs in kwargslist:
        # Add on the GSLIB call to execute asynchronously
        pool.apply_async(gslibprogram.run, (), kwargs, callback=pbar_update)
    # Close the list of processes for the pool
    if not mute:
        print('Pool assembled, asynchronously processing')
    pool.close()
    # "Join" the list of processes to execute and wait for completion
    pool.join()
    if not mute:
        print('Asynchronous execution finished.')
    if progressbar:
        pbar.close()
    # Remove the temporary parameter files if any were created
    if len(temppars) > 0:
        for temppar in temppars:
            remove(temppar)


class GennyWithAMemry:
    """
    A random seed generator object, callable, to replace gs.rseed, to generate more unique
    random seeds
    """

    def __init__(self, seed=None):
        self.seed = seed
        self._prng = np.random.RandomState(self.seed)
        self._rnd_list = []

    def _grnd(self, nseeds=100000):
        rndstart = self._prng.randint(0, 232135)
        dat = np.arange(rndstart, rndstart + nseeds)
        self._prng.shuffle(dat)
        return dat.tolist()

    def __call__(self, prng=None):
        '''Returns a ACORNI (GSLIB-suitable) random number seed'''
        if prng is not None:
            if isinstance(prng, int):
                self._prng = np.random.RandomState(prng)
            else:
                self._prng = prng
            self._rnd_list = []
        if len(self._rnd_list) == 0:
            self._rnd_list = self._grnd()
        return self._rnd_list.pop(0)


# this function replaces the rseed function to generate random seeds that are guaranteed unique
# for up to 100000 calls to the function. For more unique random seeds, use `rseed_list`
rseed = GennyWithAMemry()


def rseed_list(nseeds, seed=None):
    '''
    Returns a list of ACORNI (GSLIB-suitable) random number seeds. A initial seed can be passed
    ensureing the same list of seeds is returned to a script that is rerun.

    Parameters:
        nseeds (int): Number of seeds to return
        seed (int): Initialization seed

    Returns:
        seeds (list): List of random number seeds

    '''
    return GennyWithAMemry(seed)._grnd(nseeds)


def parallel_function(function, arglist=None, kwarglist=None, nprocess=None, returnvals=False,
                      progressbar=False):
    """
    Quickly parallelize a function with a set of arguments or keyword arguments. If the function
    returns something (as oppose to writes out values to files), set returnvals=True to get the
    dictionary that can be used to collect the results.

    Parameters:
        function (func): a callable function **DEFINED IN A .py FILE**. The function must be imported
            from a module since it has to be pickled to be parallelized. Defining the function in
            the jupyter notebook doesnt seem to work.
        arglist (list or tuples): a list of **tuple** arguments to pass to the function, see
            examples
        kwarglist (list): a list of keyword dictionaries to pass to the functions, i.e.
            [{'arg1': value, 'arg2': value}, {'arg1': value, 'arg2': value}, etc]
        nprocess (int): the number of parallel processes to run. Drawn from Parameters['config.nprocess']
            if None.
        returnvals (bool): if the function returns something, collect it in a dictionary, you can
            use the .get() method of the parallel result to collect the required data.

    Returns:
        res (dict): optionally return a dictionary of parallel processing results


    Usage:

    For a function that takes a single argument, setup arglist in this way:

        >>> arglist = []
        >>> arglist.append((arg,))

    OR:

        >>> arglist = [(arg,) for arg in range(nparallel)]

    Alternatively the argument tuple may take several arguments:

        >>> arglist = [(args, for, function), (args, for, function), (args, for, function)]

    Detailed usage:

        >>> arglist = []
        >>> for sr in series:
        ...     arglist.append(('keyout.out', rbfpath + 'keyout%s.out' % sr, griddef,
        ...                     griddefs[sr], [3]))
        >>> gs.parallel_function(rm.changegrid, arglist=arglist)

    OR:

        >>> kwarglist = []
        >>> for sr in series:
        ...     kwarglist.append({'infl': 'keyout.out',
        ...                       'outfl': rbfpath + 'keyout%s.out' % sr,
        ...                       'ingrid': griddef,
        ...                       'outgrid': griddefs[sr],
        ...                       'avmethods': [3]
        ...                       })
        >>> gs.parallel_function(rm.changegrid, kwarglist=kwarglist)
    """
    # some minor checks
    if arglist is None and kwarglist is None:
        raise ValueError('Must pass either lists of arglist or a list of kwarglist')
    if arglist is not None and kwarglist is not None:
        raise ValueError('Pass arglist OR kwarglist.. not both!')
    # start the parallel pool
    import multiprocessing as mp
    if nprocess is None:
        nprocess = Parameters['config.nprocess']
    if progressbar:
        pbar = log_progress(range(len(kwarglist if arglist is None else arglist)))

        def pbar_update(*args, **kwargs):
            pbar.update()
    else:
        pbar_update = None
    pool = mp.Pool(processes=nprocess)
    res = {}
    if arglist is not None:
        if not all([isinstance(item, tuple) for item in arglist]):
            printerr('all arguments in `arglist` must be of type `tuple`. If only 1 argument is '
                     'being passed to the parallel function, use: [(arg,)] (see docs). Attempting '
                     'to automagically convert the arguments to tuples ...')
        try:
            arglist = [item if isinstance(item, tuple) else tuple([item]) for item in arglist]
        except TypeError:
            printerr('items in `arglist` cannot be converted to type `tuple` automatically')
        # call the parallel function with the arguments
        for i, args in enumerate(arglist):
            res[i] = pool.apply_async(function, args, callback=pbar_update)
    else:
        # call the parallel function with the kwargs
        for i, kwargs in enumerate(kwarglist):
            res[i] = pool.apply_async(function, (), kwargs, callback=pbar_update)
    pool.close()
    pool.join()
    if progressbar:
        pbar.close()
    if returnvals:
        return res


def parstr_kwargs(parstr, fmt='pars'):
    """
    Print a formatted list of kwargs found in the parfile. Tested and working for the {} style
    string formatting, found in the gamsim_ave parfile found below this function (for exampel).
    This is mostly used for being lazy and not writing out the kwargs you just entered into the
    parfile.... can also be helpful if you define the parfile elsewhere and you cant remember what
    kwargs you setup !

    Parameters
    ----------
    parstr: str
        the par string with the {} formatted parameters
    fmt: str
        the output format, permissible arguments are `pars` or `dict`

    Examples
    --------
    Get the formatted parfile

        >>> parstr = '''          Parameters for GAMSIM_AVE
        >>>             *************************
        >>> START OF PARAMETERS:
        >>> {lithfl}        -file with lithology information
        >>> {lithcol}   {lithcode}                        -   lithology column (0=not used), code
        >>> {datafl}       -file with data
        >>> ....
        >>> '''

        >>> print_parfile_kwargs(parstr)
        ... lithfl=, lithcol=, lithcode=, datafl=, ...

    """
    from string import Formatter
    if fmt.lower() == 'dict':
        print(', '.join(["'%s': " % s[1] for s in Formatter().parse(parstr)][:-1]))
    else:
        print(', '.join(['%s=' % s[1] for s in Formatter().parse(parstr)][:-1]))


def dedent_parstr(indented_parstr):
    """
    Remove leading indents on each line from a parameter file. This is automatically called by
    :func:`Program.run() <pygeostat.programs.programs.Program.run>` so that parfiles may be tabbed to permit better
    structuring of python code.

    Examples
    --------
    An un-tabbed parstr:

        >>> parstr = '''          Parameters for GAMSIM_AVE
        >>>             *************************
        >>> START OF PARAMETERS:
        >>> {lithfl}        -file with lithology information
        >>> {lithcol}   {lithcode}                        -   lithology column (0=not used), code
        >>> {datafl}       -file with data
        >>> ....
        >>> '''

    A tabbed parstr:

        >>> parstr = '''          Parameters for GAMSIM_AVE
        >>>                 *************************
        >>>     START OF PARAMETERS:
        >>>     {lithfl}        -file with lithology information
        >>>     {lithcol}   {lithcode}                        -   lithology column (0=not used), code
        >>>     {datafl}       -file with data
        >>>     ....
        >>> '''
    """
    from textwrap import dedent
    dedented_parstr = ''
    for line in indented_parstr.splitlines():
        dedented_parstr += dedent(line) + '\n'
    return dedented_parstr


def temp_gslib_file(fltype='gslib'):
    """
    Generate unique temporary file names with different extensions as indicated by the ``fltype``
    parameter.

    Parameters
    ----------
        fltype: str
            one of "gslib", "par", "h5"
    """
    import tempfile
    import os
    cpath = os.path.relpath(os.getcwd())
    if fltype.lower() == 'gslib':
        return os.path.relpath(tempfile.NamedTemporaryFile(dir=cpath).name + '.out')
    elif fltype.lower() == 'par' or fltype.lower() == 'parfile':
        return os.path.relpath(tempfile.NamedTemporaryFile(dir=cpath).name + '.par')
    elif fltype.lower() == 'h5':
        return os.path.relpath(tempfile.NamedTemporaryFile(dir=cpath).name + '.h5')