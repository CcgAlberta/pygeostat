# !/usr/bin/env python
#  -*- coding: utf-8 -*-
# public
"""
Functions developed to permit compilation of ``pygeostat`` and arbitrary fortran modules. Using
``f2py`` and flexible Fortran and C compilers.

Note:
    ``compile.py`` is command line callable. Open a terminal in ``pygeostat/fortran/`` and build the
    pygeostat fortran extensions with ``python compile.py -compiler=gnu all``

.. codeauthor:: Ryan Martin - 24-04-2018
"""
from __future__ import absolute_import, division, print_function
from subprocess import call, check_output, Popen, PIPE, STDOUT
import os
import shutil
import sys
import glob
import copy
import numpy as np

# setup some default paths to consider for things relative to this file specifically
FILEDIR = os.path.dirname(os.path.realpath(__file__))
THISSCRIPT = os.path.abspath(__file__)
try:
    call("rm", stderr=PIPE, stdout=PIPE)
except FileNotFoundError:
    RMCALL = FILEDIR + '/resource/rm/rm.exe'
PYPATH = "/home/ryan/anaconda3/bin/"
PYVERSION = str(sys.version_info[0]) + str(sys.version_info[1])


def ensure_dir(f):
    """Make sure that directory(s) exists"""
    if isinstance(f, str):
        dirs = [f]
    else:
        dirs = copy.copy(f)
    for f in dirs:
        if not os.path.exists(f):
            os.makedirs(f)


def ensure_path(path):
    """
    Function ensures that all folders in a given path or list of paths are created if they
    do not exist
    """
    if isinstance(path, str):
        path = [path]
    path = [p.replace('\\', '/') for p in path]
    for p in path:
        incremental_path = ''
        dirlist = p.split('/')
        for d in dirlist:
            incremental_path += d + '/'
            if ':' not in d and d != '':
                ensure_dir(incremental_path)


def buildgnuf2py(mode, exflags, modulename, sourcefiles, includes=None, opt="O2", wraponly=None,
                 outdir="./", srcdir="./src/", tmpdir="./tmp/", env=None, verbose=True):
    """
    Compiles a single PYD using the gfortran compiler

    Parameters:
    -----------
        mode: str
            either `release` or `debug`
        exflags: str
            permissible gfortran compiler flags, space delimited,
        modulename: str
            the name of the resulting pyd, e.g. `covasubs` -> covasubs.pyd
        sourcefiles: list of str
            the source code files in order of dependency
        includes: list of str
            files that get included in the compile / link step
        opt: str
            The optimization level, possible to omit the `-`
        wraponly: list of str
            Names of the functions found in the final fortran file in ``sourcefiles`` that should
            be wrapped for python
        outdir, srcdir, tmpdir : str
            The path for output, sourcecode, tempfiles.
        env: dict
            A dictionary of environment variables to consider for the compiling process. Generated
            with ``env = assert_environ("gnu")`` or ``env = assert_environ("C:/mingw/bin")``
        verbose: bool
            Write compiler output to the terminal

    .. codeauthor:: Ryan Martin - 20-03-2018
    """
    tmppath = tmpdir + '{}'.format(modulename) + "/"
    ensure_path([outdir, tmpdir, tmppath])
    # generate the source files
    srcfiles = [os.path.join(srcdir, f) for f in sourcefiles]
    f2pycall = ['f2py', '-c', "--f90flags='-fPIC'", '-m', modulename]
    if mode == "debug":
        f2pycall.append("--debug-capi")
    if wraponly is not None:
        f2pycall.extend(["only:"] + wraponly + [":"])

    gfortbasecall = ["gfortran", "-fPIC", "-c", opt, '-J' + tmppath]

    if exflags is not None:
        exflagiter = exflags.split() if isinstance(exflags, str) else exflags
        for arg in exflagiter:
            if arg not in gfortbasecall:
                gfortbasecall.append(arg)

    objs2add = []

    def localbuild(fl):
        if '/' in fl:
            flname = fl[fl.rfind('/') + 1:fl.rfind('.')]
        else:
            flname = fl[0:fl.rfind('.')]
        localcall = copy.deepcopy(gfortbasecall)
        objfl = tmppath + flname + ".o"
        objs2add.append(objfl)
        localcall.extend([fl, "-o", objfl])
        _buildcall(localcall, env, verbose=verbose)

    for fl in srcfiles[:-1]:
        localbuild(fl)
    f2pycall.extend(objs2add + ["-I"+tmppath])
    if includes is not None: 
        if isinstance(includes, list):
            f2pycall.extend(includes)
        else:
            f2pycall.append(includes) ## this should be a string
    f2pycall.append(srcfiles[-1])
    _buildcall(f2pycall, verbose=verbose)


def _buildcall(buildstr, env=None, verbose=True, builddir=None):
    """
    Function to permit continuous output from the called building executable to the ipython
    notebook
    """
    if verbose:
        _updatestatus(buildstr)
    p = Popen(buildstr, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True,
              env=env, cwd=builddir)
    fulloutput = ''
    for line in p.stdout:
        fulloutput += line
        if verbose:
            print(line, end=" ")
    if ' error ' in fulloutput or 'Error:' in fulloutput.lower():
        if not verbose:
            print(fulloutput)
        extramessage = _checkerror(fulloutput)
        # update with the final message and kill the build process and script
        _updatestatus('**BUILD FAILED: Check errors ^^ \n' + extramessage)
        p.kill()
        sys.exit()


def _checkerror(fulloutput):
    """
    Function to check the full output for known strings and plausible fixes to the error.
    Future: add items to `edict` where the key is a unique string contained in the offending
    output, and the data is the reccomended solution to resolve the problem
    """
    edict = {'multiply': ('NOTE: you might(?) need to clean the `tmp/` folder!'),
             'already defined': ('NOTE: you probably (might?) need to clean the `tmp/` folder!'),
             'unresolved externals': ('NOTE: consider recompiling the linked libraries to'
                                      'have the correct name mangling for cl.exe:'
                                      'ifort: /names:lowercase /assume:underscore '),
             "KeyError: 'void'": ('There may be an issue with public/private function '
                                  'definitions or a missing variable definition in the last '
                                  'function listed above. For the first error consider using '
                                  'the parameter `functiondict` or checking to ensure all '
                                  'module functions are public... For the second error, check '
                                  'that all of the parameters in the subroutine are defined'),
             "No such file or directory": ('There may be a space in the path to one of the '
                                           'source code or library folders'),
             "LINK : fatal error LNK1104: cannot open file": ('The pyd is currently in use, '
                                                              'restart any kernels using it !')
             }
    # iterate through the keys in the error dictionary and see if the key is in the full output
    extramessage = ''
    for error_key in edict.keys():
        if error_key in fulloutput:
            extramessage = edict[error_key]
    return extramessage


def _updatestatus(msg, width=75):
    print('\n' + width * "*")
    print(msg)
    print('\n' + width * "*")


def build_pygeostat_fortran_linux(pygeostatmodules="all", compiler="gnu", mode="release", exflags="",
                                  opt="O2", verbose=False):
    """
    This function builds the f2py extensions with minimal interaction from the user. The goal is to
    make this function callable with ``gs.build_pygeostat_fortran()``, and have all requirements
    sorted out automatically resulting in a set of compiled pygeostat fortran extensions.

    Modules are compiled to ``pygeostat/fortran/``, regardless of where the function is called from.
    If no gnu compiling tools are found on the path (``gfortran -v``, ``gcc -v`` returns
    nothing), then ``conda install mingw -y`` is run to install MinGW, and this compiler is used.

    Parameters:
        pygeostatmodules (str): either ``"all"`` or one of the dictionary keys found in
            ``pygeostat/fortran/sourcedefs.py``
        compiler (str): either ``"intel"``, ``"gnu"`` or a path to a local folder containing gfortran
            compilers, e.g., ``C:/mingw/bin/``
        mode (str): "release" or "debug"
        exflags (str): compiler-specific permissable fortran compiler flags
        opt (str): optimization level, defaults to `O2`
        clean (bool): if `True` then the build directory is cleaned before compiling (recommended)
        verbose (bool): write all output to the terminal

    .. codeauthor:: Ryan Martin - 07-04-2018
    """
    try:
        from .sourcedefs import sources
    except:
        from sourcedefs import sources
    allpygeostat = list(sources.keys())
    assert pygeostatmodules in ["all"] + allpygeostat, (
        "{} is not in Pygeostat!".format(pygeostatmodules))
    # regardless of where the function is called the build is called here.
    compilerargs = ["-compiler={}".format(compiler), "-mode={}".format(mode), "-opt={}".format(opt)]
    if exflags != "":
        exflags = ['-exflags="{}"'.format(" ".join(exflags.split()))]
    else:
        exflags = []
    compilescript = "compile_linux.py"

    _updatestatus("Building ... ")
    # call to build the specified module
    buildstr = ["python", compilescript, pygeostatmodules] + compilerargs + exflags
    _buildcall(buildstr, None, verbose, FILEDIR)


if __name__ == '__main__':
    import argparse
    from subprocess import call
    from sourcedefs import sources, gnuincludes, wraponly

    allkeys = sources.keys()
    parser = argparse.ArgumentParser(
        description=('Process the build calls to this function. '
                     'Valid modules are: {}').format(', '.join([str(key) for key in allkeys])))
    parser.add_argument('modulename', type=str, nargs='?',
                        help=('the module to be built, use `list` to list'
                              ' permissable module names'))
    # parser.add_argument('-compiler', type=str, default='gnu',
                        # help='Set the compiler to one of `intel` or `gnu`')
    # parser.add_argument('-clean', type=str, default=None,
    #                     help='Clean the build directory, either `all` or one of `list`')
    parser.add_argument('-mode', type=str, default='release',
                        help=('Valid are default is `release`, or `debug`'))
    parser.add_argument('-fflags', nargs="+", type=str,
                        help='String of extra fortran flags for compiling, e.g for intel: `"/check:bounds '
                             '/check:stack"` or `"-flag"` for gnu compilers ')
    parser.add_argument('-tmpdir', type=str, default='./tmp/',
                        help=('Directory for temporary output'))
    parser.add_argument("-opt", type=str, default="-O2", help=("The optimization level"))
    parser.add_argument("-outdir", type=str, default="./",
                        help=("the directory where final pyds end up"))
    parser.add_argument('-forcelapack', action="store_true",
                        help="Force recompile of lapack library")
    # this is a list of strings corresponding to functions that should be wrapped with f2py
    parser.add_argument("-wraponly", nargs="+", type=str, default=None,
                        help="Space separated list of the functions that will be wrapped with f2py,"
                        "e.g. `-wraponly func1 func2` ")
    parser.add_argument("-relink", action="store_true",
                        help="Only perform the link step")
    args = parser.parse_args()
    origenv = os.environ.copy()
    # if args.clean is not None:
    #     extensioncleaner(args.clean, args.tmpdir, args.outdir)
    # elif args.modulename is None:
    #     parser.print_help()
    if args.modulename == "list":
        print('POSSIBLE KEYS: \nall  (use to build all modules)\n%s\n' %
              ' \n'.join([key for key in sorted(allkeys)]))
    else:
        if args.modulename != "all" and args.modulename not in allkeys:
            parser.print_help()
            _updatestatus("ERROR: \n" + args.modulename + ' is not in the default build dicitonary! '
                          'Try `all` or one of: \n%s\n ' % ', '.join([key for key in sorted(allkeys)]))
            sys.exit(0)
        if not os.path.isfile("./resource/lapack_solve.a"):
            call(["python", "compile_lapack_linux.py"], cwd="./resource/")
        else:
            print("Found LAPACK!")

        includes = gnuincludes
        ensure_path(args.tmpdir)
        buildkeys = allkeys if args.modulename == "all" else [args.modulename]
        # # actually call the build functions
        for buildname in buildkeys:
            sourcecode = sources[buildname]
            include = includes.get(buildname, None)
            towrap = wraponly.get(buildname, None) if args.wraponly is None else args.wraponly
            buildgnuf2py(args.mode, args.fflags, buildname, sourcecode, opt=args.opt,
                         wraponly=towrap, includes=include, tmpdir=args.tmpdir)
