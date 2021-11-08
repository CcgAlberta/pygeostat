# !/usr/bin/env python
#  -*- coding: utf-8 -*-
# public
""" Tools to build F2PY extensions for windows using Intel and Gfortran compilers,
visual studio 2012, 2015 and gcc C compilers """
from __future__ import absolute_import, division, print_function
from subprocess import call, Popen, PIPE, STDOUT
import os
import shutil
import sys
import glob
import copy
import numpy as np

#   Contributors
#   ------------
#   Version 1.00
#   Author  - Ryan Martin                    DATE: 2016-03-20
#
#   Version 1.10
#   Author  - Ryan Martin                    DATE: 2016-10-11
#           - updates for limiting the list of functions that
#               are wrapped from the Fortran source file
#

cdir = os.path.dirname(os.path.realpath(__file__))
pygeostatfolder = cdir[0:cdir.rfind('pygeostat')].replace("\\", '/')
rmfolder = pygeostatfolder + 'pygeostat/fortran/resource/rm/'


class FortranBuild:
    """
    This module can be used to generate F2PY modules using either gnu or intel compilers,
    tested for python 2.7, 3.4 and 3.5.

    This code depends on a set of paths defining compilers for the current machine. This file is
    located in the ``pygeostat/fortran`` folder, called ``userpaths.py``. If this file does not
    exist when either this script is run from the command line or is called from python, it will
    get created by making a copy of the ``defaultpaths.py`` file, found in the same folder.

    This script is setup such that as long as mingw was installed using ">conda install mingw" to
    the current python distribution (as defined in ``userpaths.py``), all required paths and
    libraries should already be defined to build modules using GNU compilers.

    If ``compiler='intel'`` is specified, make sure to define all paths and libraries in the file
    ``pygeostat/fortran/userpaths.py``

    If no sourcedict or libdict dictionaries are passed when the instance of this object is
    created, the object will collect lists of the default pygeostat sourcecodes and libraries
    defined in the ``sourcedefs.py`` file, and set the current source and out directory to
    'pygeostat' as default. With this in mind, if you do not overwrite these parameters for
    some other set of codes, your compiled modules will end up in the ``pygeostat/fortran`` folder!

    As always, take care to ensure sourcecode is error free and/or gnu compatible before
    compiling. If the \*.pyd is not generated AND there are no errors thrown by this script, there
    was most likely a fortran-compilation error in your source code, which you are able to see in
    the output on the terminal (i.e. the cmd window associated with the ipython notebook, etc)

    Parameters:
        sourcedict (dict): Dictionary of sourcecodes to build, format is {'key':['source.f90']}
            resulting in key.pyd from source.f90
        libdict (dict): Dictionary of included libraries for each matching key target in
            sourcedict. For example to include a LAPACK library on the linking step, define a
            matching key and the path to the library here.
        functiondict (dict): Dictionary function names with the corresponding list
            of functions that should get wrapped for python in the compiled module. format is
            {'key':['function_name']}. The data contained in the dictionary must be a list of
            strings or a string. function_name is case sensitive and must match the function name
            in the fortran source code *exactly*. Defaults is to wrap all functions for python.
        sourcepath (string): default = 'pygeostat', but your source code may be in the current
            folder (i.e. = './')
        outpath (string): default = 'pygeostat', but you may want output to the current folder
            (i.e. = './')
        compiler (string): supports either 'gnu' or 'intel' by making sure paths in
            pygeostat/fortran/libdefs.py are properly defined
        VisualStudioVersion (string): default = None, in which case it will first search for the
            VS12 folder if it can't find that it will search for the VS15 folder. If you pass it any
            of the following options it will set `self.vsversion` and `self.cpu` accordingly.
            Available optiosn are `vs12_64`, and `vs15_64`
        mode (string): either `Release` or `Debug`. The `Debug` flag does not allow full debugging,
            but the terminal writeout will be more verbose about where things went wrong.
        usegcc (bool): `True` is only valid for the intel compiler with VS 2012. VS2015
            changed various basic libraries that requires the C-code to be compiled by 'cl'. In
            general if compiling with Intel Fortran (i.e. `compiler` == `intel`), this parameter
            should remain unchanged.

    Note:
        Building with ``intel`` and ``py36`` requires Visual Studio version 2015 or higher with the
        **compiler tools** also installed. Check the Visual Studio directory, in VC/ folder a batch
        file called ``vcvarsall.bat`` must be present for this to work!

    Note:
        In order to build 'all' modules, the library ``lapack_solve.a`` (for GNU) or
        ``lapack_solve_intel.lib`` (for intel) must be found in the ``pygeostat/fortran/resource``
        folder. This script will check for existence of this file prior to building 'all' or any
        module depending on this library, and throw errors if this required library is not found.
        To build this library for GNU, navigate to ``pygeostat/fortran/resource`` and use:

        > compile_lapack_solve.bat

        To build this library for Intel, specify the path to 'ifort.exe' in
        ``compile_lapack_solve_intel.py``, then call:

        > python compile_lapack_solve_intel.py

    Usage From Command Line:
        This function can be called from a command prompt opened in the ``pygeostat/fortran``
        folder:

        >>> python compile.py <args>

        To start, you can simply use:

        >>> python compile.py

        which will generate a summary of the available options, most importantly the valid arguments
        are:

        >>> python compile.py <compiler> <modname> <sourcefiles>

        * where valid compilers are: ``"intel" or "gnu"``

        * ``"modname"`` is the name of the resulting \*.pyd file, can be 'all' or one of the source
          keys defined in the userpaths.py file

        * all source files are listed in dependency order, last depends on first, the f2py wrapper
          **must** be last

        This initial call will also generate a ``userpaths.py`` file in the same directory, which
        contains the a set of default paths pointing to the intel 2016 compiler and the mingw gnu
        compilers. If mingw was installed with >>> conda install mingw, the paths should be ready
        to go. If using intel compilers, make sure to specify paths as indicated in the respective
        section.

        Command line usage of this tool proceeds as follows, to build all pygeostat modules with
        the 'gnu' compiler, use:

        >>> python compile.py gnu all

        To clean all temporary build files and compiled \*.pyd modules, use:

        >>> python compile.py clean

        To clean a single pyd from the build directory, for example covasubs.pyd, use:

        >>> python compile.py clean covasubs

    Usage From a Python Script:
        You can use this function to speed up co-development of fortran and python. By editing some
        Fortran source code, you can clean, recompile and rerun your testing script in a few easy
        lines as is demonstrated in the following. Here, assume that ``example.f90`` contains some
        source code that we are testing and developing.

        Setup the imports and create the required source dictionary containing the key and the .f90
        file definition (which is assumed to be in the same folder as this script is running from):

        >>> import pygeostat as gs
        >>> sourcedict = {'example': ['example.f90']}

        To create ``example.pyd`` in the current folder, **after checking the compiler path
        definitions in** ``userpaths.py``, initialize the class:

        >>> fbuild = gs.FortranBuild(sourcedict=sourcedict,sourcepath='./',outpath='./',
        ...                          compiler='gnu')

        Since nothing is built yet, simply call:

        >>> fbuild.build('all')

        which should generate the example.pyd in the current folder. This module can now be
        imported:

        >>> import example
        >>> print(example.__doc__)

        Now if the .f90 source code is edited, it is possible to call:

        >>> fbuild.clean()
        >>> fbuild.build('all')

        To ensure a clean rebuild of the source code.

        To specify that only the function :code:`test_function(test)` found in ``example.f90``
        should be wrapped for python, use:

        >>> import pygeostat as gs
        >>> sourcedict = {'example': ['example.f90']}
        >>> functiondict = {'example': ['test_function']}
        >>> fbuild = gs.FortranBuild(sourcedict=sourcedict,sourcepath='./',outpath='./',
        ...                          compiler='gnu', functiondict=functiondict)

        This will ensure that test_function is the only fortran function that can be called from
        python. This can be useful if the module contains many subroutines with `illegal` f2py
        coding styles, but it is only necessary to call a single subroutine/function from python.

    Note:
        When calling this function from python in an interactive ipython-notebook, after rebuilding
        the fortran source code, the new \*.pyd cannot be loaded until the kernel is restarted.
        Something about python and dll's prevent reloading, even though the \*.pyd is erased and
        recompiled, and an import is run on the new \*.pyd

    .. codeauthor:: Ryan Martin - 2016-03-20

    """

    def __init__(self, sourcedict=None, libdict=None, functiondict=None, sourcepath='pygeostat',
                 outpath='pygeostat', compiler='gnu', exfortflags=None, VisualStudioVersion=None,
                 mode='Release', usegcc=False):
        # Check and ensure a userpaths file is found in the pygeostat folder, and then work from
        # defs in that file
        if not os.path.isfile(pygeostatfolder + 'pygeostat/fortran/userpaths.py'):
            shutil.copy(pygeostatfolder + 'pygeostat/fortran/defaultpaths.py',
                        pygeostatfolder + 'pygeostat/fortran/userpaths.py')
        try:
            from . import userpaths as dd
        except:
            import userpaths as dd

        # set the compiler either 'intel' or 'gnu'
        compiler = compiler.lower()
        if compiler != 'intel' and compiler != 'gnu':
            raise Exception('ERROR: incompatible compiler specified, valid options are '
                            '"gnu" or "intel" \nwith the properly specified paths ')
        self.compiler = compiler.lower()

        # Check to see if a visual studio version is passed and cpu architecture
        if VisualStudioVersion is None:
            self.vsversion = 'auto'
            self.cpu = 'auto'
        elif VisualStudioVersion.lower() == 'vs12_64':
            self.vsversion = 12
            self.cpu = 64
        elif VisualStudioVersion.lower() == 'vs15_64':
            self.vsversion = 15
            self.cpu = 64
        else:
            raise Exception('ERROR: incompatible VisualStudioVersion specified, valid options '
                            'are `vs12_64`, `vs15_64`')  # `vs12_32`, `vs15_32`,

        # Ensure that if intel + gcc are set, only vs2012 can be used
        if self.compiler == 'intel' and self.vsversion == 15 and usegcc:
            updatestatus('NOTE: compile.py can only work with the following combinations: \n'
                         'intel fortran + vs2012 + gcc\n'
                         'intel fortran + vs2012 + cl\n'
                         'intel fortran + vs2015 + cl\n'
                         'gnu fortran + gcc')
        elif self.compiler.lower() == 'intel' and usegcc:
            self.vsversion = 12
            self.usecl = False
        elif self.compiler.lower() == 'intel' and not usegcc:
            self.usecl = True
        elif self.compiler.lower() == 'gnu':
            self.usegcc = True
            self.usecl = False

        # ------ Source Dictionary
        if sourcedict is None:
            self.sourcecode = get_dictionary('sourcecode', compiler)
        else:
            self.sourcecode = sourcedict

        # ------ Library Dictionary
        if libdict is None:
            self.libdict = get_dictionary('libraries', compiler)
        else:
            self.libdict = libdict

        # ------- Path Dictionary
        # get the rest of the paths from userpaths ONLY! no other options.
        self.pathdict = get_dictionary('compilerpaths', compiler)

        # ------- Path Dictionary for Intel
        if self.compiler.lower() == 'intel':
            self.pathdict['vccomp_dir'] = self._autouserpaths('vccomp_dir')
            self.pathdict['fcomp_dir'] = self._autouserpaths('fcomp_dir')
            self.pathdict['link'] = self._autouserpaths('link')
            self.pathdict['ifort'] = self._autouserpaths('ifort')
            self.pathdict['cl'] = self._autouserpaths('cl')
            self.pathdict['h5libs'] = self._autouserpaths('h5libs_dir')
            self.pathdict['rm'] = rmfolder

        # ------- Source Code and Output Paths
        if sourcepath == 'pygeostat':
            self.pathdict['srcdir'] = dd.srcdir
        else:
            if sourcepath[-1] != '/':
                sourcepath = sourcepath + '/'
            self.pathdict['srcdir'] = sourcepath

        if outpath == 'pygeostat':
            self.outpath = pygeostatfolder + 'pygeostat/fortran/'
        else:
            self.outpath = outpath

        # -------- Extra Flags for the fortran compiler
        if exfortflags is None:
            self._exfortflags = ''
        else:
            self._exfortflags = exfortflags

        # -------- Debug Flags
        self._exf2pyflags = ''
        self._exlinkflags = ''
        if mode.lower() == 'debug':
            self._exfortflags += ' /debug /traceback /check:bounds /check:stack /warn:interfaces '
            self._exlinkflags += ' /DEBUG /pdb:link.pdb '
            self._exf2pyflags += ' --debug-capi '

        # -------- Flags to construct the list of `functiondict` functions:
        self._f2pyfunctiondict = {}
        if functiondict is not None:
            # build the key: ' only:function1,function2: ' string
            assert(isinstance(functiondict, dict) is True), "Function dict must be a dictionary!"
            for modulename in self.sourcecode.keys():
                if modulename in functiondict:
                    funcdata = functiondict[modulename]
                    funcstr = ' only: '
                    if isinstance(funcdata, list):
                        funcstr += ' '.join(functiondict[modulename])
                    elif isinstance(funcdata, str):
                        funcstr += funcdata
                    else:
                        raise ValueError('Must have functions as a string or a list of strings')
                    funcstr += ' : '
                    self._f2pyfunctiondict[modulename] = funcstr

        # The following dictionaries are used to check paths to things during
        #     compilation (or before?)
        if self.compiler.lower() == 'gnu':
            self.pathlist = [self.pathdict['srcdir'], self.pathdict['mingwdir'],
                             self.pathdict['pypath']]
            self.pathvarname = ['srcdir', 'mingwdir', 'pypath']
            self.filelist = [self.pathdict['dllpython'], self.pathdict['fort_comp'],
                             self.pathdict['cc']]
            self.filevarname = ['dllpython', 'fort_comp', 'cc']
            self.h5libs = {}
        elif self.compiler.lower() == 'intel':
            self.pathlist = [self.pathdict['srcdir'], self.pathdict['fcomp_dir'],
                             self.pathdict['vccomp_dir'], self.pathdict['pypath']]
            self.pathvarname = ['srcdir', 'fcomp_dir', 'vccomp_dir', 'pypath']
            self.filelist = [self.pathdict['libpython'], self.pathdict['ifort'],
                             self.pathdict['link'], self.pathdict['cl']]
            self.filevarname = ['libpython', 'ifort', 'link', 'cc']

        # Ensure the vcvarsall environment variables are correctly set if we are using cl.exe
        # basically from: http://stackoverflow.com/a/26168133
        if self.usecl:
            if self.cpu == 32:
                configstr = 'x86'
            else:
                configstr = 'amd64'
            with open('tempbat.bat', 'w') as p:
                p.write('call "' + self.pathdict['vccomp_dir'] + 'vcvarsall.bat" %s' % configstr)
                p.write('\nset > tempenvs.txt')
            # run the wrapper batch file, creates tempenvs.txt
            os.system('call tempbat.bat')
            # import the variables to the current ennvironment
            f = open('tempenvs.txt', 'r')
            lines = f.read().splitlines()
            f.close()
            # cleanup
            os.remove('tempenvs.txt')
            os.remove('tempbat.bat')
            self._env = copy.deepcopy(os.environ)
            for line in lines:
                pair = line.split('=', 1)
                self._env[pair[0]] = pair[1]
        else:
            # copy the environment to this object so that it can be passed to subprocess later!
            self._env = copy.deepcopy(os.environ)

    def _autouserpaths(self, path):
        """
        Utility function to try and automatically determine user system paths
        """
        try:
            from . import userpaths as dd
        except:
            import userpaths as dd

        if path == 'fcomp_dir':
            if dd.fcomp_dir == 'auto':
                fcomp_dir = None
                if os.path.isdir('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017'
                                 '/windows/'):
                    # ver17
                    fcomp_dir = ('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017'
                                 '/windows/')
                elif os.path.isdir('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_'
                                   '2016/windows/'):
                    # ver16
                    fcomp_dir = ('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016'
                                 '/windows/')
                elif os.path.isdir('C:/Program Files (x86)/Intel/Composer XE 2013 SP1/'):
                    # ver14
                    fcomp_dir = 'C:/Program Files (x86)/Intel/Composer XE 2013 SP1/'
                elif os.path.isdir('C:/Program Files (x86)/Intel/Composer XE 2013/'):
                    # ver14
                    fcomp_dir = 'C:/Program Files (x86)/Intel/Composer XE 2013/'
                if fcomp_dir is None:
                    raise ValueError("The path for `fcomp_dir` could not be determined"
                                     "automatically. Please set it manually.")
            else:
                fcomp_dir = dd.fcomp_dir.strip()
                if not fcomp_dir.endswith('/'):
                    fcomp_dir += '/'
            return fcomp_dir

        if path == 'vccomp_dir':
            if dd.vccomp_dir == 'auto':
                vccomp_dir = None
                temp_dir = 'C:/Program Files (x86)/Microsoft Visual Studio {ver}.0/VC/'.format
                # check first to see if a vs version was passed and set
                if self.vsversion == 12:
                    if os.path.isfile(temp_dir(ver=11) + 'vcvarsall.bat'):
                        vccomp_dir = temp_dir(ver=11)
                    else:
                        raise Exception("Visual Studio Version 12 passed but correct directory not "
                                        "Found.")
                elif self.vsversion == 15:
                    if os.path.isfile(temp_dir(ver=14) + 'vcvarsall.bat'):
                        vccomp_dir = temp_dir(ver=14)
                    else:
                        raise Exception("Visual Studio Version 15 passed but correct directory not "
                                        "Found.")
                # if a vs version wasn't passed then look for the directories vs12 is default
                elif os.path.isfile(temp_dir(ver=14) + 'vcvarsall.bat'):
                    vccomp_dir = temp_dir(ver=14)
                elif os.path.isfile(temp_dir(ver=12) + 'vcvarsall.bat'):
                    vccomp_dir = temp_dir(ver=12)
                elif os.path.isfile(temp_dir(ver=11) + 'vcvarsall.bat'):
                    vccomp_dir = temp_dir(ver=11)
                elif os.path.isfile(temp_dir(ver=9) + 'vcvarsall.bat'):
                    vccomp_dir = temp_dir(ver=9)
                if vccomp_dir is None:
                    raise ValueError("The path for `vccomp_dir` could not be determined"
                                     "automatically. Please set it manually.")
            else:
                vccomp_dir = dd.vccomp_dir.strip()
            if not vccomp_dir.endswith('/'):
                vccomp_dir += '/'
            return vccomp_dir

        if path == 'link':
            if dd.link == 'auto':
                try:
                    vccomp_dir = self.pathdict['vccomp_dir']
                except:
                    raise ValueError("The path for `vccomp_dir` needs to be set before trying to"
                                     " automatically determining the path for `link`")
                link = vccomp_dir + 'bin/amd64/link.exe'
                if not os.path.isfile(link):
                    raise ValueError("The path for `link` could not be determined automatically."
                                     "Please set it manually.")
            else:
                link = dd.link.strip()

            return link

        if path == 'ifort':
            if dd.ifort == 'auto':
                try:
                    fcomp_dir = self.pathdict['fcomp_dir']
                except:
                    raise ValueError("The path for `fcomp_dir` needs to be set before trying to"
                                     " automatically determining the path for `ifort`")
                ifort = fcomp_dir + 'bin/intel64/ifort.exe'
                if not os.path.isfile(ifort):
                    raise ValueError("The path for `ifort` could not be determined automatically."
                                     "Please set it manually.")
            else:
                ifort = dd.ifort.strip()

            return ifort

        if path == 'h5libs_dir':
            h5libs_dir = None
            if dd.h5libs_dir == 'auto':
                try:
                    link = self.pathdict['link']
                except:
                    h5libs_dir = None
                # Get the visual studio version number
                vsver = None
                h5libs_dir = None
                try:
                    from win32com.client import Dispatch
                    ver_parser = Dispatch('Scripting.FileSystemObject')
                    vsver = ver_parser.GetFileVersion(link).split('.')[0]
                except:
                    pass
                # Try and get the
                if vsver is not None:
                    h5libs_dir = './resource/HDF5/VS{ver}_x{bit}'.format
                    if os.path.isdir(h5libs_dir(ver=vsver, bit=64)):
                        h5libs_dir = h5libs_dir(ver=vsver, bit=64)
                    elif os.path.isdir(h5libs_dir(ver=vsver, bit=32)):
                        h5libs_dir = h5libs_dir(ver=vsver, bit=32)
                if not isinstance(h5libs_dir, str) or not (os.path.isdir(h5libs_dir + '/lib') and
                                                           os.path.isdir(h5libs_dir + '/static')):
                    h5libs_dir = None
            else:
                h5libs_dir = dd.h5libs_dir

            return h5libs_dir

        if path == 'cl':
            cl = self.pathdict['vccomp_dir'] + 'bin/amd64/cl.exe'
            if not os.path.isfile(cl):
                raise ValueError("The path for `cl` could not be determined automatically."
                                 "Please set it manually.")
            return cl

    def clean(self, opt='all'):
        """
        A function that can be used to clean tmp/ and .pyd files for the build

        Parameters:
            opt (string) : Can be 'all', 'tmp', 'pyd', or one of the keys in the source
                           dictionary in ``self.sourcedict``, to remove all files pertaining
                           to that module

        """
        def cleanpyd(obj, key=None):
            """
            Function to clean pyd's from the outfolder, by calling the rm.exe cygwin
            tool if the pyd is in use
            """
            if key is not None:
                iterlist = [key]
            else:
                iterlist = list(sorted(obj.sourcecode.keys()))
            for item in iterlist:
                if os.path.isfile(obj.outpath + item + '.pyd'):
                    try:
                        os.remove(obj.outpath + item + '.pyd')
                    except PermissionError:
                        updatestatus('Looks like the pyd is in use! attempting to use the gnu '
                                     '"rm.exe" tool!')
                        try:
                            if os.path.isfile(rmfolder + 'rm.exe'):
                                call(rmfolder + 'rm.exe ' + obj.outpath + item + '.pyd')
                                updatestatus('Success!')
                            else:
                                call('rm ' + obj.outpath + item + '.pyd')
                                updatestatus('Success!')
                        except FileNotFoundError:
                            raise Exception(obj.outpath + item + '.pyd is in use and the cygwin '
                                            'utility rm was not found in the path. Check your '
                                            'cygwin path variables or unload the .pyd by '
                                            'restarting the kernel!')

        def cleantmp(obj):
            """
            Function to clean the tmp build directories, just deletes it all !
            """
            if obj.pathdict['tdir'] != './' and \
               obj.pathdict['tdir'] != '' and obj.pathdict['tdir'] != ' ':
                if os.path.isdir(obj.outpath + obj.pathdict['tdir']):
                    shutil.rmtree(obj.outpath + obj.pathdict['tdir'])
            else:
                for item in list(sorted(obj.sourcecode.keys())):
                    if os.path.isdir(obj.outpath + obj.pathdict['tdir'] + item):
                        shutil.rmtree(obj.outpath + obj.pathdict['tdir'] + item)

        if opt.lower() == 'all':
            cleanpyd(self)
            cleantmp(self)
        elif opt.lower() == 'tmp':
            cleantmp(self)
        elif opt.lower() == 'pyd':
            cleanpyd(self)
        elif opt in self.sourcecode.keys():
            cleanpyd(self, opt)
            cleantmp(self)
        else:
            raise Exception('Invalid opt passed to clean(). "all","tmp","pyd" or a key from the '
                            'code dictionary are valid arguments')
        # check for the extra link.pdb and vc110.pdb
        cwd = os.getcwd() + '\\'
        vcfile = glob.glob(cwd + 'vc*.pdb')
        if len(vcfile) > 0:
            if os.path.isfile(cwd + 'link.pdb'):
                vcfile.append(cwd + 'link.pdb')
            for file in vcfile:
                try:
                    os.remove(file)
                except PermissionError:
                    try:
                        if os.path.isfile(rmfolder + 'rm.exe'):
                            call(rmfolder + 'rm.exe ' + file)
                        else:
                            call('rm ' + file)
                    except:
                        pass

    def checkitems(self, builditem=None):
        """
        Functions to check for existence of paths before running things, throws errors and attempts
        to let the user know where to look to fix those errors
        """
        # check for missing folders in the path definitions
        for i, item in enumerate(self.pathlist):
            if not os.path.isdir(item):
                raise Exception('NOTE: path ' + item + ' does not exist. \nCheck definitions in ' +
                                'userpaths.py for variable: ' + self.pathvarname[i])
        # check for missing files in the file definitions
        for i, item in enumerate(self.filelist):
            if not os.path.isfile(item):
                raise Exception('NOTE: path ' + item + ' does not exist. \nCheck definitions in ' +
                                'userpaths.py for variable: ' + self.pathvarname[i])
        if builditem is None or builditem.lower() == 'all':
            # check for missing source code files in the file definitions dictionary
            for i, key in enumerate(list(sorted(self.sourcecode.keys()))):
                checkfile(key, self.sourcecode[key], self.pathdict['srcdir'])
                # check for missing libraries and extra files in the library definitions dictionary
                if key in self.libdict:
                    checkfile(key, self.libdict[key], self.pathdict['srcdir'])
        else:
            checkfile(builditem, self.sourcecode[builditem], self.pathdict['srcdir'])
            if builditem in self.libdict:
                checkfile(builditem, self.libdict[builditem], self.pathdict['srcdir'])

    def build(self, modulenames, outpath=None, verbose='none'):
        """
        This is the main function called on this class to build the module once the object has been
        initialized with source code and compiler paths and whatnot

        Parameters:
            modulenames (string) : can be 'all', or must match a key in the self.sourcecode
                dictionary
            outpath (string) : Path for .pyd output, can use this as a last minute switch to change
                               the desired outpath for this module
            verbose (string) : 'all' prints all output, 'min' prints each step, 'none' prints the
                               module name only
        """
        try:
            from . import userpaths as dd
        except:
            import userpaths as dd

        # make sure things are defined correctly before continuing
        self.verbose = verbose.lower()
        build_lapack(self.compiler, verbose=self.verbose)
        self.checkitems(modulenames)

        if outpath is not None:
            if outpath == 'pygeostat':
                self.outpath = dd.srcdir
            else:
                self.outpath = outpath
        if not os.path.isdir(self.outpath + self.pathdict['tdir']):
            os.makedirs(self.outpath + self.pathdict['tdir'])
        print('Building fortran modules to ' + self.outpath)
        # Create the local build object using the internal parameter dictionaries
        localbuild = Build(self)
        if modulenames == 'all':
            for name in list(sorted(self.sourcecode.keys())):
                # for 'all' dont build hdf5_io just yet...
                if name != 'hdf5_io':
                    if self.compiler == 'intel':
                        localbuild.buildintelfortran(name)
                    else:
                        localbuild.buildgnufortran(name)
        elif isinstance(modulenames, list):
            for name in modulenames:
                if name not in self.sourcecode.keys():
                    raise Exception(name + ' is not defined in the sourcecode dictionary')
                if self.compiler == 'intel':
                    localbuild.buildintelfortran(name)
                else:
                    localbuild.buildgnufortran(name)
        else:
            if modulenames not in self.sourcecode.keys():
                raise Exception('Source code for module "' + modulenames +
                                '" are not defined in sourcecode dictionary!')
            if self.compiler == 'intel':
                localbuild.buildintelfortran(modulenames)
            else:
                localbuild.buildgnufortran(modulenames)

    def clearuserpaths(self, forreals=False):
        """
        If for some reason you would like to clear the local userpaths.py by resetting to default
        call this function by passing forreals=True
        """
        if forreals:
            # Overwrite the local paths with the defaults
            shutil.copy(pygeostatfolder + 'pygeostat/fortran/defaultpaths.py',
                        pygeostatfolder + 'pygeostat/fortran/userpaths.py')


class Build:
    def __init__(self, existingobject):
        # copy all the parameters from the calling object to this object
        # keep this object and its functions hidden from the python calling scope!
        self.__dict__.update(existingobject.__dict__)

    def build_call(self, buildstr):
        """
        Function to permit continuous output from the called building executable to the ipython
        notebook
        """
        p = Popen(buildstr, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True,
                  env=self._env)
        fulloutput = ''
        for line in p.stdout:
            fulloutput += line
            if self.verbose in ('all'):
                print(line, end='')
        if ' error ' in fulloutput or 'Error:' in fulloutput:
            if self.verbose != 'all':
                print(fulloutput)
            extramessage = self.check_error_msg(fulloutput)
            # update with the final message and kill the build process and script
            updatestatus('**BUILD FAILED: Check errors ^^ \n' + extramessage)
            p.kill()
            sys.exit()

    def check_error_msg(self, fulloutput):
        """
        Function to check the full output for known strings and plausible fixes to the error.
        Future: add items to `edict` where the key is a unique string contained in the offending
        output, and the data is the reccomended solution to resolve the problem
        """
        edict = {'multiply': ('NOTE: you probably (might?) need to clean the `tmp/` folder!'),
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

    def buildgnufortran(self, buildname):
        """
        Using the target name, build all source and link with paths and flags specific to the GNU
        compiler.
        """
        # Dependent Paths and Files
        if self.verbose in ('all', 'min'):
            updatestatus('Building ' + buildname)
        tmppath = self.outpath + self.pathdict['tdir'] + buildname
        if not tmppath.endswith('/'):
            tmppath += '/'

        if os.path.isfile(self.pathdict['pypath'] + 'Scripts/f2py.exe '):
            f2pyexe = self.pathdict['pypath'] + 'Scripts/f2py.exe '
        elif os.path.isfile(self.pathdict['pypath'] + 'Scripts/f2py.py '):
            f2pyexe = 'python ' + self.pathdict['pypath'] + 'Scripts/f2py.py '
        else:
            raise Exception('No f2py found, check paths and numpy installation!')

        fort_flags = '-c -O2 ' + self._exfortflags + ' '

        cc_flags = '-DMS_WIN64 -Ofast '
        cc_include = '-I' + self.pathdict['pypath'] + 'lib/site-packages/numpy/core/include -I' + \
                     self.pathdict['pypath'] + 'include -I' + self.pathdict['pypath'] + \
                     'PCBuild -I"' + self.pathdict['srcdir'] + '" -I' + self.pathdict['pypath'] + \
                     'lib/site-packages/numpy/f2py/src  ' #+ \
                    #  '-IC:/Anaconda3/MinGW/x86_64-w64-mingw32/include'

        # link paths need " " quotations around each entry
        link_paths = "" #'"' + self.pathdict['mingwdir'] + 'x86_64-w64-mingw32/lib/libgfortran.dll.a" '

        link_flags = '-shared -o "' + self.outpath + buildname + '.pyd" '

        # BUILD THE MODULE!
        # f2py
        if self.verbose in ('all', 'min'):
            updatestatus('Calling F2PY')
        f2pycall = (f2pyexe + '-m ' + buildname + ' --build-dir "' + tmppath + '" ' +
                    '"' + self.pathdict['srcdir'] + '"' + self.sourcecode[buildname][-1])
        # add the list of 'only' functions to the wrapper
        if buildname in self._f2pyfunctiondict:
            f2pycall += self._f2pyfunctiondict[buildname]
        if self.verbose in ('all'):
            print(f2pycall)
        self.build_call(f2pycall)

        # fortran building all separately because of the outflag differences
        if self.verbose in ('all', 'min'):
            updatestatus('Calling GNU Fortran')
        for item in self.sourcecode[buildname]:
            if '/' in item:
                itemname = item[item.rfind('/') + 1:item.rfind('.')]
            else:
                itemname = item[0:item.rfind('.')]
            fortcall = (self.pathdict['fort_comp'] + fort_flags + ' "' + self.pathdict['srcdir'] +
                        item + '" -o"' + tmppath + itemname + '.o"' + ' -J"' + tmppath + '" ')
            if self.verbose in ('all'):
                print(fortcall)
            self.build_call(fortcall)

        # check for wrapping code generated from f2py, build if it exists
        if os.path.isfile(tmppath + buildname + '-f2pywrappers.f'):
            fortcall = (self.pathdict['fort_comp'] + fort_flags + ' "' + tmppath + buildname +
                        '-f2pywrappers.f" -o"' + tmppath + buildname + '-f2pywrappers.o"')
            if self.verbose in ('all'):
                print(fortcall)
            self.build_call(fortcall)
        if os.path.isfile(tmppath + buildname + '-f2pywrappers2.f90'):
            fortcall = (self.pathdict['fort_comp'] + fort_flags + ' "' + tmppath + buildname +
                        '-f2pywrappers2.f90" -o"' + tmppath + buildname + '-f2pywrappers.o"')
            if self.verbose in ('all'):
                print(fortcall)
            self.build_call(fortcall)

        # c source
        if self.verbose in ('all', 'min'):
            updatestatus('Calling GCC first time')
        cccall1 = (self.pathdict['cc'] + cc_flags + cc_include + ' -c ' + self.pathdict['pypath'] +
                   'Lib/site-packages/numpy/f2py/src/fortranobject.c -o "' + tmppath + 'fobj.o" ')
        if self.verbose in ('all'):
            print(cccall1)
        self.build_call(cccall1)
        if self.verbose in ('all', 'min'):
            updatestatus('Calling GCC second time')
        cccall2 = (self.pathdict['cc'] + cc_flags + cc_include + ' -c "' + tmppath + buildname +
                   'module.c" -o "' + tmppath + buildname + 'module.o" ')
        if self.verbose in ('all'):
            print(cccall2)
        self.build_call(cccall2)

        extralibs = ''
        # make sure linker is pointed to any of the libraries in inc_dict
        if buildname in self.libdict:
            extralibs = (extralibs + ' ' +
                         ' '.join(str('"' + self.pathdict['srcdir'] + libitem + '"') for
                                  libitem in self.libdict[buildname]) + ' ')
        # final call to the linker
        if self.verbose in ('all', 'min'):
            updatestatus('Calling Linker')
        linkcall = ('"' + self.pathdict['fort_comp'] + '" ' + link_flags + '"' + tmppath + '*.o" ' +
                    extralibs + link_paths + self.pathdict['dllpython'])
        if self.verbose in ('all'):
            print(linkcall)
        self.build_call(linkcall)

    def buildintelfortran(self, buildname):
        """
        Using the target name, build all source and link with paths and flags specific to the Intel
        Compiler.
        """
        # Dependent Paths and Files
        if self.verbose in ('all', 'min'):
            updatestatus('Building ' + buildname)
        tmppath = self.outpath + self.pathdict['tdir'] + buildname
        if not tmppath.endswith('/'):
            tmppath += '/'

        if os.path.isfile(self.pathdict['pypath'] + 'Scripts/f2py.exe '):
            f2pyexe = '"' + self.pathdict['pypath'] + 'Scripts/f2py.exe" '
        elif os.path.isfile(self.pathdict['pypath'] + 'Scripts/f2py.py '):
            f2pyexe = 'python "' + self.pathdict['pypath'] + 'Scripts/f2py.py" '
        else:
            raise Exception('No f2py found, check paths and numpy installation!')
        f2pyexe += self._exf2pyflags

        # Get the right libraries for HDF5 mod
        if 'hdf5_io' in buildname:
            if self.pathdict['h5libs'] is None:
                raise ValueError("The path for `h5libs_dir` could not be determined automatically."
                                 " Please set it manually.")
            else:
                self.pathdict['h5libs'] = self.pathdict['h5libs'].strip()
                if not (os.path.isdir(self.pathdict['h5libs'] + '/lib') and
                        os.path.isdir(self.pathdict['h5libs'] + '/static')):
                    raise ValueError("The path `h5libs_dir` is set to does not exist or does not"
                                     " contain the folders '/lib' and '/static'.")
            libfls = glob.glob(self.pathdict['h5libs'] + '/lib/*.lib')
            libfls = [os.path.abspath(libfl).replace("\\", "/") for libfl in libfls]
            hdf5libs = ' '.join(libfls)
            hdf5mods = os.path.abspath(self.pathdict['h5libs'] + '/static').replace("\\", "/")
            hdf5mods = ' /I:"' + hdf5mods + '" '
            mdflag = ' /MD '  # everything is defaulting to /MT (static multithreaded )

            # make sure gcc is used for C-code, with the
            intelmanglingflags = ''
            gccmanglingflags = ' -DUPPERCASE_FORTRAN -DNO_APPEND_FORTRAN '
            self.usecl = False
        else:
            mdflag = ' /MT '  # Make static pyds
            hdf5libs = ''
            hdf5mods = ''
            if self.usecl:
                intelmanglingflags = ' /names:lowercase /assume:underscore '
                gccmanglingflags = ' '
            else:
                intelmanglingflags = ' '
                gccmanglingflags = ' -DUPPERCASE_FORTRAN -DNO_APPEND_FORTRAN '

        # setup 32-64bit stuff:
        # NOTE! this is untested since it seems to depend on a 32bit python library
        if self.cpu == 32:
            self.pathdict['cl'] = self.pathdict['cl'].replace('amd64/', '')
            self.pathdict['link'] = self.pathdict['link'].replace('amd64/', '')
            self.pathdict['ifort'] = self.pathdict['ifort'].replace('intel64', 'ia32')
            cfg = 'x86'
        else:
            cfg = 'x64'

        # Intel Fortran Flags
        # in this case, /names:lowercase and /assume:underscore ensure that the ifort compiled
        # name manging matches both gcc.exe and cl.exe
        ifort_flags = ' ' + mdflag + '/c /assume:buffered_io /O3 ' + \
                      intelmanglingflags + \
                      ' /I"' + self.pathdict['fcomp_dir'] + 'mkl/include/fftw" ' + \
                      self._exfortflags + hdf5mods + ' /object:"' + tmppath + '" /module:"' + \
                      tmppath + '"'

        # -DUPPERCASE_FORTRAN -DNO_APPEND_FORTRAN  gcc mangling flags that make the
        # gcc name mangling compatible with default intel name mangling... currently these are
        # only enabled when compiling hdf5_io with gcc
        cc_flags = ' -DMS_WIN64 -Ofast ' + gccmanglingflags
        cc_include = '-I"' + np.get_include() + '" -I"' + \
                     self.pathdict['pypath'] + 'include" -I"' + self.pathdict['pypath'] + \
                     'PCBuild" -I"' + self.pathdict['srcdir'] + '" -I"' + self.pathdict['pypath'] + \
                     'lib/site-packages/numpy/f2py/src"'

        # cl flags, for the VC compiler, -O2 optimization with the rest of the include
        cl_flags = ' -O2 '
        cl_include = '-I"' + np.get_include() + '" -I"' + \
                     self.pathdict['pypath'] + 'include" -I"' + self.pathdict['pypath'] + \
                     'PCBuild" -I"' + self.pathdict['srcdir'] + '" -I"' + self.pathdict['pypath'] + \
                     'lib/site-packages/numpy/f2py/src"'

        # In VS2015, various 'legacy' functions which are now contained in this static library
        if self.vsversion == 15 and self.usecl:
            vc14_exlib = '"' + self.pathdict['vccomp_dir'] + \
                         'lib/amd64/legacy_stdio_definitions.lib"'
        else:
            vc14_exlib = ''

        # a fun compiler default dir change for libs... seeing errors like cannot find ifconsol.lib
        iflibdirs = ('compiler/lib/intel64', 'compiler/lib/intel64_win', 'mkl/lib/intel64',
                     'mkl/lib/intel64_win')
        iflibpath = ''
        for ifdir in iflibdirs:
            if os.path.isdir(self.pathdict['fcomp_dir'] + ifdir):
                iflibpath += '/libpath:"%s%s" ' % (self.pathdict['fcomp_dir'], ifdir)
        # link paths need " " quotations around each entry
        link_paths = iflibpath + \
                     '/libpath:"' + self.pathdict['vccomp_dir'] + 'lib/amd64" ' + vc14_exlib + \
                     ' /libpath:"' + self.pathdict['pypath'] + 'lib/site-packages/numpy/f2py/src" '
        link_flags = '/DLL /SUBSYSTEM:console /STACK:1512000000,1512000000 ' + \
                     '/MACHINE:' + cfg + ' /MAP:"' + tmppath + buildname + \
                     '.map"  ' + hdf5libs + ' ' + self._exlinkflags
        link_out = '/out:"' + self.outpath + buildname + '.pyd"'

        # =============== Actually start building the codes =======================================
        # f2py
        if self.verbose in ('all', 'min'):
            updatestatus('Calling F2PY')
        f2pycall = f2pyexe + '-m ' + buildname + ' --build-dir --lower "' + tmppath + '" "' + \
                   self.pathdict['srcdir'] + self.sourcecode[buildname][-1] + '"'
        # add the list of 'only' functions to the wrapper
        if buildname in self._f2pyfunctiondict:
            f2pycall += self._f2pyfunctiondict[buildname]

        if self.verbose in ('all'):
            print(f2pycall)
        self.build_call(f2pycall)

        # fortran
        if self.verbose in ('all', 'min'):
            updatestatus('Calling Intel Fortran')
        fortcall = self.pathdict['ifort'] + ifort_flags + ' ' + \
                   ' '.join(str(' "' + self.pathdict['srcdir'] + item + '"')
                            for item in self.sourcecode[buildname])
        # check for wrapping code generated from f2py, build if it exists
        if os.path.isfile(tmppath + buildname + '-f2pywrappers.f') or \
           os.path.isfile(tmppath + buildname + '-f2pywrappers2.f90'):
            fortcall = fortcall + ' "' + tmppath + '*.f*" '
        # check to see if there are any libraries that need linking, see if they exist, add to call
        if self.verbose in ('all'):
            print(fortcall)
        self.build_call(fortcall)

        # C-source
        if self.verbose in ('all', 'min'):
            updatestatus('Building fortranobject.c')
        # Modify the calls to build the C-source if using CL.exe:
        if self.usecl and self.vsversion in (15, 12, 'auto'):
            cccall1 = self.pathdict['cl'] + mdflag + cl_flags + cl_include + ' -c "' + \
                      self.pathdict['pypath'] + \
                      'Lib/site-packages/numpy/f2py/src/fortranobject.c" /Fo"' + tmppath + 'fobj.obj" '
        else:
            cccall1 = self.pathdict['cc'] + cc_flags + cc_include + ' -c "' + \
                      self.pathdict['pypath'] + \
                      'Lib/site-packages/numpy/f2py/src/fortranobject.c" -o "' + tmppath + 'fobj.o" '
        if self.verbose in ('all'):
            print(cccall1)
        self.build_call(cccall1)
        if self.verbose in ('all', 'min'):
            updatestatus('Building %smodule.c' % buildname)

        if self.usecl and self.vsversion in (15, 12, 'auto'):
            cccall2 = self.pathdict['cl'] + mdflag + cl_flags + cl_include + ' -c "' + \
                      tmppath + buildname + 'module.c" -Fo"' + tmppath + buildname + 'module.obj" '
        else:
            cccall2 = self.pathdict['cc'] + cc_flags + cc_include + ' -c "' + tmppath + buildname + \
                      'module.c" -o "' + tmppath + buildname + 'module.o" '
        if self.verbose in ('all'):
            print(cccall2)
        self.build_call(cccall2)

        extralibs = ''
        # make sure linker is pointed to any of the libraries in inc_dict
        if buildname in self.libdict:
            extralibs = extralibs + ' '.join(str('"' + self.pathdict['srcdir'] + libitem + '"')
                                             for libitem in self.libdict[buildname]) + ' '

        # final call to the linker
        if self.verbose in ('all', 'min'):
            updatestatus('Calling Linker')
        linkcall = '"' + self.pathdict['link'] + '" ' + link_flags + extralibs + link_paths + \
                   '"' + tmppath + '*.obj" "' + tmppath + '*.o" "' + self.pathdict['libpython'] + \
                   '" ' + link_out
        if self.verbose in ('all'):
            print(linkcall)
        self.build_call(linkcall)

        # cleanup linker files:
        if os.path.isfile(self.outpath + buildname + '.exp'):
            try:
                os.remove(self.outpath + buildname + '.exp')
            except:
                pass
        if os.path.isfile(self.outpath + buildname + '.lib'):
            try:
                os.remove(self.outpath + buildname + '.lib')
            except:
                pass


def get_dictionary(dicttype, compiler):
    """
    Function to return a dictionary of variables obtained from the userpaths.py either in the
    current folder or from the pygeostat/fortran/userpaths.py file
    Paths are compiler dependent, as defined in userpaths.py

    Parameters:
        dicttype (string)   : type of dictionary, i.e. `sourcecode`, `libraries`, `compilerpaths`
        compiler (string)   : compiler to grab paths for. `intel` or `gnu`

    """
    def spacecontrol(inputdict, s_ns_dict):
        """
        Small function to ensure spaces or no spaces at the end of the string based on the key
        specified in the s_ns_dict.. pretty clunky way of doing this.

        Parameters:
            inputdict (dict)    : The input dictionary to mess with
            s_ns_dict (dict)    : a space/no_space after dictionary. 0 = no space, 1 = space

        """
        for key in list(sorted(inputdict.keys())):
            if s_ns_dict[key] == 1 and inputdict[key][-1] != ' ':
                inputdict[key] = inputdict[key] + ' '
            elif s_ns_dict[key] == 0 and inputdict[key][-1] == ' ':
                inputdict[key] = inputdict[key].strip()
        return inputdict

    # Start of get_dictionary Function
    try:
        from . import userpaths as dd
        from . import sourcedefs as dds
    except:
        import userpaths as dd
        import sourcedefs as dds
    if compiler.lower() == 'gnu':
        if dicttype.lower() == 'sourcecode':
            return dds.gnu_build_dict
        elif dicttype.lower() == 'libraries':
            return dds.gnu_inc_dict
        elif dicttype.lower() == 'compilerpaths':
            pathdict = {'tdir': dd.tdir, 'srcdir': dd.srcdir, 'pypath': dd.pypath, 'cc': dd.cc,
                        'mingwdir': dd.mingwdir, 'dllpython': dd.dllpython,
                        'fort_comp': dd.fort_comp}
            s_ns_dict = {'tdir': 0, 'srcdir': 0, 'pypath': 0, 'cc': 1,
                         'mingwdir': 0, 'dllpython': 1, 'fort_comp': 1}
            fin_dict = spacecontrol(pathdict, s_ns_dict)
            return fin_dict
    elif compiler.lower() == 'intel':
        if dicttype.lower() == 'sourcecode':
            return dds.intel_build_dict
        elif dicttype.lower() == 'libraries':
            return dds.intel_inc_dict
        elif dicttype.lower() == 'compilerpaths':
            pathdict = {'tdir': dd.tdir, 'srcdir': dd.srcdir, 'pypath': dd.pypath, 'cc': dd.cc,
                        'libpython': dd.libpython, 'sdkpath': dd.sdkpath}
            s_ns_dict = {'tdir': 0, 'srcdir': 0, 'pypath': 0, 'cc': 1, 'fcomp_dir': 0,
                         'vccomp_dir': 0, 'libpython': 1, 'sdkpath': 0, 'ifort': 1, 'link': 1}
            fin_dict = spacecontrol(pathdict, s_ns_dict)
            return fin_dict


def updatestatus(msg):
    print('\n******************************************************************************')
    print(msg)
    print('\n******************************************************************************')


def checkfile(buildname, flname, folder=''):
    # Not using gs.printerr to keep this independent of pygeostat
    if type(flname) is list:
        for item in flname:
            if not os.path.isfile(folder + item):
                updatestatus('WARNING "' + buildname + '" BUILD FAILED\n' + folder + item +
                             ' does not exist.\n Check source the source definitions and compiler'
                             ' paths files')
                sys.exit()
    else:
        if not os.path.isfile(folder + flname):
            updatestatus('WARNING "' + buildname + '" BUILD FAILED\n' + folder + item +
                         ' does not exist.\n Check source the source definitions and compiler'
                         ' paths files')
            sys.exit()


def replaceslash(flname):
    f = open(flname, 'r')
    paths = f.read()
    f.close()
    npaths = paths.replace("\\", "/")
    f = open(flname, 'w')
    f.write(npaths)
    f.close()


def build_lapack(compiler, force_compile=False, verbose='none'):
    """A small utility to build the lapack library that is required in some of the .pyd builds.

    Parameters:
        compiler (str): 'intel' or 'gnu'
        force_compile (bool): If `True` will force compile even if lapack library is already built
        verbose (bool) : 'all' prints all output, 'min' only prints a message if it starts compiling
            lapack, `none` prints only if error is found

    .. codeauthor:: Tyler Acorn March 08, 2017
    """
    lapack_folder = pygeostatfolder + 'pygeostat/fortran/resource/'
    print(lapack_folder)
    if compiler.lower() == 'intel':
        script_loc = lapack_folder + 'compile_lapack_solve_intel.py'
        if not os.path.isfile(lapack_folder + 'lapack_solve_intel.lib') or force_compile is True:
            if verbose.lower() in ['all', 'min']:
                print('Compiling Lapack for Intel')
            process = Popen('python ' + script_loc, cwd=lapack_folder,
                            stdout=PIPE, stdin=PIPE, stderr=PIPE)
            message = process.communicate()
            if verbose.lower() == 'all' or 'ERROR' in message:
                print(''.join([msg.decode('ascii') for msg in message]))
        # Check to make sure the library is there
        if not os.path.isfile(lapack_folder + 'lapack_solve_intel.lib'):
            updatestatus('Error: lapack library for INTEL was not build')
    elif compiler.lower() == 'gnu':
        script_loc = lapack_folder + 'compile_lapack_solve.bat'
        if not os.path.isfile(lapack_folder + 'lapack_solve.a') or force_compile is True:
            if verbose.lower() in ['all', 'min']:
                print('Compiling Lapack for GNU')
            process = Popen(script_loc, cwd=lapack_folder,
                            stdout=PIPE, stdin=PIPE, stderr=PIPE)
            message = process.communicate()
            if verbose.lower() == 'all' or 'ERROR' in message:
                print(''.join([msg.decode('ascii') for msg in message]))
        # Check to make sure the library is there
        if not os.path.isfile(lapack_folder + 'lapack_solve.a'):
            updatestatus('Error: lapack library for GNU was not build')
    else:
        raise Exception('choose either "intel" or "gnu" compiler')


# attempt to make this command line compatible also.. intended to be run from the pygeostat/fortran
# directory using:
# >>> python compile.py <args>
# details on usage in the docstring above!


if __name__ == '__main__':
    import argparse
    from sourcedefs import intel_build_dict as ibd
    from sourcedefs import gnu_build_dict as gbd
    sourcekeys = list(set((list(ibd.keys()) + list(gbd.keys()))))

    description_str = 'Process the build calls to this function. Valid modules are: %s' \
                       % ', '.join([str(key) for key in sourcekeys])
    parser = argparse.ArgumentParser(description=(description_str))

    parser.add_argument('modulename', type=str, nargs='?', default='list',
                        help=('the module to be built, use `list` to list'
                              ' permissable module names'))
    parser.add_argument('-compiler', type=str, default='intel',
                        help='Set the compiler to one of `intel` or `gnu`')
    parser.add_argument('-clean', type=str, default='None',
                        help='Clean the build directory, either `all` or one of `list`')
    parser.add_argument('-mode', type=str, default='release',
                        help=('Valid are default is `release`, or `debug`'))
    parser.add_argument('-fflags', type=str, default=None,
                        help='Extra fortran flags for compiling, e.g for intel: /check:bounds '
                             '/check:stack or  -flag for gnu compilers ')
    parser.add_argument('-vsvers', type=str, default=None,
                        help='Visual Studio Version, one of `vs12_64` or `vs15_64`')
    args = parser.parse_args()

    if not os.path.isfile('userpaths.py'):
        shutil.copy('defaultpaths.py', 'userpaths.py')

    if args.clean == 'all':
        for pyd in glob.glob('./*.pyd'):
            try:
                os.remove(pyd)
            except:
                try:
                    if os.path.isfile(rmfolder + 'rm.exe'):
                        call(rmfolder + 'rm.exe ' + pyd)
                    else:
                        call('rm ' + pyd)
                except:
                    raise Exception('The file: ' + pyd + ' is in use, and rm.exe '
                                    'was not found on the path. Close the file before cleaning '
                                    'or make sure rm.exe is callable from the system variable '
                                    'paths')
        if os.path.isdir('tmp'):
            try:
                shutil.rmtree('tmp')
            except:
                raise Exception('The folder tmp/ is in use, close command prompts or files '
                                'that might be open in this folder!')
    elif args.clean in sourcekeys:
        if os.path.isfile(args.clean + '.pyd'):
            try:
                os.remove(args.clean + '.pyd')
            except:
                try:
                    if os.path.isfile(rmfolder + 'rm.exe'):
                        call(rmfolder + 'rm.exe ' + args.clean + '.pyd')
                    else:
                        call('rm ' + args.clean + '.pyd')
                except:
                    raise Exception('The file: ' + args.clean + ' is in use, and rm.exe '
                                    'was not found on the path. Close the file before cleaning '
                                    'or make sure rm.exe is callable from the system variable '
                                    'paths')
        else:
            updatestatus(args.clean + '.pyd does not exist!')
        # clean the temp directory
        if os.path.isdir('tmp'):
            try:
                shutil.rmtree('tmp')
            except:
                raise Exception('The folder tmp/ is in use, close command prompts or files '
                                'that might be open in this folder!')
    elif args.modulename == 'all':
        modules = FortranBuild(compiler=args.compiler, mode=args.mode)
        if 'hdf5_io' in modules.sourcecode.keys():
            del modules.sourcecode['hdf5_io']
        modules.build('all', verbose='all')
    elif args.modulename in sourcekeys:
        buildall = FortranBuild(compiler=args.compiler, mode=args.mode)
        buildall.build(args.modulename, verbose='all')
    elif args.modulename == 'list':
        parser.print_help()
    else:
        raise Exception(args.modulename + ' is not in the default build dicitonary! Try one '
                        ' of: \n%s ' % '\n'.join([key for key in
                                                    sorted(sourcekeys)]))
