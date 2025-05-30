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
except Exception:
    RMCALL = FILEDIR + '/resource/rm/rm.exe'
PYPATH = sys.executable[0:sys.executable.rfind('python.exe')].replace('\\', '/')
PYVERSION = str(sys.version_info[0]) + str(sys.version_info[1])
LIBPYTHON = PYPATH + 'libs/Python' + PYVERSION + '.lib'
DLLPYTHON = PYPATH + 'python' + PYVERSION + '.dll'


def quote(string):
    " surround the string quotes "
    return '"' + string + '"'


if os.path.isfile(PYPATH + 'scripts/f2py.py'):
    F2PY = ' '.join([sys.executable, quote(PYPATH + "scripts/f2py.py")])
else:
    F2PY = 'f2py '


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


def load_ifortenv(ifortdir=None, x86=False):
    """ Load environment variables """
    from subprocess import call
    if ifortdir is None:
        if os.path.isdir('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries'
                         '/windows/bin/'):
            # ver18?? or generic?
            ifortdir = ('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries'
                        '/windows/bin/')
        elif os.path.isdir('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017'
                           '/windows/bin/'):
            # ver17
            ifortdir = ('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017'
                        '/windows/bin/')
        elif os.path.isdir('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_'
                           '2016/windows/bin/'):
            # ver16
            ifortdir = ('C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016'
                        '/windows/bin/')
        elif os.path.isdir('C:/Program Files (x86)/Intel/Composer XE 2013 SP1/'):
            # ver14
            ifortdir = 'C:/Program Files (x86)/Intel/Composer XE 2013 SP1/'
        elif os.path.isdir('C:/Program Files (x86)/Intel/Composer XE 2013/'):
            # ver14
            ifortdir = 'C:/Program Files (x86)/Intel/Composer XE 2013/'
        if ifortdir is None:
            return -1
    if x86:
        configstr = 'ia32'
    else:
        configstr = 'intel64'
    with open('tempbat.bat', 'w') as p:
        p.write('call "{}/ifortvars.bat" {}'.format(ifortdir, configstr))
        p.write('\nset > tempenvs.txt')
    # run the wrapper batch file, creates tempenvs.txt
    call('tempbat.bat')
    # import the variables to the current ennvironment
    with open('tempenvs.txt', 'r') as f:
        lines = f.read().splitlines()
    # cleanup
    os.remove('tempenvs.txt')
    os.remove('tempbat.bat')
    _env = copy.deepcopy(os.environ)
    for line in lines:
        pair = line.split('=', 1)
        _env[pair[0]] = pair[1]
    return _env


# def buildgnudll(mode, exflags, modulename, sourcefiles, includes=None,
#                 outdir="../dlls/", srcdir="./src/", tmpdir="./tmp/", env=None):
#     """
#     Compiles a single DLL using the gfortran compiler

#     Parameters:
#     -----------
#         mode: str
#             either `release` or `debug`
#         exflags: str
#             permissible gfortran compiler flags, space delimited,
#         modulename: str
#             the name of the resulting dll, e.g. `covasubs` -> covasubs.dll
#         sourcefiles: list of str
#             the source code files in order of dependency
#         includes: list of str
#             files that get included in the compile / link step

#     .. codeauthor:: Ryan Martin - 20-03-2018
#     """
#     ensure_path([outdir, tmpdir])
#     buildcall = ["gfortran", "-shared"]
#     if isinstance(exflags, str):
#         buildcall.extend([f for f in exflags.split()])
#     elif isinstance(exflags, (list, tuple)):
#         buildcall.extend(exflags)
#     buildcall.extend([os.path.join(srcdir, f) for f in sourcefiles])
#     tmppath = tmpdir + '{}'.format(modulename)
#     ensure_path(tmppath)
#     buildcall.append('-J"' + tmppath + '"')
#     if includes is not None:
#         if isinstance(includes, str):
#             includes = [includes]
#         buildcall.extend([os.path.join(srcdir, f) for f in includes])
#     buildcall.append("-static")
#     buildcall.extend(["-o", "{}lib{}.dll".format(outdir, modulename)])
#     _updatestatus("Building {}.dll".format(modulename))
#     buildcall = " ".join(buildcall)
#     print(buildcall)
#     _buildcall(buildcall, env, verbose=verbose)


# def buildinteldll(mode, exflags, modulename, sourcefiles, includes=None,
#                   outdir="../dlls/", srcdir="./src/", tmpdir="./tmpd/", env=None):
#     """
#     Compiles a single DLL using the gfortran compiler

#     Parameters:
#     -----------
#         mode: str
#             either `release` or `debug`
#         exflags: str
#             permissible gfortran compiler flags, space delimited,
#         modulename: str
#             the name of the resulting dll, e.g. `covasubs` -> covasubs.dll
#         sourcefiles: list of str
#             the source code files in order of dependency
#         includes: list of str
#             files that get included in the compile / link step

#     .. codeauthor:: Ryan Martin - 20-03-2018
#     """
#     ensure_path([outdir, tmpdir])
#     buildcall = ["ifort", "/dll"]
#     if isinstance(exflags, str):
#         [buildcall.append(f) for f in exflags.split()]
#     elif isinstance(exflags, (list, tuple)):
#         buildcall.extend(exflags)
#     buildcall.extend(["/libs:static", "/threads", "/assume:buffered_io",
#                       "/names:lowercase", "/assume:underscore", "/O3", "/MT"])
#     buildcall.extend([os.path.join('"' + srcdir, f + '"') for f in sourcefiles])
#     tmppath = tmpdir + '{}'.format(modulename)
#     ensure_path(tmppath)
#     buildcall.extend(['/object:"' + tmppath + '/"', '/module:"' + tmppath + '/"'])
#     buildcall.append("/link")
#     if includes is not None:
#         if isinstance(includes, str):
#             includes = [includes]
#         buildcall.extend([os.path.join('"' + srcdir, I + '"') for I in includes])
#     buildcall.append('/out:"{}lib{}.dll"'.format(outdir, modulename))
#     _updatestatus("Building {}.dll".format(modulename))
#     buildcall = " ".join(buildcall)
#     print(buildcall)
#     _buildcall(buildcall, env, verbose=verbose)


def buildintelf2py(mode, exflags, modulename, sourcefiles, includes=None, opt="O2", wraponly=None,
                   outdir="./", srcdir="./src/", tmpdir="./tmp/", env=None, verbose=True):
    """
    Compiles a single PYD using the intel compiler

    Parameters:
    -----------
        mode: str
            either `release` or `debug`
        exflags: str, list
            List or str of space delimited intel compiler flags
        modulename: str
            the name of the resulting pyd, e.g. `covasubs` -> covasubs.pyd
        sourcefiles: list of str
            the source code files in order of dependency
        includes: list of str
            files that get included in the compile / link step
        opt: str
            The optimization level, possible to omit the `/`
        wraponly: list of str
            Names of the functions found in the final fortran file in ``sourcefiles`` that should
            be wrapped for python
        outdir, srcdir, tmpdir : str
            The path for output, sourcecode, tempfiles.
        env: dict
            A dictionary of environment variables to consider for the compiling process. Generated
            with ``env = assert_environ("intel")``
        verbose: bool
            Write compiler output to the terminal

    .. codeauthor:: Ryan Martin - 20-03-2018
    """
    tmppath = tmpdir + '{}'.format(modulename) + "/"
    ensure_path([outdir, tmpdir, tmppath])
    # generate the source files
    srcfiles = [quote(os.path.join(srcdir, f)) for f in sourcefiles]
    f2pycall = [F2PY, '-m', modulename, "--build-dir", quote(tmppath), "--lower"]
    if mode == "debug":
        f2pycall.append("--debug-capi")
    f2pycall.append(srcfiles[-1])
    if wraponly is not None:
        f2pycall.extend(["only:"] + wraponly + [":"])

    _buildcall(" ".join(f2pycall), env, verbose=verbose)

    ifortcall = ["ifort", "/c", "/MT", "/assume:buffered_io", "/names:lowercase",
                 "/assume:underscore"]
    if mode == "debug":
        ifortcall.extend(["/debug", "/traceback", "/check:bounds", "/check:stack",
                          "/warn:interfaces"])
    else:
        opt = "/" + opt if "/" not in opt else opt
        ifortcall.append(opt)

    if exflags is not None:
        exflagiter = exflags.split() if isinstance(exflags, str) else exflags
        for arg in exflagiter:
            if arg not in ifortcall:
                ifortcall.append(arg)

    ifortcall.extend(['/object:' + quote(tmppath), '/module:' + quote(tmppath)])
    ifortcall.extend(srcfiles)
    if os.path.isfile(tmppath + modulename + '-f2pywrappers.f') or \
            os.path.isfile(tmppath + modulename + '-f2pywrappers2.f90'):
        ifortcall.append(quote(tmppath + '*.f*'))
    _buildcall(" ".join(ifortcall), env, verbose=verbose)

    clcall = ["cl", "/MT"]
    clflags = ["/O2"]
    clinclude = ['-I' + quote(np.get_include()), '-I' + quote(PYPATH + "include/"),
                 '-I' + quote(PYPATH + 'lib/site-packages/numpy/f2py/src')]
    cl_1 = clcall + clflags + clinclude + \
        ['-c', quote(os.path.join(PYPATH, "Lib/site-packages/numpy/f2py/src/fortranobject.c")),
         '/Fo' + quote(tmppath + "fobj.obj")]
    cl_2 = clcall + clflags + clinclude + \
        ['-c', quote(tmppath + modulename + 'module.c'),
         '/Fo' + quote(tmppath + modulename + 'module.obj')]
    _buildcall(" ".join(cl_1), env, verbose=verbose)
    _buildcall(" ".join(cl_2), env, verbose=verbose)

    linkflags = ['/DLL', '/SUBSYSTEM:console', '/STACK:1512000000,1512000000', '/MACHINE:x64']
    if mode == "debug":
        linkflags.extend(["/DEBUG", "/pdb:link.pdb"])
    if includes is None:
        linkincludes = []
    else:
        linkincludes = [quote(srcdir + inc) for inc in includes]
    libpaths = ['/libpath:' + quote(os.path.join(PYPATH, "lib/site-packages/numpy/f2py/src"))]
    link_out = '/out:"' + outdir + modulename + '.pyd"'

    linkcall = ["link"] + linkflags + linkincludes + libpaths + \
        [quote(LIBPYTHON), quote(tmppath + "*.obj"), quote(tmppath + "*.o"), link_out]
    _buildcall(" ".join(linkcall), env, verbose=verbose)
    if os.path.isfile(quote(outdir + modulename + ".pyd")):
        _updatestatus("{} success!".format(modulename))

    # cleanup linker files:
    if os.path.isfile(outdir + modulename + '.exp'):
        try:
            os.remove(outdir + modulename + '.exp')
        except:
            pass
    if os.path.isfile(outdir + modulename + '.lib'):
        try:
            os.remove(outdir + modulename + '.lib')
        except:
            pass


def relinkintel(mode, exflags, modulename, sourcefiles, includes=None, opt="O2", wraponly=None,
                outdir="./", srcdir="./src/", tmpdir="./tmp/", env=None, verbose=True):
    """
    Relinks a single PYD assuming the required code is already compiled and only the link step needs
    to be performed.

    See ``buildintelf2py()`` for parameter defininitions.

    .. codeauthor:: Ryan Martin - 20-03-2018
    """
    tmppath = tmpdir + '{}'.format(modulename) + "/"
    ensure_path([outdir, tmpdir, tmppath])

    # intel's compiler embeds dependencies on the python library, so these need to be recompiled
    # to link
    clcall = ["cl", "/MT"]
    clflags = ["/O2"]
    clinclude = ['-I' + quote(np.get_include()), '-I' + quote(PYPATH + "include/"),
                 '-I' + quote(PYPATH + 'lib/site-packages/numpy/f2py/src')]
    cl_1 = clcall + clflags + clinclude + \
        ['-c', quote(os.path.join(PYPATH, "Lib/site-packages/numpy/f2py/src/fortranobject.c")),
         '/Fo' + quote(tmppath + "fobj.obj")]
    cl_2 = clcall + clflags + clinclude + \
        ['-c', quote(tmppath + modulename + 'module.c'),
         '/Fo' + quote(tmppath + modulename + 'module.obj')]
    _buildcall(" ".join(cl_1), env, verbose=verbose)
    _buildcall(" ".join(cl_2), env, verbose=verbose)

    linkflags = ['/DLL', '/SUBSYSTEM:console', '/STACK:1512000000,1512000000', '/MACHINE:x64']
    if mode == "debug":
        linkflags.extend(["/DEBUG", "/pdb:link.pdb"])
    if includes is None:
        linkincludes = []
    else:
        linkincludes = [quote(srcdir + inc) for inc in includes]
    libpaths = ['/libpath:' + quote(os.path.join(PYPATH, "lib/site-packages/numpy/f2py/src"))]
    link_out = '/out:"' + outdir + modulename + '.pyd"'

    linkcall = ["link"] + linkflags + linkincludes + libpaths + \
        [quote(LIBPYTHON), quote(tmppath + "*.obj"), quote(tmppath + "*.o"), link_out]
    _buildcall(" ".join(linkcall), env, verbose=verbose)
    if os.path.isfile(quote(outdir + modulename + ".pyd")):
        _updatestatus("{} success!".format(modulename))

    # cleanup linker files:
    if os.path.isfile(outdir + modulename + '.exp'):
        try:
            os.remove(outdir + modulename + '.exp')
        except:
            pass
    if os.path.isfile(outdir + modulename + '.lib'):
        try:
            os.remove(outdir + modulename + '.lib')
        except:
            pass


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
    f2pycall = [F2PY, '-m', modulename, "--build-dir", quote(tmppath)]
    if mode == "debug":
        f2pycall.append("--debug-capi")
    f2pycall.append(quote(srcfiles[-1]))
    if wraponly is not None:
        f2pycall.extend(["only:"] + wraponly + [":"])
    _buildcall(" ".join(f2pycall), env, verbose=verbose)

    if mode == "debug":
        opt = "-g -fbacktrace -fcheck=all"
    else:
        opt = "-" + opt if "-" not in opt else opt
    gfortbasecall = ["gfortran", "-c", opt, '-J' + quote(tmppath)]

    if exflags is not None:
        exflagiter = exflags.split() if isinstance(exflags, str) else exflags
        for arg in exflagiter:
            if arg not in gfortbasecall:
                gfortbasecall.append(arg)

    def localbuild(fl):
        _, flname = os.path.split(fl)
        flname, _ = os.path.splitext(flname)  # get the name of the file without ext
        localcall = copy.deepcopy(gfortbasecall)
        localcall.extend([quote(fl), "-o", quote(tmppath + flname + ".o")])
        _buildcall(" ".join(localcall), env, verbose=verbose)

    for fl in srcfiles:
        localbuild(fl)

    if os.path.isfile(tmppath + modulename + '-f2pywrappers.f'):
        localbuild(tmppath + modulename + '-f2pywrappers.f')
    if os.path.isfile(tmppath + modulename + '-f2pywrappers2.f90'):
        localbuild(tmppath + modulename + '-f2pywrappers2.f90')

    clcall = ["gcc"]
    clflags = ["-DMS_WIN64", "-Ofast"]
    clinclude = ['-I' + quote(np.get_include()), '-I' + quote(PYPATH + "include"),
                 '-I' + quote(os.path.join(PYPATH, "lib/site-packages/numpy/f2py/src"))]
    cl_1 = clcall + clflags + clinclude + \
        ['-c', quote(os.path.join(PYPATH, "Lib/site-packages/numpy/f2py/src/fortranobject.c")),
         '-o', quote(tmppath + "fobj.o")]
    cl_2 = clcall + clflags + clinclude + \
        ['-c', quote(tmppath + modulename + "module.c"),
            '-o', quote(tmppath + modulename + "module.o")]
    _buildcall(" ".join(cl_1), env, verbose=verbose)
    _buildcall(" ".join(cl_2), env, verbose=verbose)

    linkincludes = []
    if includes is not None:
        linkincludes.extend([quote(os.path.join(srcdir, inc)) for inc in includes])
    outname = ["-o", quote(os.path.join(outdir, modulename + ".pyd"))]
    libpaths = ['-L' + quote(os.path.join(PYPATH, "lib/site-packages/numpy/f2py/src"))]
    linkflags = ["-static", "-static-libgfortran", "-static-libgcc"]
    if mode == "debug":
        linkflags.insert(0, "-debug")
    files = [quote(fl) for fl in glob.glob(tmppath + "*.o")]
    linkcall = ["gfortran", "-shared"] + outname + libpaths + [quote(DLLPYTHON)] + \
        files + linkincludes + linkflags
    _buildcall(" ".join(linkcall), env, verbose=verbose)
    if os.path.isfile(quote(outdir + modulename + ".pyd")):
        _updatestatus("{} success!".format(modulename))


def relinkgnu(mode, exflags, modulename, sourcefiles, includes=None, opt="O2", wraponly=None,
              outdir="./", srcdir="./src/", tmpdir="./tmp/", env=None, verbose=True):
    """
    Relinks a single PYD assuming the required code is already compiled and only the link step needs
    to be performed.

    See ``buildgnuf2py()`` for parameter defininitions.

    .. codeauthor:: Ryan Martin - 20-03-2018
    """
    tmppath = tmpdir + '{}'.format(modulename) + "/"
    ensure_path([outdir, tmpdir, tmppath])

    linkincludes = []
    if includes is not None:
        linkincludes.extend([quote(srcdir + inc) for inc in includes])
    outname = ["-o", quote(outdir + modulename + ".pyd")]
    libpaths = ['-L' + quote(os.path.join(PYPATH, "lib/site-packages/numpy/f2py/src"))]
    linkflags = ["-static", "-static-libgfortran", "-static-libgcc"]
    if mode == "debug":
        linkflags.insert(0, "-debug")
    files = [quote(fl) for fl in glob.glob(tmppath + "*.o")]
    linkcall = ["gfortran", "-shared"] + outname + libpaths + [quote(DLLPYTHON)] + \
        files + linkincludes + linkflags
    _buildcall(" ".join(linkcall), env, verbose=verbose)
    if os.path.isfile(quote(outdir + modulename + ".pyd")):
        _updatestatus("{} success!".format(modulename))


def build_lapack(compiler, force_compile=False, verbose=True):
    """A small utility to build the lapack library that is required in some of the .pyd builds.

    Parameters
    ----------
        compiler: str
            'intel' or 'gnu'
        force_compile: bool
            If `True` will force compile even if lapack library is already built
        verbose: bool
            'all' prints all output, 'min' only prints a message if it starts compiling
            lapack, `none` prints only if error is found

    .. codeauthor:: Tyler Acorn March 08, 2017
    """
    lapack_folder = FILEDIR + '/resource/'
    script = 'compile_lapack.py'
    if compiler.lower() == 'intel':
        if not os.path.isfile(lapack_folder + 'lapack_solve_intel.lib') or force_compile:
            if verbose.lower() in ['all', 'min']:
                print('Compiling Lapack for Intel')
            lapackcall = ["python", script, compiler.lower()]
            _buildcall(lapackcall, os.environ, verbose=verbose, builddir=lapack_folder)
        # Check to make sure the library is there
        if not os.path.isfile(lapack_folder + 'lapack_solve_intel.lib'):
            _updatestatus('Error: lapack library for INTEL was not build')
            sys.exit(1)
    elif compiler.lower() == 'gnu' or os.path.isdir(compiler):
        if not os.path.isfile(lapack_folder + 'lapack_solve.a') or force_compile:
            if verbose.lower() in ['all', 'min']:
                print('Compiling Lapack for GNU')
            lapackcall = ["python", script, compiler.lower()]
            _buildcall(lapackcall, os.environ, verbose=verbose, builddir=lapack_folder)
        # Check to make sure the library is there
        if not os.path.isfile(lapack_folder + 'lapack_solve.a'):
            _updatestatus('Error: lapack library for GNU was not build')
            sys.exit(1)
    else:
        raise Exception('choose either "intel" or "gnu" compiler')


def _buildcall(buildstr, env, verbose=True, builddir=None):
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
    if ' error ' in fulloutput or 'Error:' in fulloutput:
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


def extensioncleaner(nametoremove, tmpdir, outdir):
    """ Clean dll's, pyd's and temp directories """
    for ext in ["dll", "pyd"]:
        if nametoremove == 'all':
            if ext == "dll":
                built_ext = "{}lib*.dll".format(outdir)
            else:
                built_ext = "{}*.pyd".format(outdir)
            for built_ext in glob.glob(built_ext):
                try:
                    os.remove(built_ext)
                except:
                    try:
                        call([RMCALL, built_ext])
                    except:
                        raise Exception('The file: ' + built_ext + ' is in use, and rm.exe '
                                        'was not found on the path. Close the file before cleaning '
                                        'or make sure rm.exe is callable from the system variable '
                                        'paths')
            if os.path.isdir(tmpdir):
                try:
                    shutil.rmtree(tmpdir)
                except:
                    raise Exception('The folder {} is in use, close command prompts or files '
                                    'that might be open in this folder!'.format(tmpdir))
        else:
            if ext == "dll":
                built_ext = "{}lib{}.dll".format(outdir, nametoremove)
            else:
                built_ext = "{}{}.pyd".format(outdir, nametoremove)
            if os.path.isfile(built_ext):
                try:
                    os.remove(built_ext)
                except:
                    try:
                        call([RMCALL, built_ext])
                    except:
                        raise Exception('The file: ' + nametoremove + ' is in use, and rm.exe '
                                        'was not found on the path. Close the file before cleaning '
                                        'or make sure rm.exe is callable from the system variable '
                                        'paths')
            if os.path.isdir(tmpdir):
                if tmpdir == "./" or tmpdir == "":
                    print("Current directory `./` will not be cleared!")
                    return
                try:
                    shutil.rmtree(tmpdir)
                except:
                    raise Exception('The folder tmp/ is in use, close command prompts or files '
                                    'that might be open in this folder!')


def assert_environment(compiler):
    """
    Ensure that the specified compiler(s) can be found. If compiler == `'gnu'`, first check if a
    gfortran compiler is found on the path. If yes, then we're using that one (and associated gcc
    that should also be found on the path).. If no, check a standard location for MinGW with the
    Anaconda python distribution. If that fails then an attempt at installing the required compilers
    is made. Helpful output is given if nothing passes.
    """
    if compiler == "gnu":
        env = None
        try:
            call("gfortran", stderr=PIPE, stdout=PIPE)
            if PYVERSION == "27":
                p = Popen("gfortran -v", stderr=PIPE, stdout=PIPE)
                gfortinfo = p.communicate()[0]
            else:
                from subprocess import getoutput
                gfortinfo = getoutput("gfortran -v")
            if "cygwin" in gfortinfo:  # currently binaries built with cygwin gfortran are no good..
                print("WARNING: found `cygwin` compilers which are not currently supported.. "
                      "looking for MinGW from Anaconda")
                raise FileNotFoundError
        except FileNotFoundError:
            if not os.path.isdir(PYPATH + "MinGW/bin"):
                print("WARNING: could not find `mingw` ... installing ... ")
                try:
                    _buildcall("conda install -c msys2 m2w64-gcc-fortran -y", os.environ)
                except FileNotFoundError:
                    _buildcall("pip install mingw -y", os.environ)
                except FileNotFoundError:
                    raise ValueError("Could not automatically install mingw. Consider installing "
                                     "a gfortran compiler with tools such as mingw, pip "
                                     "or conda")
                # if we made it here, mingw is installed, but the bin is not on the path, so we
                # inject it to the top of the path so it is found first
                os.environ["PATH"] = PYPATH + "MinGW/bin;" + os.environ["PATH"]
            else:
                print("FOUND: `mingw` compilers at " + PYPATH + "MinGW/bin \n"
                      "consider adding this folder to the path!")
                os.environ["PATH"] = PYPATH + "MinGW/bin;" + os.environ["PATH"]
            # one final assert to make sure gfortran is callable
            try:
                call("gfortran", stderr=PIPE, stdout=PIPE)
            except FileNotFoundError:
                raise ValueError("Could not automatically install mingw. Consider installing "
                                 "a gfortran compiler with tools such as cygwin, mingw")
        else:
            try:
                call("gcc", stderr=PIPE, stdout=PIPE)
            except FileNotFoundError:
                raise ValueError("Found a fortran compiler, but couldnt automatically find an "
                                 "associated `gcc` compiler... Consider installing "
                                 "compiler tools with cygwin or mingw")
        env = os.environ.copy()
    elif compiler == "intel":  # the compiler was 'intel'
        try:
            env = load_ifortenv()
        except ValueError:
            raise ValueError("Could not automatically configure environment variables for `intel`"
                             "compilers! Consider `gfortran`")
    elif os.path.isdir(compiler):  # a path to a gfortran compiler was passed
        os.environ["PATH"] = compiler + ";" + os.environ["PATH"]
        env = os.environ.copy()
        # assert that this indeed provided a gfortran compiler
        try:
            call("gfortran", stderr=PIPE, stdout=PIPE)
        except FileNotFoundError:
            raise ValueError("ERROR: the specified path {} did not provide the required gfortran "
                             "compilers. Consider using one of `gnu` or `intel` instead ")
    else:
        raise ValueError("ERROR: did not understand compiler `{}`".format(compiler))
    return env


def build_custom_fortran(sources, includes={}, wraponly={}, name=None, srcdir="./",
                         outdir="./", tmpdir="./tmp/", compiler="gnu", mode="release", exflags="",
                         opt="O2", clean=True, verbose=True):
    """
    This function is intended to allow arbitrary f2py extensions to be constructed using the tooling
    provided in this module. This function replaces `FortranBuild` as a more succinct methodology
    to compile the requred sources.

    Parameters:
        sources (list or dict): either a list of fortran files where the last depends on the first
            and the last file in the list contains the f2py wrapping code, or a dictionary where
            keys in the dictionary indicate the name of the module to be built and the values are the
            lists of associated sources to generate that module. See ``pygeostat/fortran/sourcedefs.py``
            for inspiration on the structure of these dictionaries
        includes (list or dict): a matching item to `sources` that contains a list or dict of extra
            items to include on the link step of compilation. See ``pygeostat/fortran/sourcedefs.py``.
        wraponly (list or dict): matching item to ``sources`` and ``includes`` that contains a list or
            dictionary of functions that get wrapped for python, other functions in the f2py fortran
            code are ignored. See ``pygeostat/fortran/sourcedefs.py``.
        name (str): if ``sources`` is a list, a name must be provided for the resulting ``<name>.pyd``
        srcdir, tmpdir, outdir (str): various locations to consider
        compiler (str): either ``"intel"``, ``"gnu"`` or a path to a compiler ``bin`` directory
        mode (str): ``"release"`` or ``"debug"``
        exflags (str): compiler-specific permissable fortran compiler flags
        opt (str): optimization level, defaults to ``O2``
        clean (bool): if ``True`` then the build directory is cleaned before compiling (recommended)
        verbose (bool): write all output to the terminal

    .. codeauthor:: Ryan Martin - 07-04-2018
    """
    buildfunc = buildintelf2py if compiler == "intel" else buildgnuf2py
    env = assert_environment(compiler)
    workingdir = os.getcwd()

    if isinstance(sources, (list, tuple)):
        assert name is not None, "Must pass a `name` if sources is a list"
        sources = {name: sources}
    else:
        assert isinstance(sources, dict)
    if isinstance(includes, (list, tuple)):
        assert name is not None, "Must pass a `name` if includes is a list"
        includes = {name: includes}
    else:
        assert isinstance(includes, dict)
    if isinstance(wraponly, (list, tuple)):
        assert name is not None, "Must pass a `name` if wraponly is a list"
        wraponly = {name: wraponly}
    else:
        assert isinstance(wraponly, dict)

    tmpdir = os.path.join(workingdir, tmpdir)
    outdir = os.path.join(workingdir, outdir)

    _updatestatus("Building ... ")
    for modname, codelist in sources.items():
        if clean:
            extensioncleaner(modname, tmpdir, outdir)
        if srcdir is not None:
            codelist = [os.path.join(srcdir, cf) for cf in codelist]
        incl = includes.get(modname, None)
        if srcdir is not None and incl is not None:
            incl = [os.path.join(srcdir, ifl) for ifl in incl]
        wraps = wraponly.get(modname, None)
        buildfunc(mode, exflags, modname, codelist, includes=incl, opt=opt, wraponly=wraps,
                  srcdir=srcdir, outdir=outdir, tmpdir=tmpdir, env=env, verbose=verbose)

    if all([os.path.isfile(os.path.join(outdir, m + ".pyd")) for m in sources.keys()]):
        _updatestatus("Successfully built required fortran extensions!")


def build_pygeostat_fortran(pygeostatmodules="all", compiler="gnu", mode="release", exflags="",
                            opt="O2", clean=True, verbose=False):
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
    from .sourcedefs import sources
    allpygeostat = list(sources.keys())
    assert pygeostatmodules in ["all"] + allpygeostat, (
        "{} is not in Pygeostat!".format(pygeostatmodules))
    assert_environment(compiler)
    # regardless of where the function is called the build is called here.

    compilerargs = ["-compiler={}".format(compiler), "-mode={}".format(mode), "-opt={}".format(opt)]
    if exflags != "":
        exflags = ['-exflags="{}"'.format(" ".join(exflags.split()))]
    else:
        exflags = []
    compilescript = "compile.py"
    # first clean the build
    if clean:
        _updatestatus("Cleaning ... ")
        cleanstr = " ".join(["python", compilescript, "-clean", pygeostatmodules])
        _buildcall(cleanstr, None, verbose, FILEDIR)

    _updatestatus("Building ... ")
    # call to build the specified module
    buildstr = " ".join(["python", compilescript, pygeostatmodules] + compilerargs + exflags)
    _buildcall(buildstr, None, verbose, FILEDIR)

    tocheck = allpygeostat if pygeostatmodules == "all" else [pygeostatmodules]
    isfound = [os.path.isfile(os.path.join(FILEDIR, modname + ".pyd")) for modname in tocheck]
    if all(isfound):
        _updatestatus("Successfully built required fortran extensions!")
    else:
        failedstr = ["FAILED BUILDING: {}".format(name)
                     for name, found in zip(tocheck, isfound) if not found]


if __name__ == '__main__':
    import argparse
    from subprocess import call
    from sourcedefs import sources, intelincludes, gnuincludes, wraponly

    allkeys = sources.keys()
    parser = argparse.ArgumentParser(
        description=('Process the build calls to this function. '
                     'Valid modules are: {}').format(', '.join([str(key) for key in allkeys])))
    parser.add_argument('modulename', type=str, nargs='?',
                        help=('the module to be built, use `list` to list'
                              ' permissable module names'))
    parser.add_argument('-compiler', type=str, default='gnu',
                        help='Set the compiler to one of `intel` or `gnu`')
    parser.add_argument('-clean', type=str, default=None,
                        help='Clean the build directory, either `all` or one of `list`')
    parser.add_argument('-mode', type=str, default='release',
                        help=('Valid are default is `release`, or `debug`'))
    parser.add_argument('-fflags', nargs="+", type=str,
                        help='String of extra fortran flags for compiling, e.g for intel: `"/check:bounds '
                             '/check:stack"` or `"-flag"` for gnu compilers ')
    parser.add_argument('-tmpdir', type=str, default='tmp/',
                        help=('Directory for temporary output'))
    parser.add_argument("-opt", type=str, default="O2", help=("The optimization level"))
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
    if args.clean is not None:
        extensioncleaner(args.clean, args.tmpdir, args.outdir)
    elif args.modulename is None:
        parser.print_help()
    elif args.modulename == "list":
        print('POSSIBLE KEYS: \nall  (use to build all modules)\n%s\n' %
              ' \n'.join([key for key in sorted(allkeys)]))
    else:
        if args.modulename != "all" and args.modulename not in allkeys:
            parser.print_help()
            _updatestatus("ERROR: \n" + args.modulename + ' is not in the default build dicitonary! '
                          'Try `all` or one of: \n%s\n ' % ', '.join([key for key in sorted(allkeys)]))
            sys.exit(0)
        # parse the dictionaries for intel or gnu compilers
        if args.compiler == "intel":
            includes = intelincludes
            buildfunc = relinkintel if args.relink else buildintelf2py  # buildinteldll if args.out == "dll" else
            _ = assert_environment("gnu")
            env = assert_environment("intel")
        elif args.compiler == "gnu" or os.path.isdir(args.compiler):
            includes = gnuincludes
            buildfunc = relinkgnu if args.relink else buildgnuf2py  # buildgnudll if args.out == "dll" else
            env = assert_environment(args.compiler)
        else:
            print("specfied compiler or compiler path not understood: {}".format(args.compiler))
        # assert lapack is built after asserting the environment
        build_lapack(args.compiler, args.forcelapack, verbose="all")

        buildkeys = allkeys if args.modulename == "all" else [args.modulename]
        # actually call the build functions
        for buildname in buildkeys:
            sourcecode = sources[buildname]
            include = includes.get(buildname, None)
            towrap = wraponly.get(buildname, None) if args.wraponly is None else args.wraponly
            buildfunc(args.mode, args.fflags, buildname, sourcecode, opt=args.opt,
                      wraponly=towrap, includes=include, tmpdir=args.tmpdir, env=env)
        os.environ = origenv
