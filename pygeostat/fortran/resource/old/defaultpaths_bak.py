# !/usr/bin/env python
#  -*- coding: utf-8 -*-
# public
"""
This script contains the compiler paths and things needed for the intel and gnu compilers
to build all the specified modlues
"""
import sys
import os

# ------------------------------------------------------------------------------------------------#
# THINGS THAT ARE 'AUTOMATIC'....
# ------------------------------------------------------------------------------------------------#
# Path to the current directory where this script resides
cdir = os.path.dirname(os.path.realpath(__file__))
# Path to the pygeostat folder
pygeostatfolder = cdir[0:cdir.rfind('pygeostat')].replace("\\", '/')
# Automatically determine the python version from the calling script (compile.py)
pyversion = str(sys.version_info[0]) + str(sys.version_info[1])
# Path to Python distribution
pypath = sys.executable[0:sys.executable.rfind('python.exe')].replace('\\', '/')
# Path to the Python[ver].lib for the intel compiler
libpython = pypath + 'libs/Python' + pyversion + '.lib'
# Path tot he Python[ver].dll for the gnu compiler
dllpython = pypath + 'python' + pyversion + '.dll'
# Exact path to the gcc.exe from MinGW
cc = pypath + 'MinGW/bin/gcc.exe'

# ------------------------------------------------------------------------------------------------#
# SEMI-AUTOMATIC STUFF? #
# ------------------------------------------------------------------------------------------------#
# Local path where local build files are generated (can be defined in python)
tdir = 'tmp/'
# Path to where source fortran codes are stored (can be defined in python)
srcdir = pygeostatfolder + 'pygeostat/fortran/src/'
# ------------------------------------------------------------------------------------------------#

# ------------------------------------------------------------------------------------------------#
# INTEL SECTION - defaults for intel compiler 14 or 16, vs2012 and sdk7.1
# ------------------------------------------------------------------------------------------------#
# Path to the SDK that is correct for the VS version, with kernel32.lib and some other msflib!
sdkpath = 'C:/Program Files/Microsoft SDKs/Windows/v7.1/Lib/x64/'

# -------------------- #
# SEMI-AUTOMATIC STUFF #
# -------------------- #
# The following three paths can be set to a value of `auto`; however, if it doesn't work you'll
# have to set them manually
fcomp_dir = 'auto'   # Default intel compiler directories, where /bin and /mkl folder are
vccomp_dir = 'auto'  # VC directory, where /bin folder is. Leave as a
ifort = 'auto'       # Exact path to the ifort.exe within the `fcomp_dir` directory
link = 'auto'        # Exact path to the link.exe within the `vccomp_dir` directory
h5libs_dir = 'auto'  # Directory containing /lib and /static h5 folders. Looks for /VS{ver}_x{bit}

# ------------------------------------------------------------------------------------------------#
# GNU SECTION - defaults tested for MinGW installed through Anaconda
# ------------------------------------------------------------------------------------------------#
# MinGW file path
mingwdir = pypath + 'MinGW/'
# Exact path to the gnufortran.exe
fort_comp = mingwdir + 'bin/gfortran.exe'
# ------------------------------------------------------------------------------------------------#
