Building Pygeostat Fortran Extensions (windows only):
-----------------------------------------------------
From Python:
 >  ipython
 $  import pygeostat as gs
 $  gs.build_pygeostat_fortran()

From commandline:
- open a cmd prompt in `/pygeostat/fortran/`
- run:
 >  python compile.py -compiler=gnu pygsb


Alternatively:
--------------
- download updated mingw w64 compilers (e.g., `x86_64-7.3.0-release-posix-seh-rt_v5-rev0.7z`)
  from sourceforge (here: https://goo.gl/U228Q3)
- copy contents of `x86_64-7.3.0-release-posix-seh-rt_v5-rev0` somewhere, e.g. create 
  a folder `C:/mingw` and extract the contents into `mingw/` so that `C:/mingw/bin` exists
- open a cmd prompt in `pygeostat/fortran/` (shift + right click > open cmd here)

 >  python compile.py all -compiler="C:/mingw/bin"

OR start python:

 > ipython 

and run:

 $ import pygeostat; pygeostat.build_pygeostat_fortran(compiler="C:/mingw/bin")


Recommended:
------------
Add `C:/mingw/bin` to your environment variables. Confirm this compiler is found by 
inspecting the output of `where [or which] gfortran` from a terminal. (Note: restart
open terminals after adding `C:/mingw/bin` to the path). Original methods will utilize 
the updated compilers: 

 >  python compile.py all -compiler=gnu  :: from cmd in /pygeostat/fortran/

 $  import pygeostat; pygeostat.build_pygeostat_fortran()  # from python, anywhere


NOTE:
-----
Cygwin compilers are able to compile the pygeostat fortran dll's, however, something in
the linking toolchain is fundamentally incompatible with Anaconda Python. Thus, if a 
`cygwin` compiler is detected, the build script reverts to the Anaconda MinGW installed 
using `conda install mingw`. Using the updated compiler set above is recommended instead 
of cygwin.