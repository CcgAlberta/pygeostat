# !/usr/bin/env python
#  -*- coding: utf-8 -*-
# public

import sys
import glob
sys.path.append("..")  # for some reason I cannot relative import this
from compile_linux import _updatestatus, _buildcall


def build_lapack_gnu():
    # Compiler Flags
    sourcefiles = glob.glob("./LAPACK/ALL/*.f")
    wrappers = 'LAPACK/eig.f90 LAPACK/solve.f90 '
    allfiles = sourcefiles + wrappers.split()
    lapackcall = "gfortran -c -O2".split() + allfiles
    _buildcall(lapackcall, verbose=False)
    # Make the .lib file:
    _buildcall("ar rvs lapack_solve.a".split() + allfiles)
    # Cleanup:
    _buildcall(["rm"] + glob.glob("./*.o"))


# def build_lapack_intel():
#     _ = assert_environment("gnu")  # this is now required since ar.exe is needed
#     env = assert_environment("intel")
#     # Compiler Flags
#     flflags = 'LAPACK/all/*.f LAPACK/eig.f90 LAPACK/solve.f90 '
#     lapackcall = "ifort /c /names:lowercase /assume:underscore /O2".split() + flflags.split()
#     _buildcall(" ".join(lapackcall), env=env, verbose=True)
#     # Make the .lib file:
#     _buildcall("ar rvs lapack_solve_intel.lib *.obj", env=env)
#     # Cleanup:
#     _buildcall("./rm/rm.exe *.obj", env=env)


if __name__ == '__main__':
    if any("intel" in a for a in sys.argv):
        _updatestatus("ERROR: this extension cannot compile with intel on linux!")
        # build_lapack_intel()
    else:
        _updatestatus("Building LAPACK for gnu")
        build_lapack_gnu()
