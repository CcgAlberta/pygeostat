import sys
sys.path.append("..")  # for some reason I cannot relative import this
from compile import assert_environment, _updatestatus, _buildcall


def build_lapack_gnu():
    env = assert_environment("gnu")
    # Compiler Flags
    flflags = 'LAPACK/all/*.f LAPACK/eig.f90 LAPACK/solve.f90 '
    lapackcall = "gfortran -c -O2".split() + flflags.split()
    _buildcall(" ".join(lapackcall), env=env, verbose=True)
    # Make the .lib file:
    _buildcall("ar rvs lapack_solve.a *.o", env=env)
    # Cleanup:
    _buildcall("./rm/rm.exe *.o", env=env)


def build_lapack_intel():
    _ = assert_environment("gnu")  # this is now required since ar.exe is needed
    env = assert_environment("intel")
    # Compiler Flags
    flflags = 'LAPACK/all/*.f LAPACK/eig.f90 LAPACK/solve.f90 '
    lapackcall = "ifort /c /names:lowercase /assume:underscore /O2".split() + flflags.split()
    _buildcall(" ".join(lapackcall), env=env, verbose=True)
    # Make the .lib file:
    _buildcall("ar rvs lapack_solve_intel.lib *.obj", env=env)
    # Cleanup:
    _buildcall("./rm/rm.exe *.obj", env=env)


if __name__ == '__main__':
    if any("intel" in a for a in sys.argv):
        _updatestatus("Building LAPACK for intel")
        build_lapack_intel()
    else:
        _updatestatus("Building LAPACK for gnu")
        build_lapack_gnu()
