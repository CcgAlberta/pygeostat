:: This batch file compiles Fortran code used by pygeostat into a pyd code
:: Requires - f2py, gcc and gfortran/gnu95 compilers

:: Note that on 64 bit Windows machines, the following lines in:
:: Updated version of {python-dir}\Lib\site-packages\numpy\distutils\fcompiler\gnu.py
:: must be changed - note the addition of a pass and commenting out
::            # XXX: fix this mess, does not work for mingw
::            if is_win64():
::                c_compiler = self.c_compiler
::                if c_compiler and c_compiler.compiler_type == "msvc":
::                    return []
::                else:
::                    pass #raise NotImplementedError("Only MS compiler supported with gfortran on win64")

call python compile.py -clean all
call python compile.py -compiler=gnu all
:: call f2py -c -m --compiler=mingw32 --fcompiler=gnu95 covasubs covasubs.for
:: call f2py -c -m --compiler=mingw32 --fcompiler=gnu95 dgm dgm.for
:: call f2py -c -m --compiler=mingw32 --fcompiler=gnu95 histcorrect random.f90 normaldist.f90 sortem.f90 histcorrect.f90
:: call f2py -c -m --compiler=mingw32 --fcompiler=gnu95 gausskde gausskde.f90
:: call f2py -c -m --compiler=mingw32 --fcompiler=gnu95 spbusim random.f90 covasubs.for sortem.f90 normaldist.f90 gausskde.f90 linearinterp.f90 spbusim.f90 ../../resource/lapack_solve.a
:: call f2py -c -m --compiler=mingw32 --fcompiler=gnu95 --opt='-O2' varmodel covasubs.for random.f90 varmodel.f90
:: call f2py -c -m --compiler=mingw32 --fcompiler=gnu95 --opt='-O2' varcalc sortem.f90 varsubs.f90 varcalc.f90
:: call f2py -c -m --compiler=mingw32 --fcompiler=gnu95 getcollocated getcollocated.f90 sortem.f90