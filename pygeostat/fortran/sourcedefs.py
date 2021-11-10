# !/usr/bin/env python
#  -*- coding: utf-8 -*-

# Dictionary of source code -- path is relative to /pygeostat/fortran/src/

# Format:  { key (resulting .pyd filename) : [list of sourcecodes in order of dependence (last
#                                             depends on first)] }
# i.e. : { 'modulename' : ['list.f90','of.f90','source.f90','f2pycode.f90'] } -> generates
# modulename.pyd
sources = {
    "covasubs": ["covasubs.for"],
    "dgm": ["dgm.for"],
    "histcorrect": [
        "subs/random.f90",
        "subs/normaldist.f90",
        "subs/sortem.f90",
        "histcorrect.f90",
    ],
    "gausskde": ["gausskde.f90"],
    "varcalc": [
        "subs/sortem.f90",
        "subs/random.f90",
        "subs/varsubs.f90",
        "varcalc.f90",
    ],
    "varmodel": ["subs/random.f90", "covasubs.for", "varmodel.f90"],
    # "spbusim": [
    #     "subs/random.f90",
    #     "subs/sortem.f90",
    #     "subs/normaldist.f90",
    #     "linearinterp.f90",
    #     "covasubs.for",
    #     "gausskde.f90",
    #     "spbusim.f90",
    # ],
    "subsample": ["subs/random.f90", "subsample.f90"],
    "varsim": [
        "subs/fileprocessing.f90",
        "subs/gslib_binary.f90",
        "subs/varsim.for",
        "varsim_wrap.f90",
    ],
    "getcollocated": ["subs/sortem.f90", "getcollocated.f90"],
    "fgslib_io": ["fgslib_io.f90"],
    "pygsb": ["subs/fileprocessing.f90", "subs/gslib_binary.f90", "pygsb.f90"],
    "nscore": [
        "subs/quicksort.f90",
        "subs/acorni.f90",
        "subs/normaldist.f90",
        "subs/nscore_module.f90",
        "nscore.f90",
    ],
    "backtr": ["backtr.f90"],
    "esri_io": ["esri_io.f90"],
    "supersec": ["supersec.f90"],
    "sgvertices": ["sgvertices.f90"],
}
# Dictionary of Linked Libraries, make sure key matches dictionary above
# TODO: unify the naming such that gnu / intel compiler doesnt matter here
intelincludes = {"spbusim": ["./resource/lapack_solve_intel.lib"]}
gnuincludes = {"spbusim": ["./resource/lapack_solve.a"]}

# a dictionary indicating which functions from the f2py wrapper (last file in the build) should be
# exposed to python.
wraponly = {
    "histcorrect": ["histcorrect"],
    "sgvertices": ["calc_sgvertices"],
}
