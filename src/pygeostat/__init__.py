'''
Python interface for geostatistics
==================================

pygeostat is a Python module for geostatistics including interfacing
with GSLIB (Geostatistical Software Library, Copyright 1992-2009
Clayton Deutsch and Andre Journel) and CCG (Centre for Computational
Geostatistics, Copyright 1998-2014 individual authors) software.
It aims to provide a very efficient method to run, script and
interact with GSLIB and CCG software using Python.
'''

__version__ = '1.2.0'

# parameter configurations
from . pygeostat_parameters import Parameters, PlotStyle

# Load necessary modules
from .data import *
from .programs import *
from .plotting import *
from .datautils import *
from .transformations import *
from .multivariate import *
from .statistics import *
from .utility import *
