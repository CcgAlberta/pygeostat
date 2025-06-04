#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Function that alters the matplotlib rc parameters as desired with some predefined styles"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from functools import wraps
import matplotlib as mpl
from matplotlib import pyplot as plt
from ..pygeostat_parameters import PlotStyle


def set_plot_style(plot_function):
    """
    Applies pygeostat plotting style before plotting function is called based on the provided style name and custom dictionary then
    retore it back to the original plot style.

    **Example:**

    .. code-block:: Python

        import pygeostat as gs

        @gs.set_plot_style
        def coolnewplot(data, griddef, **otherkwargs, plot_style=None, custom_style=None):
            # plot some things
            # apply cool formatting
            return things
    """
    @wraps(plot_function) # Using the wraps for invoking update_wrapper(). This keeps the name and doc string of the original function
    def wrapper_plot_function(*args, **kwargs):
        # set the passed style
        PlotStyle.update_mplrcParams(kwargs.get("plot_style", None), 
                                       kwargs.get("custom_style", None))
        # call the plotting function
        returnargs = plot_function(*args, **kwargs)
        # plt.show()
        # set the style back to whatever the default was
        PlotStyle.restore_mplrcParams()
        return returnargs
    
    return wrapper_plot_function
