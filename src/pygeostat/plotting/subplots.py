#!/usr/bin/env python
# -*- coding: utf-8 -*-
# public
"""imagegrid.py: Provides a means to easily generate grids of plots"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from . set_style import set_plot_style


@set_plot_style
def subplots(nrows=1, ncols=1, figsize=None, axes_pad=(0.9, 0.4), aspect=True, label_mode='L',
             cbar_mode=None, cbar_location='right', cbar_pad=.03, cbar_size='3%', cbar_set_cax=True,
             share_all=True, fig=None, axes_class=None, rect=111):
    '''
    A wrapper of the Matplotlib ImageGrid class, providing additional ease of use for the following
    reasons:

        - The format of the arguments is modified to mimic the more commonly used
          matplotlib.pyplot.subplots
        - A matplotlib.figure.Figure no longer needs to be initiated seperately
        - Kwarg locations and defaults are modified to emphasize more commonly used
          and modified settings. Rarely used kwargs such as direction and ngrids are removed.
        - The involved import statement of ImageGrid is avoided
        - Project defaults are integrated via pygeostat.Parameters (future work)

    ImageGrid:

        "A class that creates a grid of Axes. In matplotlib, the axes location (and size) is specified
        in the normalized figure coordinates. This may not be ideal for images that needs to be
        displayed with a given aspect ratio.  For example, displaying images of a same size with some
        fixed padding between them cannot be easily done in matplotlib. ImageGrid is used in such case."

    Parameters:
        nrows(int): number of of rows in the grid
        ncols(int): number of of columns in the grid

    Optional keyword arguments (from ImageGrid):

        ================  ========  =========================================
        Keyword           Default   Description
        ================  ========  =========================================
        axes_pad          0.02      float| pad between axes given in inches
                                    or tuple-like of floats,
                                    (horizontal padding, vertical padding)
        aspect            True      [ True | False ] If True, each the length
                                    of the x and y axes are based on their
                                    absolute range. Generally used for plots
                                    such as maps, but should be False for
                                    plots such as variograms.
        label_mode        "L"       [ "L" | "1" | "all" ]
        cbar_mode         None      [ "each" | "single" | "edge" ]
        cbar_location     "right"   [ "left" | "right" | "bottom" | "top" ]
        cbar_pad          None
        cbar_size         "5%"
        cbar_set_cax      True      [ True | False ] if True, each axes in
                                    the grid has a cax attribute that is bind to
                                    associated cbar_axes.
        share_all         False     [ True | False ]
        axes_class        None      a type object which must be a subclass
                                    of axes_grid's subclass of
                                    :class:`~matplotlib.axes.Axes`
        ================  ========  =========================================

    Keyword Arguments:
        figure(matplotlib.figure.Figure): a figure is created if one is not
            provided (based on rcParams['figure.figsize'] if None)
        figsize(matplotlib.figure.Figure): size of the figure, if one must
            be created. Based on rcParams['plotting.figsize'] if None
        rect(int): "[left, bottom, width, height]* (in:class:`
            ~matplotlib.figure.Figure` coordinates) or the subplot position code
            (e.g., "121"). The default setting occupies the entire figure.

    Returns:
        (fig, axes): mirroring the output of matplotlib.subplots()

    .. seealso::
        :func:`Subplots Clean <pygeostat.plotting.subplots.subplots_clean>`


    **Example:**

    Plot 4 realizations in each panel with a common colorbar.

    .. plot::

        import pygeostat as gs

        # Initialize gs.subplots
        fig, axes = gs.subplots(2, 2, figsize=(8, 8), cbar_mode='single')

        # Default grid definition
        gs.Parameters['data.griddef'] = gs.GridDef('120 5 10\\n110 1205 10\\n1 0.5 1.0')

        # Iterate over the axes/realizations
        for i, ax in enumerate(axes):
            sim = gs.ExampleData('grid2d_surf_real'+str(i+1))
            gs.slice_plot(sim, var='Top Elevation', ax=ax, vlim=(375, 385),
                        title='Realization '+str(i+1))


    .. codeauthor:: pygeostat development team 2018-04-05
    '''
    from mpl_toolkits.axes_grid1 import ImageGrid
    from matplotlib.pyplot import figure
    if fig is None:
        fig = figure(figsize=figsize)
    if rect is None:
        rect = 111
    axes = ImageGrid(fig, rect, (nrows, ncols),
                     axes_pad=axes_pad, share_all=share_all,
                     aspect=aspect, label_mode=label_mode,
                     cbar_mode=cbar_mode, cbar_location=cbar_location,
                     cbar_pad=cbar_pad, cbar_size=cbar_size,
                     cbar_set_cax=cbar_set_cax, axes_class=axes_class)
    return (fig, axes)


def subplots_clean(axes, nused, ncols=None):
    '''Remove unused axes. Add labels where space is made available due to
       the removed axes.

       Parameters:
        axes (list of axes from gs.subplots or plt.subplots): where the
            number of axes exceeds the number used (nused)
        nused (int):number of used axes, where it assumed that axes[:nused]
            are used.
        ncols(int): number of of columns in the grid, which is required to
            add labels.
    '''
    nax = len(axes)
    if nused >= nax:
        raise ValueError('nused >= len(axes)! nothing can be done!')
    # Remove the unrequired axes
    for ax in axes[nused:]:
        ax.remove()
    if ncols is not None:
        if nax - nused > ncols:
            raise ValueError('len(axes) - nused > ncols! use a smaller grid!')
    # Add xlabels where required
    for i, ax in enumerate(axes[:nused]):
        if i + ncols >= nused:
            ax.axis['bottom'].toggle(ticklabels=True, label=True)
