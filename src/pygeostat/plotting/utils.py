#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A collection of utility functions for plotting"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import copy
import gc

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

from .. pygeostat_parameters import Parameters

def titleoverlap(ax, title):
    """
    A small utility function called from many of the plotting functions. This returns a boolean if
    the title overlaps with x-axis tick labels that have been placed on the
    top spline.

    Parameters:
        ax (ax): A matplotlib axis object with x axis tick labels
        title (text object): A text object created by plt.title(), plt.suptitle(), or plt.text()
            (e.g., titletxt = plt.suptitle('Title', y=0.98))

    Returns:
        overlap (bool): Indicator if the x-axis tick labels overlap the plot title.

    .. codeauthor:: pygeostat development team - 2015-10-07
    """
    # Set-up a function that completes the check
    def ticktitlecheck(labels, title):
        """Check to see if the x-axis tick labels overlap the plot title"""
        fig = plt.gcf()
        renderer = fig.canvas.get_renderer()
        fig.draw(renderer)
        if not labels or not title:
            return False
        else:
            try:
                xboxes = [label.get_window_extent() for label in labels]
                titlebox = title.get_window_extent()
                overlaps = titlebox.count_overlaps(xboxes)
                return overlaps >= 1
            except:
                print("There was an error checking for overlapping x axis tick labels and plot",
                      "title. You'll have to set the argument manually to rotate them if necessary")
    overlaps = ticktitlecheck(ax.get_xticklabels(), title)

    return overlaps


def supaxislabel(axis, label, rotation=None, label_prop=None, labelpad=0, ha='center', va='center',
                 fig=None):
    """
    A small utility function called from many of the plotting functions. This adds super ylabel or
    xlabel to the figure similar to mpl.suptitle when using subplots.
    From: http://stackoverflow.com/a/29107972

    Parameters:
        axis (str): Indicator as to which axis to place a label. Allows 'x' or 'y' as parameters
        label (str): Label to place on the selected axis
        label_prop: Keyword dictionary to be set to plt.text
        labelpad (float): Padding from the axis (default: 0)
        ha (str): horizontal alignment (default: 'center')
        va (str): vertical alignment (default: 'center')

    .. codeauthor:: pygeostat development team - 2015-09-30
    """
    if fig is None:
        fig = plt.gcf()
    xmin = []
    ymin = []
    for ax in fig.axes:
        xmin.append(ax.get_position().xmin)
        ymin.append(ax.get_position().ymin)
    xmin, ymin = min(xmin), min(ymin)
    dpi = fig.dpi
    if axis.lower() == "y":
        if rotation is None:
            rotation = 90
        x = xmin - float(labelpad) / dpi
        y = 0.5
    elif axis.lower() == 'x':
        if rotation is None:
            rotation = 0
        x = 0.5
        y = ymin - float(labelpad) / dpi
    else:
        raise Exception("Unexpected axis: x or y")
    if label_prop is None:
        label_prop = dict()
    plt.text(x, y, label, rotation=rotation, transform=fig.transFigure, ha=ha, va=va, **label_prop)


def getminmax(vlim, sigfigs):
    """
    A small utility function called from many of the plotting functions. This determines some
    logical minimum and maximum values to bound a dataset based on a specified number of
    significant figures.

    Parameters:
        vlim (float tuple): Minimum and maximum values to honor if  parameter of "None" is not
            passed (vmin, vmax)
        sigfigs (int): Number of significant figures to consider

    Returns:
        valrng (float tuple): Calculated logical (min, max) bounds

    .. codeauthor:: pygeostat development team 2015-10-13
    """
    # Set-up some parameters
    valrng = [copy.deepcopy(vlim[0]), copy.deepcopy(vlim[1])]
    valabs = [False, False]
    flip = False
    # Get absolute values if required
    if valrng[0] < 0:
        valabs[0] = True
        valrng[0] = abs(valrng[0])
    if valrng[1] < 0:
        valabs[1] = True
        valrng[1] = abs(valrng[1])
    # Flip the value range if required
    if valrng[0] > valrng[1]:
        flip = True
        valrng = valrng[::-1]
    # Get logical min/max values
    try:
        precision = int(('{:.%ie}' % sigfigs).format(valrng[1]).partition('+')[2])
    except ValueError:
        precision = -1 * int(('{:.%ie}' % sigfigs).format(valrng[1]).partition('-')[2])
    valrng[1] = valrng[1] / (10**(precision - sigfigs + 1))
    valrng[0] = valrng[0] / (10**(precision - sigfigs + 1))
    valrng[1] = valrng[1] + 0.5
    if valabs[0]:
        valrng[0] = valrng[0] + 0.5
    else:
        valrng[0] = valrng[0] - 0.5
    valrng[1] = round(valrng[1], 0) * (10.0**(precision - sigfigs + 1)) + 0
    valrng[0] = round(valrng[0], 0) * (10.0**(precision - sigfigs + 1)) + 0
    # Flip back if required
    if flip:
        valrng = valrng[::-1]
    # Undo absolute values if required
    if valabs[0]:
        valrng[0] = valrng[0] * -1 + 0
    if valabs[1]:
        valrng[1] = valrng[1] * -1 + 0

    return (valrng[0], valrng[1])


def get_contcbarargs(data, sigfigs, vlim=None, nticks=5, catdata=None):
    """
    A small utility function called from many of the plotting functions. This determines logical
    continuous colorbar arguments if value limit values are not provided.
    If they are, determine some intermediate values.

    Parameters:
        data: A numpy array containing the data that is being plotted and required a colorbar
        sigfigs (int): The number of significant figures to consider
        vlim (float tuple): Value limits of the data being plotted as a tuple (vmin, vmax)
        nticks (int): The number of tick labels that need to be generated
        catdata (bool): Basically assume the codes are integers

    Returns
    -------
        vlim: float tuple
            A tuple of the new limits to use for plotting
        ticklocs: list
            A list of tick locations for the colorbar
        ticklabels: list
            A list of tick labels for each tickloc value

    .. codeauthor:: pygeostat development team 2015-09-30
    """
    vmin = np.nanmin(np.array(data))
    vmax = np.nanmax(np.array(data))
    # If limits weren't provided, get them from the data and set-up a boolean for later
    if vlim is None:
        vlim = (vmin, vmax)
        # Move the vmin and vmax values to something logical
        vlim = getminmax(vlim, sigfigs)
        vlims = False
    else:
        tvlim = []
        if hasattr(vlim[0], 'values'):
            tvlim.append(vlim[0].values[0])
        else:
            tvlim.append(vmin if vlim[0] is None else vlim[0])
        if hasattr(vlim[1], 'values'):
            tvlim.append(vlim[1].values[0])
        else:
            tvlim.append(vmax if vlim[1] is None else vlim[1])
        vlim = tuple(tvlim)
        vlims = True
    # Determine the other tick locations from highest to lowest
    step = float(1 / (nticks - 1))  # step percentage of each tick
    diff = abs(vlim[1] - vlim[0])  # difference between the limits
    if vlims:
        ticklocs = [vlim[1]]
    else:
        ticklocs = [float(('{:.%ig}' % sigfigs).format(vlim[1]))]
    for i in range(nticks - 2, 0, -1):
        value = diff * (i * step) + vlim[0]
        value = float(('{:.%ig}' % sigfigs).format(value))
        ticklocs.append(value)
    if vlims:
        ticklocs.append(vlim[0])
    else:
        ticklocs.append(float(('{:.%ig}' % sigfigs).format(vlim[0])))
    # Check to make sure that the data doesn't contain any values beyond the specified limits, if
    # so, add a greater than or less than sign to the respective label.
    if catdata:
        ticklabels = ["%i" % v for v in np.unique(data)]
    else:
        ticklabels = ticklocs
        # Convert the labels to strings
        ticklabels = [('{:.%ig}' % sigfigs).format(label) for label in ticklabels]
    # Add on the greather than or less than sign if there are values within the dataset that are
    # outside of the provided limits
    if vlim[0] > vmin:
        ticklabels[nticks - 1] = r'$\leq$ %s' % round(vlim[0], sigfigs)
    if vlim[1] < vmax:
        ticklabels[0] = r'$\geq$ %s' % round(vlim[1], sigfigs)

    return vlim, ticklocs, ticklabels


def get_palette(palette, ncat, cmap=True):
    """
    Return a colour palette with the required number of colours. The available palettes are
    ``cat_pastel``, ``cat_vibrant``, and ``cat_dark``. Please refer to the colour palette
    documentation for more information.

    The returned data can be either a list of colours or a Matplotlib ListColormap object.

    Parameters:
        palette (str): Desired colour palette
        ncat (int): The number of categories being plotted
        cmap (bool): Indicate if the retruned pallet is a list of colors or a Matplotlib
            ListColormap

    Returns
    -------
        palette: Matplotlib colormap, list
            if ``cmap`` is ``True`` or a list if ``cmap`` is ``False`` with the required number of
            colours

    Examples:
        To get a list of colour palettes within pygeostat in python, call the following function:

        >>> gs.cmaps.avail_palettes

    .. codeauthor:: pygeostat development team 2015-10-08
    """
    from . cmaps import cat_palettes

    # Check to make sure the number of colours needed doesn't exceed the number the palette selected
    # can provide.
    if palette in ['cat_pastel', 'cat_vibrant'] and ncat > 12:
        print("The colour palette '%s' can only accomodate up to 12 catagories" % palette)
        return
    if palette in ['cat_dark'] and ncat > 6:
        print("The colour palette '%s' can only accomodate up to 6 catagories" % palette)
        return
    # Get the required colour palette
    colours = cat_palettes[palette][:ncat]
    if cmap:
        palette = ListedColormap(colours, name=palette)
    else:
        palette = colours

    return palette


def get_cmap(colormap):
    """
    Return a colormap that is saved within pygeostat. Available colormaps are ``topo1``,
    and ``topo1``. These colormaps are used for hillshading.

    Parameters:
        colormap (str): Desired colormap

    Returns:
        cmap (cmap): Matplotlib colormap

    Examples:
        To get a list of colormaps within pygeostat in python, call the following function:

        >>> gs.get_cmap('topo1')

    .. codeauthor:: pygeostat development team - 2015-10-08
    """
    from .cmaps import cont_cmaps

    # Grab the desired colormap
    cmap = ListedColormap(cont_cmaps[colormap], name=colormap)

    return cmap


def catcmapfromcontinuous(contcmap, ncats, offset=0.5):
    """
    Create a categorical colormap from a continuous colormap.
    
    Takes a continuous colormap and samples it at regular intervals to create 
    a categorical colormap with the specified number of categories.
    
    Parameters
    ----------
    contcmap : str or matplotlib.colors.Colormap
        Valid matplotlib continuous colormap string or colormap object
    ncats : int
        Number of categories to create in the resulting colormap
    offset : float, default=0.5
        Offset value to avoid using extreme end colors of the colormap
        
    Returns
    -------
    matplotlib.colors.ListedColormap
        A categorical colormap with ncats colors
        
    Notes
    -----
    The offset parameter (default 0.5) helps avoid the extreme colors at the 
    ends of the colormap, which are often too light or too dark.
    """
    from . cmaps import avail_cmaps
    from matplotlib.colors import ListedColormap
    import matplotlib as mpl

    if contcmap in avail_cmaps:
        # Use custom pygeostat colormap
        cm = get_cmap(contcmap)
    else:
        try:
            cm = mpl.colormaps[contcmap]
        except (KeyError, TypeError):
            raise ValueError(f"Colormap '{contcmap}' not found in matplotlib or pygeostat colormaps")
    
    colors = []
    for num in range(ncats):
        position = (num + offset) / ncats
        colors.append(cm(position))

    return ListedColormap(colors, name=f"{contcmap}_categorical_{ncats}")







def get_label(data):
    """Small utility to grab label information from different types of data sources

    Parameters:
        data: Tidy (long-form) dataframe where each column is a variable and each row is an
            observation. Pandas dataframe or numpy array

    Returns:
        label (str): the name of the column

    .. codeauthor:: pygeostat development team
    """
    label = None
    if hasattr(data, 'columns'):
        if len(data.columns) == 1:
            label = data.columns[0]
    if hasattr(data, 'name'):
        label = data.name

    return label


def clrmplmem():
    """
    Clear the current plot from memory.

    When looping the generation of plots, matplotlib does not remove data from memory until it is
    needed, even if the namespace used is repeated. Instead of clearing the data, the namespace is
    simply remaped. This function uses a series of commands that will dump physical memory.

    Example:

    >>> for var in variables:
    >>>     gs.slice_plot(data[var], griddef=griddef)
    >>>     gs.clrmplmem()

    .. codeauthor:: pygeostat development team - 2016-03-22
    """
    fig = plt.gcf()
    fig.clf()
    plt.close()
    gc.collect()


def addxticks(fig, nrow, ncol, lastax, axextents=None):
    """
    A small utility function called from many of the plotting functions. This adds xticklabels to
    the last subplot in a column when it is not in the last row.

    .. codeauthor:: Ryan Martin and pygeostat development team - 2016-04-18
    """
    import pygeostat as gs

    #  Add xticklabels to the last subplot in a column when it is not in the last row
    renderer = fig.canvas.get_renderer()
    fig.draw(renderer)
    xaxis = fig.axes[ncol * (nrow - 1)].get_xaxis()
    xticklocs = []
    xticklabels = []
    for text in xaxis.get_ticklabels():
        xticklocs.append(text.get_position()[0])
        xticklabels.append(text.get_text())
    xticklocs = np.array(xticklocs)
    if axextents is not None:
        inds = (xticklocs >= axextents[0]) & (xticklocs <= axextents[1])
        xticklabels = np.array(xticklabels)
        xticklocs = xticklocs[inds]
        xticklabels = xticklabels[inds]
    for i in range(lastax, nrow * ncol):  # iax of axes that need xticklabels
        iax = i - ncol
        xaxis = fig.axes[iax].get_xaxis()
        xaxis.tick_bottom()
        xaxis.set_ticks(xticklocs)
        xaxis.set_ticklabels(xticklabels)


def get_supaxislocs(fig, nrow, ncol, figsize, pad):
    """
    A small utility function called from many of the plotting functions. This gets the required
    coordinates for super axis labels.

    Parameters:
        fig (mpl.fig): Matplotlib figure
        nrow (int): Number of rows
        ncol (int): Number of Columns
        figsize (tuple): Figure size (width, height)
        pad (float): Separation between items and the output coordinates

    Returns
    -------
        xmin: float
            Minimum coordinate of figure in x-axis
        xmid: float
            Center coordinate of figure in x-axis
        ymin: float
            Minimum coordinate of figure in y-axis
        ymid: float
            Center coordinate of figure in y-axis
        ymax: float
            Maximum coordinate of figure in y-axis

    .. codeauthor:: pygeostat development team 2016-04-18
    """
    import pygeostat as gs

    renderer = fig.canvas.get_renderer()
    fig.draw(renderer)
    dpi = fig.dpi
    axes = fig.axes
    #  Get the middle coordinate in the x-axis of the grid
    if ncol % 2 == 0:
        #  Get the max x value from the mid left axis
        iax = int(ncol / 2 - 1)
        bbox = axes[iax].get_tightbbox(renderer)
        xmax = bbox.xmax / figsize[0] / dpi
        #  Get the min x value from the mid right axis
        bbox = axes[iax + 1].get_tightbbox(renderer)
        xmin = bbox.xmin / figsize[0] / dpi
        #  Calculate the xmid coordinate
        xmid = ((xmax - xmin) / 2 + xmin)
    else:
        #  Get the middle coordinate from the middle axis
        iax = int((ncol / 2) - 0.5)
        bbox = axes[iax].get_tightbbox(renderer)
        xmin = bbox.xmin / figsize[0] / dpi
        xmax = bbox.xmax / figsize[0] / dpi
        #  Calculate the xmid coordinate
        xmid = ((xmax - xmin) / 2 + xmin)
    #  Get the middle coordinate in the y-axis of the grid
    if nrow % 2 == 0:
        #  Get the max y coord from the middle bottom axis
        iax = int((nrow / 2)) * ncol
        bbox = axes[iax].get_tightbbox(renderer)
        ymax = bbox.ymax / figsize[1] / dpi
        #  Get the min y value from the middle top axis
        iax = (int((nrow / 2)) - 1) * ncol
        bbox = axes[iax].get_tightbbox(renderer)
        ymin = bbox.ymin / figsize[1] / dpi
        #  Calculate the ymid coordinate
        ymid = ((ymax - ymin) / 2 + ymin)
    else:
        #  Get the middle coordinate from the middle axis
        iax = int((nrow / 2) - 0.5) * ncol
        bbox = axes[iax].get_tightbbox(renderer)
        ymin = bbox.ymin / figsize[1] / dpi
        ymax = bbox.ymax / figsize[1] / dpi
        #  Calculate the ymid coordinate
        ymid = ((ymax - ymin) / 2 + ymin)
    #  Get the top of the axes for the suptitle:
    bbox = axes[0].get_tightbbox(renderer)
    ymax = bbox.ymax / figsize[1] / dpi
    ymax = ymax + pad
    #  Get the limit of the xticklabels
    iax = ncol * (nrow - 1)
    #bbox = axes[iax].xaxis.get_ticklabel_extents(renderer)[0]
    #ymin = (bbox.ymin / figsize[1] / dpi) - pad
    ##  Get the x limit of the yticklabels
    #bbox = axes[0].yaxis.get_ticklabel_extents(renderer)[0]
    #xmin = (bbox.xmin / figsize[0] / dpi) - pad
    
    # New implementation due to deprecation in matplotlib
    bbox = axes[iax].get_tightbbox(renderer)
    ymin = (bbox.ymin / figsize[1] / dpi) - pad
    
    # For left margin
    bbox = axes[0].get_tightbbox(renderer)
    xmin = (bbox.xmin / figsize[0] / dpi) - pad

    return xmin, xmid, ymin, ymid, ymax

def _get_cmap(cmap, catdata, ncat):
    '''
    Determine the colormap based on categorical data flag, provided colormap 
    name, and parameters.
    Used by location_plot and other plotting functions.
    
    Parameters
    ----------
    cmap : str, matplotlib.colors.Colormap, or None
        Name of colormap or colormap object
    catdata : bool or None
        Flag indicating if data is categorical
    ncat : int
        Number of categories (used for categorical data)
        
    Returns
    -------
    matplotlib.colors.Colormap
        The selected colormap
    '''
    from . cmaps import avail_palettes, avail_cmaps
    from matplotlib.colors import LinearSegmentedColormap
    
    if catdata:
        # Get default categorical colormap if none provided
        if cmap is None:
            cmap = Parameters['plotting.cmap_cat']
            if isinstance(cmap, dict):
                cmap = list(cmap.values())
                cmap = LinearSegmentedColormap.from_list('cmap_cat', cmap, N=len(cmap))
        
        if not isinstance(cmap, mpl.colors.Colormap):
            if cmap in avail_palettes:
                # Use pygeostat's custom palette
                cmap = get_palette(cmap, ncat)
            else:
                cmap = catcmapfromcontinuous(cmap, ncat)
    
    # Handle continuous data case
    else:
        if cmap is None:
            cmap = Parameters['plotting.cmap']
            
        if not isinstance(cmap, mpl.colors.Colormap):
            if cmap in avail_cmaps:
                # Use pygeostat's custom colormap
                cmap = get_cmap(cmap)
            else:
                try:
                    cmap = mpl.colormaps[cmap]
                except (KeyError, TypeError):
                    pass
    
    return cmap







def setup_plot(ax, cbar=None, figsize=None, cax=None, aspect=None):
    '''A small utility function called from many of the plotting functions. This will set up a
    matplotlib plot instance based on whether an axis is passed or not.

    Parameters:
        ax (mpl.axis): Matplotlib axis to plot the figure
        cbar (bool): Indicate if a colorbar should be plotted or not
        figsize (tuple): Figure size (width, height)
        cax: Matplotlib.ImageGrid.cbar_axes object
        aspect (bool, str, float): Bool for creating axes, str or float
            for existing axes

    Returns:
        fig (mpl.plt.fig): Matplotlib figure
        ax (mpl.axis): Matplotlib axis to plot the figure
        cax: Matplotlib.ImageGrid.cbar_axes object
    '''
    from mpl_toolkits.axes_grid1 import ImageGrid

    if ax is None:
        # Setup up a new plot
        fig = plt.figure(figsize=figsize)
        cbar_mode = None
        if cax is None:
            if cbar:
                cbar_mode = 'single'
        if aspect is None:
            aspect = True
        imggrid = ImageGrid(fig, 111, nrows_ncols=(1, 1), axes_pad=0.07,
                            cbar_mode=cbar_mode, cbar_size=0.075, aspect=aspect)
        ax = imggrid[0]
        if cax is None:
            cax = imggrid.cbar_axes[0]
    elif hasattr(ax, "cax"):
        cax = ax.cax
        fig = plt.gcf()
    elif cbar:
        try:
            fig, cax = get_cbar_axis(ax, cax)
        except:
            fig = plt.gcf()
            if hasattr(ax, 'cax'):
                cax = ax.cax
            if cax is None:
                raise ValueError("A colorbar axes `cax` must be passed as the passed `ax` cannot be"
                                " divided.")
    else:
        fig = plt.gcf()
    return fig, ax, cax

_custom_cbar_width = 0.008
_custom_cbar_pad = 0.008


def _get_mpl_cbar_callback(ax):
    def cax_pos_callback(event):
        ax_rect = ax.get_position()
        ax.cax.set_position([ax_rect.x1 + _custom_cbar_pad,
                             ax_rect.y0, _custom_cbar_width,
                             ax_rect.height])
    return cax_pos_callback


def get_cbar_axis(ax, cax=None):
    """Add a cax (replaces make_axes_locatable)

    related to: https://github.com/matplotlib/matplotlib/issues/4282

    .. codeauthor:: RMS, 2019
    """
    fig = plt.gcf()
    if cax is None:
        cax = getattr(ax, 'cax', None)
    if cax is None:
        ax_rect = ax.get_position()
        cax = fig.add_axes([ax_rect.x1 + _custom_cbar_pad,
                            ax_rect.y0, _custom_cbar_width,
                            ax_rect.height])
        ax.cax = cax
        fig.canvas.mpl_connect('draw_event', _get_mpl_cbar_callback(ax))
    return fig, cax
# end RMS
# ----------------------------------------------------------------------------------------------


def format_subplot_axis(fig, ax, plotformat_dict, cbar=None):
    """
    A small utility function that will allow you to pass a dictionary for formatting such things as
    the fontsize and color of subplot axis

    Parameters:
        fig (matplotlib.fig): the figure for the subplots
        ax (matplotlib.ax): the matplotlib axis for the subplots
        cbar (matplotlib.fig): the cbar figure used in pygeostat plotting functions
        plotformat_dict (dict): dictionary used to format the axis

    Returns:
        fig, ax
        fig, ax, cbar

    Dictionary Keys:

        'fontsize': Dictionary

            * 'title' (int): change title fontsize
            * 'cbar_title' (int): change cbar title fontsize
            * 'xaxis_lable' (int): change xaxis label fontsize
            * 'xaxis_ticklabels' (int): change xaxis tick labels fontsize
            * 'yaxis_label' (int): change yaxis label fontsize
            * 'yaxis_ticklabels' (int): change taxis tick labels fontsize
            * 'cbar_ticklabels' (int): change cbar tick labels fontsize

        'color': Dictionary

            * 'title' (str): change title fontsize
            * 'cbar_title' (str): change cbar title fontsize
            * 'xaxis_lable' (str): change xaxis label fontsize
            * 'xaxis_ticklabels' (str): change xaxis tick labels fontsize
            * 'yaxis_label' (str): change yaxis label fontsize
            * 'yaxis_ticklabels' (str): change taxis tick labels fontsize
            * 'cbar_ticklabels' (str): change cbar tick labels fontsize

    .. codeauthor:: Tyler Acorn - Jan 03, 2018
    """
    if 'fontsize' in plotformat_dict:
        if isinstance(plotformat_dict['fontsize'], dict):
            fs_dict = plotformat_dict['fontsize']
            if 'title' in fs_dict:
                ax.title.set_fontsize(fs_dict['title'])
            if cbar and 'cbar_title' in fs_dict:
                cbar.ax.yaxis.label.set_fontsize(fs_dict['cbar_title'])
            if 'xaxis_label' in fs_dict:
                ax.xaxis.label.set_fontsize(fs_dict['xaxis_label'])
            if 'xaxis_ticklabels' in fs_dict:
                ax.tick_params(axis='x', labelsize=fs_dict['xaxis_ticklabels'])
            if 'yaxis_label' in fs_dict:
                ax.yaxis.label.set_fontsize(fs_dict['yaxis_label'])
            if 'yaxis_ticklabels' in fs_dict:
                ax.tick_params(axis='y', labelsize=fs_dict['yaxis_ticklabels'])
            if cbar and 'cbar_ticklabels' in fs_dict:
                cbar.ax.tick_params(axis='y', labelsize=fs_dict['yaxis_ticklabels'])
    if 'color' in plotformat_dict:
        if isinstance(plotformat_dict['color'], dict):
            c_dict = plotformat_dict['color']
            if 'title' in c_dict:
                ax.title.set_color(c_dict['title'])
            if cbar and 'cbar_title' in fs_dict:
                cbar.ax.yaxis.label.set_color(c_dict['cbar_title'])
            if 'xaxis_label' in c_dict:
                ax.xaxis.label.set_color(c_dict['xaxis_label'])
            if 'xaxis_ticklabels' in c_dict:
                ax.tick_params(axis='x', color=c_dict['xaxis_ticklabels'])
            if 'yaxis_label' in c_dict:
                ax.yaxis.label.set_color(c_dict['yaxis_label'])
            if 'yaxis_ticklabels' in c_dict:
                ax.tick_params(axis='y', color=c_dict['yaxis_ticklabels'])
            if cbar and 'cbar_ticklabels' in c_dict:
                cbar.ax.tick_params(axis='y', color=c_dict['yaxis_ticklabels'])

    if cbar:
        return fig, ax, cbar
    else:
        return fig, ax


def scalebar(x, y, scale, length, ax, img, height=0.005, units='km', pad=0.02):
    """
    A small utility function called from many of the plotting functions. This will Add a scale bar
    to a plot. Requires an Image object such as imshow.

    Parameters:
        x (float): X coordinate of the scale bar in figure space
        y (float): Y coordinate of the scale bar in figure space
        scale (float): Side of each pixel in the units specified by the parameter ``units``
        length (float): Length of the scale bar in the units specified by the parameter ``units``
        ax (mpl.ax): Matplotlib axis to plot the scale bar onto
        img (mpl.AxesImage): Matplotlib AxesImage object. Produced by from imshow()
        height (float): Height of the scale bar as a percentage of plot in y direction in figure
            space
        units (str): Units of the scale bar
        pad (float): Padding between the scale bar and the annotation

    .. codeauthor:: pygeostat development team - 2016-06-13
    """
    # Get some extents
    extents = img.get_extent()
    xsiz = (abs(extents[1]) + abs(extents[0]))
    ysiz = (abs(extents[2]) + abs(extents[3]))
    xmax = max(extents[1], extents[0])
    ymax = max(extents[2], extents[3])
    xmid = xmax - (xsiz * (1 - x))
    pad = ysiz * pad
    height = ysiz * height
    if extents[2] > extents[3]:
        ymid = ymax - (ysiz * y)
        ymin = ymid + (height / 2)
        txty = ymin + pad
    else:
        ymid = ymax - (ysiz * (1 - y))
        ymin = ymid - (height / 2)
        txty = ymin - pad
    length_px = length / scale
    xmin = xmid - (length_px / 2)
    # Make the scale bar rectangle
    rectangle = mpl.patches.Rectangle((xmin, ymin), length_px, height)
    rectangle.set(edgecolor='k', facecolor='k')
    ax.add_patch(rectangle)
    # Add tannotatio
    ax.text(xmid, txty, '%s %s' % (length, units), ha='center', va='top')


def smart_annotate(ax, x, y, labels, k=0.002):
    """
    A small utility function that can be called from a plotting funtion. This annotates a scatter
    plot and moves the annotations around so they do not overlap. Works okay...

    Will need to be made more flexible in the future. Works for gs.mdsplt() at the moment.

    Retrieved and modified from: http://stackoverflow.com/a/34697108

    Parameters:
        ax (mpl.axis): Matplotlib axis object the scatter plot is plotted on
        x: 1-D array with the x-coordinates used for the data in the scatter plot
        y: 1-D array with the y-coordinates used from the data in the scatter plot
        labels (list): Labels to annotate the scatter plot with
        k (float): Scalar for the distance from the data the annotation is placed

    .. codeauthor:: pygeostat development team - 2016-05-30
    """
    import networkx as nx
    G = nx.DiGraph()
    data_nodes = []
    init_pos = {}
    for xi, yi, label in zip(x, y, labels):
        data_str = 'data_{0}'.format(label)
        G.add_node(data_str)
        G.add_node(label)
        G.add_edge(label, data_str)
        data_nodes.append(data_str)
        init_pos[data_str] = (xi, yi)
        init_pos[label] = (xi, yi)
    pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=k)
    # undo spring_layout's rescaling
    pos_after = np.vstack([pos[d] for d in data_nodes])
    pos_before = np.vstack([init_pos[d] for d in data_nodes])
    scale, shift_x = np.polyfit(pos_after[:, 0], pos_before[:, 0], 1)
    scale, shift_y = np.polyfit(pos_after[:, 1], pos_before[:, 1], 1)
    shift = np.array([shift_x, shift_y])
    for key, val in pos.items():
        pos[key] = (val * scale) + shift

    for label, data_str in G.edges():
        ax.annotate(label,
                    xy=pos[data_str], xycoords='data',
                    xytext=pos[label], textcoords='data',
                    arrowprops=dict(arrowstyle="-",
                                    shrinkA=0, shrinkB=0,
                                    connectionstyle="arc3",
                                    color='k'))
    # expand limits
    all_pos = np.vstack(pos.values())
    x_span, y_span = np.ptp(all_pos, axis=0)
    mins = np.min(all_pos - x_span * 0.02, 0)
    maxs = np.max(all_pos + y_span * 0.02, 0)
    ax.set_xlim([mins[0], maxs[0]])
    ax.set_ylim([mins[1], maxs[1]])


def get_statblk(stat_blk, statsets, statlist, stat_xy):
    """
    A small utility function called from many of the plotting functions. Extracts the required
    statistics and return a string to plot and alignment settings.

    If automatic alignment is not desired, pass the string ``'noalign'`` through ``stat_blk``.

    Note done

    Returns:
        txtstats, stat_xy, ha, va

    .. codeauthor:: pygeostat development team - 2016-07-22
    """
    # Handle input
    if isinstance(stat_blk, str):
        stat_blk = [stat_blk]
    presets = list(statsets.keys())
    stats = list(statlist.keys())
    if 'noalign' in stat_blk:
        align = False
        stat_blk.remove('noalign')
    elif stat_xy in ['bottom', 'top']:
        align = False
    else:
        align = True
    # Sanity checks
    if len(stat_blk) == 0:
        raise ValueError("No parameter was passed to `stat_blk`")
    for stat in stat_blk:
        if stat not in (presets + stats):
            raise ValueError("The following value is not a valid `stat_blk` parameter: %s\n "
                             "Try one of: %s" % (stat, (presets + stats)))
    # Get the right statistics
    if stat_blk[0] in presets:
        stat_blk = [statlist[stat] for stat in statsets[stat_blk[0]]]
    else:
        stat_blk = [statlist[stat] for stat in stat_blk]
    txtstats = '\n'.join(stat_blk)
    # Set-up the alignment
    if align:
        if stat_xy[0] > 0.5:
            ha = 'right'
        else:
            ha = 'left'
        if stat_xy[1] > 0.5:
            va = 'top'
        else:
            va = 'bottom'
    elif stat_xy == 'bottom':
        va = 'bottom'
        ha = 'left'
        stat_xy = (1, 0.05)
    elif stat_xy == 'top':
        va = 'top'
        ha = 'left'
        stat_xy = (1, 0.95)
    else:
        ha = 'left'
        va = 'top'

    return (txtstats, stat_xy, ha, va)


def _set_stat_fontsize(stat_fontsize):
    '''
    Determine the stat_fontsize based on the provided parameter, Parameters
    '''
    if stat_fontsize is None:
        stat_fontsize = Parameters['plotting.stat_fontsize']
    if stat_fontsize is None:
        stat_fontsize = mpl.rcParams['font.size']
    else:
        try:
            if stat_fontsize <= 0:
                raise ValueError('stat_fontsize must be greater than 0!')
            elif stat_fontsize < 1:
                stat_fontsize = mpl.rcParams['font.size'] * stat_fontsize
            # Otherwise, stat_fontsize remains the stat_fontsize
        except:
            raise ValueError('stat_fontsize must be a float or integer!')
    return stat_fontsize


def applytickbins(ax, ntickbins=None):
    """
    Quick function to deal with tuple or integer input

    Parameters:
        ax (mpl.axes): the axis to apply the tickbins too
        ntickbins (int, tuple): either an int (applied to both) or tuple for x and y respectively
    """
    if ntickbins is not None:
        xbins, ybins = None, None
        try:
            iter(ntickbins)
            if len(ntickbins) == 2:
                xbins, ybins = ntickbins
            else:
                xbins, ybins = 5, 5
        except Exception:
            pass
        if isinstance(ntickbins, int):
            xbins, ybins = ntickbins, ntickbins
        if xbins is not None:
            ax.locator_params(axis='x', nbins=xbins)
        if ybins is not None:
            ax.locator_params(axis='y', nbins=ybins)


def format_plot(ax, xlabel=None, ylabel=None, title=None, grid=None, axis_xy=None,
              xlim=None, ylim=None, logx=None, logy=None):
    '''
    Format a plot with common properties such as labels and grid/axis visibility. Not used by the
    spatial plotting function.

    Parameters:
        ax(Matplotlib axis handle): axis handle to modify

    Keyword arguments:
        xlabel(str): label of the x-axis
        ylabel(str): label of the y-axis
        title(str): plot title
        grid(bool): plots the major grid lines if True. Based on Parameters['plotting.grid']
            if None.
        axis_xy(bool): converts the axis to GSLIB-style axis visibility (only left and bottom
            visible) if axis_xy is True. Based on Parameters['plotting.axis_xy'] if None.
        xlim(tuple): set x-axis limits to xlim[0] and xlim[1]
        ylim(tuple): set y-axis limits to ylim[0] and ylim[1]
        logx, logy (bool): if true, convert each axis to log-scale

    Return:
        ax(Matplotlib axis handle): modified axis handle
    '''
    from matplotlib.ticker import StrMethodFormatter
    # Format the axis labels
    if xlabel is not None:
        if xlabel is not False:
            ax.set_xlabel(xlabel, fontsize=None)
    if ylabel is not None:
        if ylabel is not False:
            ax.set_ylabel(ylabel, fontsize=None)
    if title is not None:
        if title is not False:
            ax.set_title(title, fontsize=None)
    # Format the axis splines
    _format_axis_xy(ax, axis_xy)
    # Format the grid lines
    _format_grid(ax, grid)
    # deal with log axes and ticks
    sigfigs = Parameters.get('plotting.sigfigs')
    if sigfigs is None:
        sigfigs = 3
    fmtstr = '{{x:.{}g}}'.format(sigfigs)
    if logx:
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(StrMethodFormatter(fmtstr))
    if logy:
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(StrMethodFormatter(fmtstr))
    # Set the x and y limits
    if xlim is not None:
        if logx and xlim[0] <= 0:
            xlower = 1e-6
        else:
            xlower = xlim[0]
        ax.set_xlim(xlower, xlim[1])
    if ylim is not None:
        if logy and ylim[0] <= 0:
            ylower = 1e-6
        else:
            ylower = ylim[0]
        ax.set_ylim(ylower, ylim[1])

    return ax


def _format_axis_xy(ax, axis_xy=None):
    '''
    Private function for the handling axis_xy, which converts the axis to
        GSLIB-style axis visibility (only left and bottom visible) if
        axis_xy is True.

    Parameters:
        ax(matplotlib axis handle): axis to modify
        axis_xy(bool): convert to axis_xy if True

    Return:
        ax(Matplotlib axis handle): modified axis handle
    '''
    if axis_xy is None:
        axis_xy = Parameters['plotting.axis_xy']
    if axis_xy:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)


def _format_grid(ax, grid=None, below=False):
    '''
    Private function for the handling of grid lines. Matplotlib kwargs
        could be added in the future for permitting custom grid lines,
        as was done in the previous pygeostat version (commented out below)

    Parameters:
        ax(matplotlib axis handle): axis to append grid lines to (or not)
        grid(bool): add grid lines, or not...
        below(bool): place the grid lines below (Ture) or above (False)
            the plot content

    Return:
        ax(Matplotlib axis handle): modified axis handle
    '''
    if grid is None:
        grid = Parameters['plotting.grid']
    if grid:
        # The grid lines won't be seen in many plots if not above...
        # Could consider having this as an option
        ax.set_axisbelow(below)
        # The linestyle, width and color are available in rcParams and should, perhaps,
        # not be mangled with. That said, the pygeostat setting below looks way better than
        # the matplotlib default so leaving for now
        ax.grid(True, which='major', color='k', linestyle=':', lw=0.2)


def _spatial_labels(ax, orient, griddef=None, slice_number=0, title=None, xlabel=None, ylabel=None,
                    unit=None, sigfigs=3):
    '''
    A small utility function called from plotting functions that involve maps and cross-sections.
    It provides the appropriate labels and title according to a variety of parameter options.

    Parameters:
        ax (mpl.axis): Matplotlib axis to plot the figure
        orient (str): Orientation to slice data. ``'xy'``, ``'xz'``, ``'yz'`` are the only accepted
            values
        griddef (GridDef): A pygeostat GridDef class created using
            :class:`gs.GridDef <pygeostat.data.grid_definition.GridDef>`
        slice_number (int): Grid cell location along the axis not plotted to take the slice of data to
            plot
        title (str): Title for the plot. If left to it's default value of ``None`` or is set to
            ``True``, a logical default title will be generated for 3-D data. Set to ``False`` if
            no title is desired.
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        unit (str): Unit to place inside the axis label parentheses

    Returns:
        ax (mpl.axis): Matplotlib axis to plot the figure
    '''
    from ..data.grid_definition import GridDef
    if unit is None:
        unit = Parameters['plotting.unit']
    if unit is None or unit == '':
        unit = ''
    else:
        unit = ' ({})'.format(unit)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=None)
    elif xlabel is None:
        if orient == 'xy':
            ax.set_xlabel(Parameters['plotting.xname'] + unit, fontsize=None)
        elif orient == 'xz':
            ax.set_xlabel(Parameters['plotting.xname'] + unit, fontsize=None)
        elif orient == 'yz':
            ax.set_xlabel(Parameters['plotting.yname'] + unit, fontsize=None)
    if ylabel:
        ax.set_ylabel(ylabel)
    elif ylabel is None:
        if orient == 'xy':
            ax.set_ylabel(Parameters['plotting.yname'] + unit, fontsize=None)
        elif orient == 'xz':
            ax.set_ylabel(Parameters['plotting.zname'] + unit, fontsize=None)
        elif orient == 'yz':
            ax.set_ylabel(Parameters['plotting.zname'] + unit, fontsize=None)
    if unit != '':
        unit = unit[2:-1]
    if isinstance(title, str):
        ax.set_title(title, fontsize=None)
    elif (title is None or title is True) and isinstance(griddef, GridDef):
        crossd = list(set(['x', 'y', 'z']) - set(list(orient)))[0]
        d3 = False
        if crossd == 'x' and griddef.nx > 1:
            d3 = True
        if crossd == 'y' and griddef.ny > 1:
            d3 = True
        if crossd == 'z' and griddef.nz > 1:
            d3 = True
        if d3:
            if orient in ['xz', 'yz']:
                cross = griddef.get_slice_coordinate(orient, slice_number)
                if crossd == 'y':
                    crossd = Parameters['plotting.yabbrev']
                else:
                    crossd = Parameters['plotting.xabbrev']
                title = 'Cross Section: %s %s %s' % (round(cross, sigfigs), unit, crossd)
            else:
                plan = griddef.get_slice_coordinate(orient, slice_number)
                title = 'Plan View: %s %s' % (round(plan, sigfigs), unit)
        else:
            title = None
    else:
        title = None
    if title:
        ax.set_title(title, fontsize=None)
    return ax


def _spatial_pointdata(data, orient, x=None, y=None, z=None, griddef=None):
    '''Used by slice_plot and location_plot'''
    import pandas as pd
    from ..data.data import DataFile
    from ..data.grid_definition import GridDef
    import pygeostat as gs
    if not any([isinstance(data, pd.DataFrame), isinstance(data, DataFile)]):
        raise ValueError("The parameter data must be either a gs.DataFile or a pd.DataFrame")
    if isinstance(data, DataFile):
        if not griddef:
            if isinstance(data.griddef, GridDef):
                griddef = data.griddef
        if x is None:
            x = data.x
        if y is None:
            y = data.y
        if z is None:
            z = data.z
        data = data.data
    if orient == 'xy':
        if any([x is None, y is None]):
            raise ValueError("`orient` is set to `'xy'` yet the column ID for `x` and `y` has"
                             " not been set and could not be automatically retrieved.")
    elif orient == 'xz':
        if any([x is None, z is None]):
            raise ValueError("`orient` is set to `'xz'` yet the column ID for `x` and `z` has"
                             " not been set and could not be automatically retrieved.")
    elif orient == 'yz':
        if any([y is None, z is None]):
            raise ValueError("`orient` is set to `'yz'` yet the column ID for `y` and `z` has"
                             " not been set and could not be automatically retrieved.")
    else:
        raise ValueError("An `orient` setting of %s is not permissible" % orient)
    return data, x, y, z, griddef


def _spatial_griddata(data, var, griddef):
    '''
    Private function that parses the data, var and griddef, returning the data
    and griddef to apply. Used by slice_plot and contourplt.
    '''
    from .. data.data import DataFile

    # Handle DataFile input
    if isinstance(data, DataFile):
        if griddef is None:
            griddef = data.griddef
        data = data.data
    # Ensure that we have a griddef
    if griddef is None:
        raise Exception('Error: griddef must be provided if data is not a DataFile object\n'
                        '       with the associated grid definition (data.griddef)')
    # Handle panda input
    if hasattr(data, 'values'):
        if isinstance(var, str) and data.shape[1] > 1:
            data = data[var].values
        else:
            data = data.values
            if data.ndim > 1:
                if isinstance(var, int):
                    # If var is provided as a scalar, then assume it's the column number
                    data = data[:, var]
                else:
                    data = data[:, 0]
    if data.ndim > 1:
        if data.shape[1] > 1:
            raise ValueError('could not coerce data into a 1-D array!')
    return data, griddef


def _spatial_slice(data, var, x, y, z, griddef, orient, slice_number, slicetol, dhid=None):
    '''Used by slice_plot and location_plot'''
    from .. datautils.utils import slicescatter
    from .. data.grid_definition import GridDef
    if not isinstance(griddef, GridDef) and slicetol is not None:
        raise ValueError("`griddef` must be set if using `slicetol`")
    # If a slice_number is passed but no slice tolerence. set slicetol to griddef cell size
    try:
        if slice_number >= 0:
            if slicetol is None:
                if orient == 'xy':
                    if griddef.nz > 1:
                        slicetol = griddef.zsiz / 2
                    else:
                        slicetol = None
                if orient == 'xz':
                    if griddef.ny > 1:
                        slicetol = griddef.ysiz / 2
                    else:
                        slicetol = None
                    slicetol = griddef.ysiz / 2
                if orient == 'yz':
                    if griddef.nx > 1:
                        slicetol = griddef.xsiz / 2
                    else:
                        slicetol = None
    except:
        pass
    # if no slice_number is passed plot everything. i.e. slice_number = 1, slicetol = None
    if slice_number is None:
        slice_number = 1
        slicetol = None
    # If slicetol is negative set it to plot everything
    try:
        if slicetol < 0:
            slicetol = None
    except:
        pass
    # Only slice if needed
    if slicetol is None:
        pointview = data
    else:
        pointview = slicescatter(data, orient, slice_number, slicetol, griddef, x, y, z)
        if len(pointview) == 0:
            #raise ValueError("No data was found within the specified slice tolerance")
            if dhid is None:
                return None, None, None
            else:
                return None, None, None, None
                
    figx, figy, figz = _spatial_orient2fig(orient, x, y, z)
    if figz is not None:
        pointsort = np.argsort(pointview[figz].values)
        pointx = pointview[figx].values[pointsort]
        pointy = pointview[figy].values[pointsort]
        try:
            pointvar = pointview[var].values[pointsort]
        except:
            pointvar = None
        try:
            point_dh = pointview[dhid].values[pointsort]
        except:
            point_dh = None

    else:
        pointx = pointview[figx]
        pointy = pointview[figy]
        try:
            pointvar = pointview[var]
        except:
            pointvar = None
        
        try:
            point_dh = pointview[dhid]
        except:
            point_dh = None

    if dhid is None:
        return pointx, pointy, pointvar
    else:
        return pointx, pointy, pointvar, point_dh


def _spatial_orient2fig(orient, x, y, z):
    '''
    Used by a variety of functions. Given orient and the x, y, z cols,
    return figx and figy, which is used for indexing.
    '''
    if orient == 'xy':
        figx = x
        figy = y
        figz = z
    elif orient == 'xz':
        figx = x
        figy = z
        figz = y
    elif orient == 'yz':
        figx = y
        figy = z
        figz = x
    else:
        raise ValueError('orient must be xy, xz or yz!')
    return figx, figy, figz


def _spatial_aspect(xlim, ylim):
    '''
    Simple utility for determining what aspect should be used for spatial plots.
    '''
    xrang = xlim[1] - xlim[0]
    yrang = ylim[1] - ylim[0]
    aspect = 'equal'
    if xrang / yrang > 5:
        aspect = 'auto'
    elif yrang / xrang > 5:
        aspect = 'auto'
    return aspect


def _format_tick_labels(ax, rotateticks=None, nticks=None):
    '''
    A small utility function called from many of the plotting functions. This will format the
    tick labels and rotate them if needed.

    Note - currently mangles padding and alignment according to the original pygeostat settings.

    Parameters:
        ax (mpl.axis): Matplotlib axis to plot the figure
        rotateticks (bool or float tuple): Bool or float values to control axis label rotations.
            See above for more info (x, y).

    '''
    if rotateticks is None:
        rotateticks = Parameters['plotting.rotateticks']
    if rotateticks is None:
        rotateticks = (None, None)
    if nticks is None:
        nticks = Parameters['plotting.nticks']
    try:
        applytickbins(ax, nticks)
    except Exception:
        pass
    ax.tick_params(axis='x', pad=2)
    ax.tick_params(axis='y', pad=1)
    if rotateticks is None:
        # The plots tick labels will not be properly accessible until the figure is "drawn", once
        # the command below is run, ax.get_ticklabel() will actually work properly.
        plt.draw()
        rotateticks = _tickoverlap(ax)
    if isinstance(rotateticks, tuple):
        rotateticks = list(rotateticks)
    # Rotate if required
    rotatex = False
    if isinstance(rotateticks[0], bool) and rotateticks[0]:
        rotateticks[0] = 45
        rotatex = True
    elif rotateticks[0] is not None:
        rotatex = True
    if rotatex:
        xlabels = ax.get_xticklabels()
        for xlabel in xlabels:
            xlabel.set_ha('center')
            xlabel.set_va('top')
            xlabel.set_rotation(rotateticks[0])
    rotatey = False
    if isinstance(rotateticks[1], bool) and rotateticks[1]:
        rotateticks[1] = -45
        rotatey = True
    elif rotateticks[1] is not None:
        rotatey = True
    if rotatey:
        if abs(rotateticks[1]) == 90:
            ax.tick_params(axis='y', pad=0)
        ylabels = ax.get_yticklabels()
        for ylabel in ylabels:
            ylabel.set_ha('right')
            ylabel.set_va('center')
            ylabel.set_rotation(rotateticks[1])
    return ax


def _tickoverlap(ax):
    """
    A small utility function called from many of the plotting functions. This returns booleans
    for whether the x and y tick labels on an Axes overlap as a tuple.

    Parameters:
        ax (ax): A matplotlib axis object with axis tick labels.

    Returns:
        overlaps (tuple bool): Indicators for both axes if the tick labels need to be rotated.

    .. codeauthor:: pygeostat development team - 2015-09-30
    """
    # Set-up a function that completes the check
    def tickcheck(labels):
        """Check to see if labels overlap, returns a boolean value"""
        fig = plt.gcf()
        renderer = fig.canvas.get_renderer()
        fig.draw(renderer)
        if not labels:
            return False
        else:
            try:
                boxes = [label.get_window_extent() for label in labels]
                overlaps = [box.count_overlaps(boxes) for box in boxes]
                return max(overlaps) > 1
            except:
                print("There was an error checking for overlapping axis tick labels. You'll have",
                      "to set the argument manually to rotate them if necessary")
    overlaps = [tickcheck(ax.get_xticklabels()), tickcheck(ax.get_yticklabels())]

    return overlaps
