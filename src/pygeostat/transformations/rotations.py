#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Contains utilities for GSLIB rotations'''

__author__ = 'pygeostat development team'
__date__ = '2016-02-14'
__version__ = '1.000'

import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt

from ..data import DataFile

origin, xaxis, yaxis, zaxis = [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]

def get_rotation_matrix(azm, dip, tilt):
    '''Returns the GSLIB rotation matrix given
       GSLIB angles for azimuth, dip and tilt'''
    from . import transformations as tf
    # Azimuth rotation about Z axis
    gamma = np.radians(azm - 90.0)
    Rz = tf.rotation_matrix(gamma, zaxis)
    # Dip rotation about Y' axis
    beta = np.radians(dip)
    Ry = tf.rotation_matrix(beta, yaxis)
    # Tilt rotation about X' axis
    alpha = np.radians(-1.0*tilt)
    Rx = tf.rotation_matrix(alpha, xaxis)
    # Combination of matrices
    R = tf.concatenate_matrices(Rx, Ry, Rz)
    # Skip the translation piece, just return the rotation matrix
    return(R[0:3,0:3])

def principalvectors(azm, dip, tilt, lefthandrule=True,
                     majorlen=1.0, minorlen=1.0, vertlen=1.0):
    '''Returns the GSLIB principal vectors:
    NOTE: we typically think of things/visualize using the left hand rule in
    which case the original axes are:
    major_axis, minor_axis, vert_axis = [1, 0, 0], [0, -1, 0], [0, 0, 1]
    although the GSLIB rotation matrix actually defines them as:
    major_axis, minor_axis, vert_axis = [1, 0, 0], [0, 1, 0], [0, 0, 1]
    This is only an issue in practice for an asymmetric minor axis'''
    # Rotation matrix - uses the transpose of the GSLIB rotation matrix
    rotmat = get_rotation_matrix(azm, dip, tilt).transpose()
    # Defined axes
    if lefthandrule:
        major_axis, minor_axis, vert_axis = [1, 0, 0], [0, -1, 0], [0, 0, 1]
    else:
        major_axis, minor_axis, vert_axis = [1, 0, 0], [0, 1, 0], [0, 0, 1]

    # Major - oriented + along y for azm=dip=tilt=0
    major_head = np.dot(rotmat, major_axis) * majorlen
    major_xs, major_ys, major_zs = (0, major_head[0]), (0, major_head[1]), (0, major_head[2])

    # Minor - oriented + along y for azm=dip=tilt=0 with Left Hand Rule
    #       - oriented - along y for azm=dip=tilt=0 with Right Hand Rule
    minor_head = np.dot(rotmat, minor_axis) * minorlen
    minor_xs, minor_ys, minor_zs = (0, minor_head[0]), (0, minor_head[1]), (0, minor_head[2])

    # Vert - oriented + along z for azm=dip=tilt=0
    vert_head = np.dot(rotmat, vert_axis) * vertlen
    vert_xs, vert_ys, vert_zs = (0, vert_head[0]), (0, vert_head[1]), (0, vert_head[2])

    return((major_xs, major_ys, major_zs),
           (minor_xs, minor_ys, minor_zs),
           (vert_xs, vert_ys, vert_zs))

def azmdip(point_vector):
    '''Calculates the azimuth/dip to a point (length 3 vector [x,y,z]) from the origin'''
    # Unpack and get angles using atan2 to avoid sign issues
    x, y, z = point_vector
    azm = 90.0 - np.degrees(np.arctan2(y, x))
    dip = np.degrees(np.arctan2(z, np.sqrt(x * x + y * y)))
    return(azm, dip)

def principaldirs(azm, dip, tilt):
    '''Calculate the principal directions from a plane of major continuity'''
    major, minor, vert = principalvectors(azm, dip, tilt)
    major_azm, major_dip = azmdip([major[0][1], major[1][1], major[2][1]])
    minor_azm, minor_dip = azmdip([minor[0][1], minor[1][1], minor[2][1]])
    vert_azm, vert_dip = azmdip([vert[0][1], vert[1][1], vert[2][1]])
    return((major_azm, major_dip), (minor_azm, minor_dip), (vert_azm, vert_dip))

class Arrow3D(FancyArrowPatch):
    '''Arrow3D from http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector'''
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

def plotprincipalvectors(axes, azm, dip, tilt, lefthandrule=True,
                         majorlen=1.0, minorlen=1.0, vertlen=1.0):
    '''Add principal direction unit arrows to plot'''
    import matplotlib.lines as mlines
    # Principal dirs
    major, minor, vert = principalvectors(azm, dip, tilt, lefthandrule,
                                          majorlen=majorlen,
                                          minorlen=minorlen,
                                          vertlen=vertlen)

    # Major
    axes.add_artist(Arrow3D(major[0], major[1], major[2], mutation_scale=20, lw=2,
                            arrowstyle="-|>", color="g"))
    # Minor
    axes.add_artist(Arrow3D(minor[0], minor[1], minor[2], mutation_scale=20, lw=2,
                            arrowstyle="-|>", color="r"))
    # Vert
    axes.add_artist(Arrow3D(vert[0], vert[1], vert[2], mutation_scale=20, lw=2,
                            arrowstyle="-|>", color="b"))

    # Legend parameters
    majorline = mlines.Line2D([], [], color='g', marker=None,
                              markersize=0, lw=2.5, label='Major')
    minorline = mlines.Line2D([], [], color='r', marker=None,
                              markersize=0, lw=2.5, label='Minor')
    vertline = mlines.Line2D([], [], color='b', marker=None,
                             markersize=0, lw=2.5, label='Vert')

    return([majorline, minorline, vertline])

def drawgsvectorwidget(legend=True, sidetext=True, majortitle=False):
    """Draws an interactive GSLIB vector widget in an IPython notebook.

    Warning:
        Must be run with either:

            %matplotlib inline or %matplotlib notebook activated

    .. codeauthor:: Jared Deutsch 2016-02-21"""
    import ipywidgets

    def view_gslib_axes(azm=0.0, dip=0.0, tilt=0.0):
        '''Some function that is called view gslib axes'''

        fig = plt.figure(figsize=(7, 7))
        ax = plt.subplot(111, projection='3d', aspect='equal')

        mlines = plotprincipalvectors(ax, azm, dip, tilt)

        ax.set_xlim((-1.3, 1.3))
        ax.set_ylim((-1.3, 1.3))
        ax.set_zlim((-1.3, 1.3))
        ax.set_xlabel('Easting - X')
        ax.set_ylabel('Northing - Y')
        ax.set_zlabel('Elevation - Z')

        if legend:
            ax.legend(handles=mlines, frameon=False, loc='best')

        if sidetext:
            (major_azm, major_dip), (minor_azm, minor_dip), (vert_azm, vert_dip) = principaldirs(azm, dip, tilt)
            angles = ('Plane of major continuity (azm,dip,tilt) = {:.1f}, {:.1f}, {:.1f}\n\n'
                      'Major (azm,dip) = {:.1f}, {:.1f}\n'
                      'Minor (azm,dip) = {:.1f}, {:.1f}\n'
                      'Vert (azm,dip) = {:.1f}, {:.1f}\n').format(azm, dip, tilt, major_azm,
                                                                  major_dip, minor_azm, minor_dip,
                                                                  vert_azm, vert_dip)

            fig.text(1.0, 0.5, angles, fontsize='large')

        if majortitle:
            ax.set_title('azm = {} ({}), dip = {} ({}), tilt = {}'.format(azm, major_azm, dip,
                                                                          major_dip, tilt))
        plt.show()
        plt.close()
    ipywidgets.interact(view_gslib_axes,
                        azm=(-360.0, 360.0, 1.0),
                        dip=(-180.0, 180.0, 1.0),
                        tilt=(-90.0, 90.0, 1.0))

def drawellipsoid(ax, hmax, hmin, vert, azm=0.0, dip=0.0, tilt=0.0,
                  color='#87CEFA', alpha=0.5):
    """Draws an orientatable ellipsoid with GSLIB conventions.

    Spherical coordinate calculation from:

        http://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib

    Rotation uses GSLIB rotation matrix definition

    .. codeauthor:: Jared Deutsch 2016-03-06"""
    # Coefficients
    rx = hmax
    ry = hmin
    rz = vert

    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))

    # Rotation matrix - uses the transpose of the GSLIB rotation matrix
    rotmat = get_rotation_matrix(azm, dip, tilt).transpose()

    # Flatten and rotate, then reshape
    xyzflat = np.array([x.reshape((100*100)),
                        y.reshape((100*100)),
                        z.reshape((100*100))])
    x, y, z = np.dot(rotmat, xyzflat)
    x = x.reshape((100, 100))
    y = y.reshape((100, 100))
    z = z.reshape((100, 100))

    ax.plot_surface(x, y, z, color=color, alpha=alpha,
                    rstride=4, cstride=4, lw=0.25,
                    antialiased=True)


def drawgsaniswidget(legend=True, sidetext=True, majortitle=False):
    """Draws an interactive anisotropy widget in an IPython notebook.

    Warning:
        Must be run with either:

            %matplotlib inline or %matplotlib notebook activated

    .. codeauthor:: Jared Deutsch 2016-03-06"""
    import ipywidgets

    def view_gslib_axes(azm=0.0, dip=0.0, tilt=0.0,
                        hmax=30.0, hmin=10.0, vert=5.0):
        '''Some function that is called view gslib axes'''

        fig = plt.figure(figsize=(9, 7))
        ax = plt.subplot(111, projection='3d', aspect='equal')

        mlines = plotprincipalvectors(ax, azm, dip, tilt,
                                      majorlen=1.3*max(hmax, hmin, vert),
                                      minorlen=1.3*max(hmax, hmin, vert),
                                      vertlen=1.3*max(hmax, hmin, vert))

        drawellipsoid(ax, hmax, hmin, vert, azm, dip, tilt)

        plotlim = max(ax.get_xlim()[1], ax.get_ylim()[1], ax.get_zlim()[1])
        ax.set_xlim((-1.3 * plotlim, 1.3 * plotlim))
        ax.set_ylim((-1.3 * plotlim, 1.3 * plotlim))
        ax.set_zlim((-1.3 * plotlim, 1.3 * plotlim))
        ax.set_xlabel('Easting - X')
        ax.set_ylabel('Northing - Y')
        ax.set_zlabel('Elevation - Z')

        if legend:
            ax.legend(handles=mlines, frameon=False, loc='best')

        if sidetext:
            (major_azm, major_dip), (minor_azm, minor_dip), (vert_azm, vert_dip) = principaldirs(azm, dip, tilt)
            angles = ('Plane of major continuity (azm,dip,tilt) = {:.1f}, {:.1f}, {:.1f}\n\n'
                      'Major (azm,dip) = {:.1f}, {:.1f}\n'
                      'Minor (azm,dip) = {:.1f}, {:.1f}\n'
                      'Vert (azm,dip) = {:.1f}, {:.1f}\n').format(azm, dip, tilt, major_azm,
                                                                  major_dip, minor_azm, minor_dip,
                                                                  vert_azm, vert_dip)

            fig.text(0.9, 0.5, angles, fontsize='large')

        if majortitle:
            ax.set_title('azm = {} ({}), dip = {} ({}), tilt = {}'.format(azm, major_azm, dip, major_dip, tilt))
        plt.show()
        plt.close()

    ipywidgets.interact(view_gslib_axes,
                        azm=(-360.0, 360.0, 1.0),
                        dip=(-180.0, 180.0, 1.0),
                        tilt=(-90.0, 90.0, 1.0),
                        hmax=(5.0, 1000.0, 5.0),
                        hmin=(5.0, 1000.0, 5.0),
                        vert=(5.0, 1000.0, 5.0))
