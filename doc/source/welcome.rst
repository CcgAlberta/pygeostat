
.. _welcome:

.. image:: ./figures/pygeostat_logo.png
   :align: center

|

|integration_test| |docs| |PyPi| |Python36| |Python37| |Python38|

Welcome
=======

Welcome to pygeostat, a Python package for geostatistical modeling. Pygeostat is aimed at preparing spatial data, scripting geostatistical workflows, modeling using tools developed at the Centre for Computational Geostatistics `(CCG) <http://www.ccgalberta.com>`_, and constructing visualizations to study spatial data, and geostatistical models.

**Features:**

* Configurations of persistent project parameters and plotting style parameters
* Data file management functions for interacting with CSV, GeoEAS (GSLIB), and VTK formats
* Utilities for managing GSLIB-style grid definitions
* Export gridded and point data files for visualization with Paraview (VTK format)
* Simplified scripting of gslib programs with parallelization and crash detection/notification
* Linear desurveying and compositing methods including automatic composite detection
* Vast library of :ref:`plotting functions <plotting>`

Indices and tables
++++++++++++++++++

* :ref:`genindex`
* :ref:`mastertoc`
* :ref:`search`

.. General Package Overview
.. ++++++++++++++++++++++++

.. The pygeostat package is designed with a flat methodology that uses wrappers to tie some of modules and functions together
.. The following figure shows a general layout of the pygeostat package.

..  .. image:: ./figures/gs_overview.jpg

Plotting Gallery
++++++++++++++++

.. raw:: html

    <style type="text/css">
    .rst-content

    .figure {
        position: relative;
        float: left;
        margin-right: 24px;
        width: 170px;
        height: 170px;
        background-position: center;
        background-size: 250px ;
        border: 2px solid black;
        background-repeat: no-repeat;
    }

    .figure img {
        position: absolute;
        display: inline;
        left: 0;
        width: 170px;
        height: 170px;
        opacity:1.0;
        margin-right: 24px;
        filter:alpha(opacity=100); /* For IE8 and earlier */
    }

    .figure:hover background-image {
        -webkit-filter: blur(3px);
        -moz-filter: blur(3px);
        -o-filter: blur(3px);
        -ms-filter: blur(3px);
        filter: blur(3px);
        opacity:1.0;
        filter:alpha(opacity=100); /* For IE8 and earlier */
    }

    .figure span {
        position: absolute;
        display: inline;
        left: 0;
        width: 170px;
        height: 170px;
        background: #000;
        color: #fff;
        visibility: hidden;
        opacity: 0;
        z-index: 100;
    }

    .figure p {
        position: absolute;
        top: 45%;
        width: 170px;
        font-size: 110%;
    }

    .figure:hover span {
        visibility: visible;
        opacity: .4;
    }

    .caption {
        position: absolue;
        width: 180px;
        top: 170px;
        text-align: center !important;
    }
    </style>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.correlation_matrix_plot>
    <div style="background-image: url(./_images/plotting-24.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>correlation_matrix_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.slice_plot>
    <div style="background-image: url(./_images/plotting-13.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center' margin=10px>
    <span class='figure-label'>
    <p>slice_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.scatter_plot>
    <div style="background-image: url(./_images/plotting-15_00.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>scatter_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.scatter_plots>
    <div style="background-image: url(./_images/plotting-22.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>scatter_plots</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.scatter_plots_lu>
    <div style="background-image: url(./_images/plotting-23_01.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>scatter_plots_lu</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.pit_plot>
    <div style="background-image: url(./_images/pitplot_mr.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>pit_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.accuracy_plot>
    <div style="background-image: url(./_images/plotting-1.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>accuracy_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.variogram_plot>
    <div style="background-image: url(./_images/plotting-37.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>variogram_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.drill_plot>
    <div style="background-image: url(./_images/plotting-28.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>drill_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.qq_plot>
    <div style="background-image: url(./_images/plotting-30.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>qq_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.validation_plot>
    <div style="background-image: url(./_images/plotting-33.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>validation_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.histogram_plot>
    <div style="background-image: url(./_images/plotting-9.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>histogram_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.histogram_plot_simulation>
    <div style="background-image: url(./_images/plotting-12.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>histogram_plot_simulation</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.location_plot>
    <div style="background-image: url(./_images/plotting-5.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>location_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.probability_plot>
    <div style="background-image: url(./_images/plotting-31.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>probability_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.loadings_plot>
    <div style="background-image: url(./_images/plotting-38.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>loadings_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <a href=./plotting.html#pygeostat.plotting.contour_plot>
    <div style="background-image: url(./_images/plotting-40.png); background-repeat:no-repeat; background-size: contain"
         class='figure align-center'>
    <span class='figure-label'>
    <p>contour_plot</p>
    </span>
    </div>
    </a>

.. raw:: html

    <div style="clear: both"></div>

Terms of Use
++++++++++++
pygeostat is licensed under the CCG Terms of Use, which may be found at the below link.
http://www.ccgalberta.com/software-terms-of-use/





.. |integration_test| image:: https://github.com/MHadavand/pygeostat_public/workflows/IntegrationCheck/badge.svg?branch=master
    :alt: Integration
    :target: https://github.com/MHadavand/pygeostat_public

.. |docs| image:: https://github.com/MHadavand/pygeostat_public/workflows/Documentation/badge.svg?branch=master
    :alt: Documentation Status
    :target: https://github.com/MHadavand/pygeostat_public

.. |PyPi| image:: https://badge.fury.io/py/pygeostat.svg
    :target: https://badge.fury.io/py/pygeostat

.. |Python36| image:: https://img.shields.io/badge/python-3.6-cyan.svg
    :target: https://www.python.org/downloads/release/python-360
    :alt: Python Version

.. |Python37| image:: https://img.shields.io/badge/python-3.7-red.svg
    :target: https://www.python.org/downloads/release/python-370
    :alt: Python Version

.. |Python38| image:: https://img.shields.io/badge/python-3.8-black.svg
    :target: https://www.python.org/downloads/release/python-380
    :alt: Python Version