PyGeoStat
=========

.. raw:: html

   <div class="jumbotron">
     <div class="container">
       <div class="row">
         <div class="col-md-6">
           <p>Welcome to pygeostat, a Python package for geostatistical modeling. Pygeostat is aimed at preparing spatial data, scripting geostatistical workflows, modeling using tools developed at the Centre for Computational Geostatistics (CCG), and constructing visualizations to study spatial data, and geostatistical models.</p>
           <p>
             <a class="btn btn-primary btn-lg" href="user_guide/getting_started.html" role="button">Getting Started</a>
             <a class="btn btn-secondary btn-lg" href="release_highlights.html" role="button">Release Highlights</a>
           </p>
         </div>
         <div class="col-md-6">
           <ul class="feature-list">
             <li>Configurations of persistent project parameters and plotting style</li>
            <li>Data file management functions for interacting with CSV, GeoEAS (GSLIB), and VTK formats</li>
            <li>Utilities for managing GSLIB-style grid definitions</li>
            <li>Export grid and point data files for visualization with Paraview (VTK format)</li>
            <li>Simplified scripting of GSLIB programs with parallelization and crash detection/notification</li>
            <li>Linear desurveying and compositing methods including automatic composite detection</li>
            <li>A collection of example data files included in pygeostat installation <code>gs.ExampleData()</code></li>
            <li>Library of plotting functions</li>
           </ul>
         </div>
       </div>
     </div>
   </div>

Getting Started
--------------

Information to install, test, and contribute to the package.

.. raw:: html

   <div class="container">
     <div class="row">
       <div class="col-md-4">
         <div class="card">
           <div class="card-header">
             <h3>Installation</h3>
           </div>
           <div class="card-body">
             <p>How to install PyGeoStat</p>
             <a class="btn btn-primary" href="user_guide/installation.html">Learn More</a>
           </div>
         </div>
       </div>
       <div class="col-md-4">
         <div class="card">
           <div class="card-header">
             <h3>User Guide</h3>
           </div>
           <div class="card-body">
             <p>The main documentation</p>
             <a class="btn btn-primary" href="user_guide/index.html">Learn More</a>
           </div>
         </div>
       </div>
       <div class="col-md-4">
         <div class="card">
           <div class="card-header">
             <h3>API</h3>
           </div>
           <div class="card-body">
             <p>The detailed API reference</p>
             <a class="btn btn-primary" href="api/index.html">Learn More</a>
           </div>
         </div>
       </div>
     </div>
   </div>

Tutorials Gallery
--------------

A collection of tutorials illustrating the use of PyGeoStat for various geostatistical modeling tasks.

.. raw:: html

   <div class="container">
     <div class="row">
       <div class="col-12">
         <div class="card">
           <div class="card-body">
             <p>Explore our examples to see PyGeoStat in action and learn how to apply it to various geostatistical problems.</p>
             <a class="btn btn-primary" href="examples/index.html">Example Gallery</a>
           </div>
         </div>
       </div>
     </div>
   </div>

.. toctree::
   :maxdepth: 2
   :hidden:

   user_guide/index
   api/index
   examples/index
   contribution/contribution
