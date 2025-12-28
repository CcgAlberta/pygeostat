PyGeoStat
=========

A Python package for geostatistical modeling and spatial uncertainty analysis.

.. admonition:: What is PyGeoStat?
   :class: note

   PyGeoStat is aimed at preparing spatial data, scripting geostatistical
   workflows, modeling using tools developed at the Centre for Computational
   Geostatistics (CCG), and constructing visualizations to study spatial data
   and geostatistical models.

Key Features
------------

- Persistent project configuration and plotting styles
- Data file management for CSV, GeoEAS (GSLIB), and VTK formats
- Utilities for GSLIB-style grid definitions
- Export to Paraview (VTK)
- Parallelized GSLIB workflow scripting with crash detection
- Desurveying and compositing utilities
- Example datasets via ``gs.ExampleData()``
- Comprehensive plotting utilities

Getting Started
---------------

.. grid:: 3
   :gutter: 3

   .. grid-item-card:: 🚀 Installation
      :link: user_guide/installation
      :link-type: doc

      How to install PyGeoStat and verify your setup.

   .. grid-item-card:: 📘 User Guide
      :link: user_guide/index
      :link-type: doc

      Core concepts and geostatistical workflows.

   .. grid-item-card:: 🧩 API Reference
      :link: api/index
      :link-type: doc

      Detailed documentation of all classes and functions.

Tutorials
---------

.. grid:: 1

   .. grid-item-card:: 🧪 Examples Gallery
      :link: examples/index
      :link-type: doc

      Hands-on tutorials illustrating PyGeoStat in real-world workflows.

.. toctree::
   :maxdepth: 2
   :hidden:

   user_guide/index
   api/index
   examples/index
   contribution/contribution

