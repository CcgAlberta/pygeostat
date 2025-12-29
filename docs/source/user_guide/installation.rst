
.. _installation:
.. public

Installation
============

Note that pygeostat requires Python 3.10 or higher. PyGeoStat follows modern Python standards
and requires features available in Python 3.10+.

.. important::
   **Minimum Python Version:** 3.10+

   For historical versions supporting older Python, see pygeostat version 1.1.1.


Python Installation
-------------------

The `Anaconda <https://www.anaconda.com>`_ distribution by Continuum Analytics is recommended. If unfamiliar with Python and virtual environments,
then accepting the defaults of registering with the system and installing to the path are
recommended.

Pygeostat Installation
--------------------------

`Pygeostat <https://pypi.org/project/pygeostat/>`_ can be installed from python package index (PyPI) repository.

.. code-block:: shell

   pip install pygeostat

Additional Software (`CCG <https://CCGAlberta.com/>`_ Members)
****************************************************************

`CCG <https://CCGAlberta.com/>`_ members have the option to install CCG/GSLIB software (executable files) to enable pygeostat :ref:`scripting features <programs>`. This can be done using the following function. While GSLIB executable files are available through a public repository, for `CCG <https://CCGAlberta.com/>`_ software a valid access token is required to download executables files from a private repository. The access token is available for CCG members at CCG knowledge base website. 

.. autofunction:: pygeostat.utility.get_executable


Create Virtual Environments
----------------------------------

Python environments can be used to isolate different projects with different dependencies. The following code snippets can be used to setup a python environment using `conda <https://docs.conda.io/projects/conda/en/latest/index.html/>`_ package manager. 

1. Create the environment

.. code-block:: shell

   conda create -n <environment_name> python=3.10

2. Activate the environment

.. code-block:: shell

   conda activate <environment_name>

3. Install pygeostat in the activated environment

.. code-block:: shell

   pip install pygeostat

The next step is to launch your favorite IDE for python (e.g. Jupyter).

Optional: Using R from Python
*******************************

.. note::
   **This is optional.** rpy2 is only required for specific functions like ``gs.nearpd()``.
   Most pygeostat functionality does not require R.

If you need to interface between Python and R:

**Requirements:**
   * R installed on your system
   * rpy2 - Python interface to R
   * Matrix - R library (install via R)

**Installation:**

.. code-block:: shell

   pip install rpy2

For platform-specific requirements and troubleshooting, see the `rpy2 documentation <https://rpy2.github.io/>`_.

