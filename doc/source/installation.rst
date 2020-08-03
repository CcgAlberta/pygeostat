
.. _installation:
.. public

Installation
============

Note that pygeostat requires Python 3.6+. Pygeostat will likely require modification to work with
any other version of python. Pygeostat is also dependent on the suite of curated python
packages provided with the Anaconda Python distribution. It is recommended to install this
prior to pygeostat installation.


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

`CCG <https://CCGAlberta.com/>`_ memebers have the option to install CCG/GSLIB software (executable files) to enable pygeostat :ref:`scripting features <programs>`. This can be done using the following function. While GSLIB executable files are avilable through a public repository, for `CCG <https://CCGAlberta.com/>`_ software a vlid access token is required to download executables files from a private repository. The access token is avilable for CCG members at CCG knowledge base website. 

.. autofunction:: pygeostat.utility.get_executable


Create Virtual Environments
----------------------------------

Python environments can be used to isolate different projects with different dependencies. The following code snippets can be used to setup a python environment using `conda <https://docs.conda.io/projects/conda/en/latest/index.html/>`_ package manager. 

1. Creat the environment

.. code-block:: shell

   conda env create -n <environment_name> -f environment.yml

2. Activate the environment

.. code-block:: shell

   conda activate <environment_name>

3. Install pygeostat on the new environment that is activate

.. code-block:: shell

   pip install pygeostat

The nexs step is to launch your favorite IDE for python (e.g. Jupyter).

Using R from Python
*******************
Consider:

    * rpy2 - a Python library providing a low-level interface to R from Python -
      http://rpy.sourceforge.net/
    * Matrix - a R library - https://stat.ethz.ch/R-manual/R-devel/library/base/html/matrix.html


To allow interfacing between Python and R, rpy2 is required. Currently, only the utility gs.nearpd()
uses rpy2 which also requires the R library "Matrix" to be installed through R. Installing rpy2 can
be temperamental. Use the following installation instructions for rpy2:

   1. Download the rpy2 wheel file (e.g., rpy2-2.6.0-cp34-none-win_amd64) from:
      http://www.lfd.uci.edu/~gohlke/pythonlibs/#rpy2
   2. From the command prompt, install the wheel file by using the command

   >>> pip install rpy2-2.6.0-cp34-none-win_amd64.whl

   3. Add the path to the folder containing R.dll to the environment variable PATH
      (e.g., C:\Program Files\R\R-3.1.2/bin\i386)
   4. Add an environment variable R_HOME pointing to R (e.g., C:\Program Files\R\R-3.1.2)
   5. Add an environment variable R_USER that is your windows user name

