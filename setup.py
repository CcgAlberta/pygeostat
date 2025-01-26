'''pygeostat: Python interface for geostatistics
 Copyright (C) 2013-2020 pygeostat core team'''

import os
import subprocess
import setuptools
from setuptools import setup

with open("README.rst", "r") as fh:
    long_description = fh.read()

with open('requirements.txt', 'r') as file:
    install_requires = file.readlines()
install_requires = [x.strip() for x in install_requires]


setup(name='pygeostat',
      version='1.1.1',
      description='Python interface for geostatistics',
      long_description=long_description,
      long_description_content_type="text/x-rst",
      url='http://www.ccgalberta.com/pygeostat/welcome.html',
      project_urls={
          'Documentation': 'http://www.ccgalberta.com/pygeostat/welcome.html', 'Source': 'https://github.com/CcgAlberta/pygeostat'
      },
      maintainer='pygeostat development team',
      maintainer_email='hadavand@ualberta.ca',
      author='Jared Deutsch, Matthew Deutsch, Ryan Martin, Warren Black, Tyler Acorn, Mostafa Hadavand, Ryan Barnett',
      license='MIT License',
      packages=setuptools.find_packages()+ ['pygeostat.fortran.src',
                                                 'pygeostat.fortran.src.subs',
                                                 'pygeostat.fortran.resource',
                                                 'pygeostat.fortran.resource.LAPACK',
                                                 'pygeostat.fortran.resource.LAPACK.ALL',
                                                 'pygeostat.fortran.resource.rm'],
      python_requires='>=3.10',
      install_requires=install_requires,
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: MIT License",
          "Operating System :: Microsoft",
      ],  
	  package_data={'': ['*.for', '*.bat', '*.f90', '*.f', 'rm.exe',
                    'libintl3.dll', 'libiconv2.dll',
                    'License.txt'], 'pygeostat': ['data/example_data/*', 'executable/*.txt']},
      include_package_data=True,
      zip_safe=False)


# To install, run python with the following command line:
# >>> pip install .
# Alternative installation, run python with the following command line:
# python setup.py install
# To uninstall, run the command:
# >>> pip uninstall pygeostat
# To create distribution files
# >>> python setup.py bdist_wheel
