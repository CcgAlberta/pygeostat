'''pygeostat: Python interface for geostatistics
 Copyright (C) 2013-2020 pygeostat core team'''

import os
import subprocess
import setuptools
from setuptools import setup

with open("README.rst", "r") as fh:
    long_description = fh.read()

with open ('requirements.txt', 'r') as file:
    install_requires = file.readlines()
install_requires = [x.strip() for x in install_requires]


setup(name='pygeostat',
      version='1.0.3',
      description='Python interface for geostatistics',
      long_description = long_description,
      long_description_content_type="text/x-rst",
      url='http://www.ccgalberta.com/pygeostat_public/welcome.html',
      project_urls ={
          'Documentation': 'http://www.ccgalberta.com/pygeostat_public/welcome.html'
          ,'Source': 'https://github.com/CcgAlberta/pygeostat'
      },
      maintainer='pygeostat development team',
      maintainer_email='hadavand@ualberta.ca',
      author='Jared Deutsch, Matthew Deutsch, Warren Black, Tyler Acorn, Ryan Martin, Mostafa Hadavand, Ryan Barnett',
      license='CCG Software Terms and Conditions',
      packages=setuptools.find_packages(),
      python_requires='>3.5',
      install_requires = install_requires,
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3", 
        "Operating System :: Microsoft", 
       ],
      package_data={'pygeostat': ['data/example_data/*', 'executable/*.exe']},
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

# Release Steps
# Upgrade setuptools
# >>> python -m pip install --user --upgrade setuptools wheel
# Create the distribution files
# >>> python setup.py sdist bdist_wheel
# Install Twine
# >>> python -m pip install --user --upgrade twine
# Upload the package tp test PyPi
# python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* (e.g. https://test.pypi.org/project/pygeostat/1.0.0/)
# Install the uploaded packeg to test PyPi and test on an environment
# pip install -i https://test.pypi.org/simple/ pygeostat
# After testing the test version upload it to PyPi
# twine upload dist/*