# Pygeostat

<picture align="center">
  <source media="(prefers-color-scheme: dark)" srcset="http://www.ccgalberta.com/pygeostat/_images/pygeostat_logo.png">
  <img alt="Pygeostat Logo" src="http://www.ccgalberta.com/pygeostat/_images/pygeostat_logo.png">
</picture> 

[![PyPI](https://badge.fury.io/py/pygeostat.svg)](https://badge.fury.io/py/pygeostat)
[![Python](https://img.shields.io/pypi/pyversions/pygeostat.svg)](https://pypi.org/project/pygeostat/)
[![CI](https://github.com/CcgAlberta/pygeostat/workflows/IntegrationCheck/badge.svg?branch=master)](https://github.com/CcgAlberta/pygeostat)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

:warning: This package has been updated. Expect breaking changes when migrating
from last stable [version](https://github.com/CcgAlberta/pygeostat/releases/tag/v1.1.1)

## Installation

### Quick Install

```bash
pip install pygeostat
```

### Requirements
- Python 3.10 or higher

### Verify Installation

```{python}
import pygeostat as gs
print(f"Pygeostat version: {gs.__version__}")
print("Basic import successful!")
```

# Test with example data

```{python}
# Load example data (this tests data loading functionality)
data = gs.ExampleData('oilsands')
print(f"Data loaded: {data.shape}")
print(f"Columns: {list(data.columns)}")
print(f"First few rows:\n{data.head()}")
```

## Introduction

This is a Python package for geostatistical modeling. Pygeostat is aimed at preparing spatial data, scripting geostatistical workflows, modeling using tools developed at the Centre for Computational Geostatistics ([CCG](http://www.ccgalberta.com)), and constructing visualizations to study spatial data, and geostatistical models. More information about installing and using pygeostat can be found in the [documentation](http://www.ccgalberta.com/pygeostat/welcome.html).

For lessons on geostatistics visit [Geostatistics Lessons](http://geostatisticslessons.com/).

For a full featured commercial alternative to pygeostat, see [RMSP](https://resourcemodelingsolutions.com/rmsp/) from [Resource Modeling Solutions](https://resourcemodelingsolutions.com).

<a href="https://resourcemodelingsolutions.com/">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://resourcemodelingsolutions.com/static/93cbdd9bde3a60780e21d4ae1c501d18/127cf/Resource-Modeling-Solutions-Home-Page-Logo-KO.webp" width="500px">
  <img alt="Adpative by theme." src="https://resourcemodelingsolutions.com/static/ec83077e0259aa9925747a9199614def/127cf/Resource-Modeling-Solutions-Home-Page-Logo-RGB.webp" width="500px">
</picture>
</a>

Contact [Resource Modeling Solutions](https://resourcemodelingsolutions.com/contact/) about a commercial or academic license.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for detailed release notes.

# Contact:
Refer to [www.ccgalberta.com](http://www.ccgalberta.com).
