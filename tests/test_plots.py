#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for plotting functionality in pygeostat"""

import os
import sys
import warnings
import tempfile
import pytest
import numpy as np
import pandas as pd
import pygeostat as gs
from matplotlib import pyplot as plt

@pytest.fixture
def plot_output_dir():
    """Create temporal output directory for plots"""
    warnings.simplefilter('ignore', category=ImportWarning)
    out_dir = tempfile.NamedTemporaryFile().name
    gs.mkdir(out_dir)
    yield out_dir
    gs.clrmplmem()
    gs.rmdir(out_dir)

@pytest.fixture
def location_plot_data():
    """Provide plot sample data"""
    return gs.ExampleData('point3d_ind_mv')


def test_location_plot_simple(plot_output_dir, location_plot_data):
    """Test capability of plotting by passing just the data file"""
    output_file = os.path.join(plot_output_dir, "location_plot_test.png")
    try:
        _ = gs.location_plot(location_plot_data, output_file = output_file)
    except Exception as ex:
        pytest.fail(f'Unable to test the location plot with simple settings \n{str(ex)}')

    assert os.path.isfile(output_file)


def test_location_plot_cmap(plot_output_dir, location_plot_data):
    """Test capability of plottin by passing argument variable name to get color map"""
    output_file = os.path.join(plot_output_dir, "location_plot_cmap_test.png")
    _, ax = plt.subplots(1,1)

    try:
        _ = gs.location_plot(location_plot_data, var='Phi', output_file=output_file)
        _ = gs.location_plot(location_plot_data, var='Phi', output_file=output_file, ax=ax)

    except Exception as ex:
        pytest.fail(f'Unable to test the location plot with passing variable name\n{str(ex)}')

    assert os.path.isfile(output_file)

def test_location_plot_lines(plot_output_dir, location_plot_data):
    """Test plotting along xz/yz orientation and use the line plots based on drill hole id."""
    output_file = os.path.join(plot_output_dir, "location_plot_lines_test.png")

    try:
        _ = gs.location_plot(
        location_plot_data,
        var = 'Phi',
        orient = 'yz',
        aspect=5,
        plot_collar=True,
        output_file=output_file
    )

    except Exception as ex:
        pytest.fail(f'Unable to test the location plot along yz/xz orientation \n{str(ex)}')

    assert os.path.isfile(output_file)


