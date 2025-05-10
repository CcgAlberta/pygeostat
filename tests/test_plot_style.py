#!/usr/bin/env python
# -*- coding: utf-8 -*-

from _pytest.capture import capfd
import pytest
import warnings

from pygeostat.plotting.set_style import PlotStyle

# Fixture component used in the tests
@pytest.fixture(scope="module")
def plot_style():
    """
    Fixture providing the PlotStyle object from gs for tests

    Returns the PlotStyle class to test its static methods and properties
    """

    warnings.simplefilter('ignore', category=ImportWarning)
    return PlotStyle

def test_update_key(plot_style):
    """Test to update plot style parameters with valid keys."""

    # Original so it can be restored later
    original_value = plot_style.get('axes.grid', None)

    try:
        # Update with valid key
        plot_style.update({'axes.grid': True})
        assert plot_style['axes.grid'] is True
    finally:
        # Restore orignal value to avoid affect rest of tests
        if original_value is not None:
            plot_style.update({'axes.grid': original_value})

def test_update_invalid_key(plot_style):
    """Test updating plot style parameters with invalid keys raises KeyError."""
    with pytest.raises(KeyError):
        plot_style.update({'FakeParameters': 0})

def test_set_default_plot_style(plot_style):
    """Test setting system default plot styles."""
    try:
        plot_style.set_systemdefault()
    except Exception as e:
        pytest.fail(f"Setting system defaults failed with: {str(e)}")

def test_get_default_plot_style(plot_style, capfd):
    """Test getting the default plot style."""
    result = plot_style.get_systemdefault()
    captured = capfd.readouterr()
    assert "Loading default Pygeostat Parameters from" in captured.out
