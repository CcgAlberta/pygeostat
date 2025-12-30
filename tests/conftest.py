#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Pytest configuration for pygeostat tests.

This file configures the test environment before any tests run.
"""

import matplotlib
import os

# Set matplotlib to use non-GUI backend for CI/headless environments
# This prevents TkInter crashes on Windows CI runners
matplotlib.use('Agg')

# Suppress matplotlib GUI warnings
os.environ['MPLBACKEND'] = 'Agg'
