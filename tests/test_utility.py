#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for utility functions in pygeostat"""

import os
import sys
import pytest
import pygeostat as gs

def test_token_error():
    """Test if the proper error handling is in place for wrong access token"""
    wrong_token = 'wrong_token'

    with pytest.raises(Exception):
        gs.get_executable(source='ccg', acess_token=wrong_token)


