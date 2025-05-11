#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from pygeostat.utility import get_executable

def test_token_error():
    """Test if the proper error handling is in place for wrong access token"""
    wrong_token = 'wrong_token'

    with pytest.raises(Exception):
        get_executable(source='ccg', acess_token=wrong_token)


