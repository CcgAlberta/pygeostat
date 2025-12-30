#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for Parameters class and validation functions in pygeostat.

This module contains comprehensive tests for:
- Parameters class (get, set, update, defaults)
- All validation functions (especially those using six library)
- Dict operations (iteritems usage)
- String type checking (string_types usage)
"""

import pytest
from pygeostat.pygeostat_parameters import (
    Parameters,
    _validate_string,
    _validate_string_or_None,
    _validate_float,
    _validate_float_or_None,
    _validate_int,
    _validate_int_or_None,
    _validate_bool,
    _validate_dict,
    _validate_dict_or_None,
    _validate_list,
    _validate_list_or_None,
    _validate_dict_or_string,
    _validate_color,
    _validate_kde_or_color,
)

# Fixture to be use in the tests
@pytest.fixture
def param_test_values():
    """Provide test values for parameter testing."""
    return {
        'key': 'data.tmin',
        'value': -789,
        'wrong_value': 'wrong'
    }

def test_parameter_descriptions():
    """
    Test that parameters have proper descriptions.
    
    This test verifies:
    1. All parameters exist in the _default_descriptions dictionary
    2. All descriptions are non-empty strings

    #TODO: Ref #96
    """
    # Get all parameters
    all_params = list(Parameters)
    
    # Get all parameters that have descriptions
    described_params = set(Parameters._default_descriptions.keys())
    
    # Find parameters without descriptions
    params_without_desc = set(all_params) - described_params
    
    # List parameters with empty descriptions
    params_with_empty_desc = [
        param for param in described_params 
        if not Parameters._default_descriptions[param] or 
           not isinstance(Parameters._default_descriptions[param], str)
    ]
    
    # Create detailed error message
    error_message = ""
    if params_without_desc:
        error_message += f"\nParameters missing from _default_descriptions ({len(params_without_desc)}):\n"
        error_message += ", ".join(sorted(params_without_desc))
    
    if params_with_empty_desc:
        error_message += f"\n\nParameters with empty descriptions ({len(params_with_empty_desc)}):\n"
        error_message += ", ".join(sorted(params_with_empty_desc))
    
    assert not params_without_desc and not params_with_empty_desc, \
        f"Parameter description issues found:{error_message}\n\n" \
        f"Total parameters: {len(all_params)}, With descriptions: {len(described_params)}"

def test_parameter_update(param_test_values):
    """Test parameter update functionality with valid and invalid values."""

    key = param_test_values['key']
    value = param_test_values['value']
    
    # Save original value to restore later
    original_value = Parameters[key]
    
    try:
        # Update parameter with valid value
        Parameters.update({key: value})
        assert Parameters[key] == value, f"Parameter {key} was not updated correctly"
        
        # Test correct key but wrong value type
        with pytest.raises(ValueError):
            Parameters.update({key: param_test_values['wrong_value']})
    finally:
        # Restore original value
        Parameters.update({key: original_value})

def test_system_default_operations():
    """Test saving and retrieving system default parameters."""
    # Get current values to restore later
    original_values = {key: Parameters[key] for key in Parameters}
    
    try:
        # Test setting system defaults
        Parameters.set_systemdefault()
        
        # Test getting system defaults
        Parameters.get_systemdefault()
        # No need for assertion as we just want to ensure it doesn't error
    finally:
        # Restore original values
        Parameters.update(original_values)

@pytest.mark.parametrize("param_name,expected_type", [
    ('data.tmin', float),
    ('plotting.cmap', str),
    ('config.verbose', bool),
])

def test_parameter_types(param_name, expected_type):
    """Test that key parameters have the expected types."""
    assert isinstance(Parameters[param_name], expected_type), \
        f"Parameter {param_name} should be of type {expected_type}"

def test_find_all():
    """Test the find_all method for searching parameters by pattern."""
    plot_params = Parameters.find_all('plotting')
    assert plot_params, "Should find plotting parameters"
    assert all('plotting' in key for key in plot_params.keys()), \
        "All keys should contain 'plotting'"

def test_set_get_system_default():
    """Test setting and getting system defaults."""
    # Store original value
    original_value = Parameters['plotting.cmap']
    temp_value = 'plasma' if original_value != 'plasma' else 'viridis'

    try:
        # Set a new value and save to system defaults
        Parameters['plotting.cmap'] = temp_value
        Parameters.set_systemdefault()

        # Reset to original, then load from system defaults
        Parameters['plotting.cmap'] = original_value
        Parameters.get_systemdefault()

        # The value should now be the temp_value we set
        assert Parameters['plotting.cmap'] == temp_value, \
            "System default values were not correctly saved or loaded"
    finally:
        # Restore original value
        Parameters['plotting.cmap'] = original_value
        # Save original value back to system defaults
        Parameters.set_systemdefault()


# Tests for validation functions (especially those using six library)

class TestValidateString:
    """Tests for _validate_string function (uses six.text_type)."""

    def test_validate_string_basic(self):
        """Test basic string validation."""
        assert _validate_string('test') == 'test'
        assert _validate_string('') == ''
        assert _validate_string('123') == '123'

    def test_validate_string_int(self):
        """Test string validation with integers."""
        assert _validate_string(123) == '123'
        assert _validate_string(0) == '0'

    def test_validate_string_float(self):
        """Test string validation with floats."""
        assert _validate_string(3.14) == '3.14'

    def test_validate_string_list(self):
        """Test string validation with lists (uses six.text_type)."""
        result = _validate_string(['a', 'b', 'c'])
        assert result == ['a', 'b', 'c']
        assert all(isinstance(item, str) for item in result)

        # Test list with integers
        result = _validate_string([1, 2, 3])
        assert result == ['1', '2', '3']

    def test_validate_string_dict(self):
        """Test string validation with dicts (uses six.text_type)."""
        result = _validate_string({'key': 'value', 'num': 123})
        assert isinstance(result, dict)
        assert all(isinstance(k, str) for k in result.keys())
        assert all(isinstance(v, str) for v in result.values())

    def test_validate_string_none_not_accepted(self):
        """Test that None converts to 'None' string when not accepted."""
        # When accept_none=False, six.text_type(None) returns 'None' string
        result = _validate_string(None, accept_none=False)
        assert result == 'None'
        assert isinstance(result, str)

    def test_validate_string_none_accepted(self):
        """Test that None is accepted when specified."""
        assert _validate_string(None, accept_none=True) is None

    def test_validate_string_or_none(self):
        """Test _validate_string_or_None wrapper."""
        assert _validate_string_or_None(None) is None
        assert _validate_string_or_None('test') == 'test'


class TestValidateFloat:
    """Tests for _validate_float function."""

    def test_validate_float_basic(self):
        """Test basic float validation."""
        assert _validate_float(3.14) == 3.14
        assert _validate_float(0.0) == 0.0
        assert _validate_float(-1.5) == -1.5

    def test_validate_float_from_int(self):
        """Test float validation from integers."""
        assert _validate_float(5) == 5.0
        assert isinstance(_validate_float(5), float)

    def test_validate_float_from_string(self):
        """Test float validation from string."""
        assert _validate_float('3.14') == 3.14
        assert _validate_float('-1.5') == -1.5

    def test_validate_float_invalid(self):
        """Test float validation with invalid input."""
        with pytest.raises(ValueError):
            _validate_float('not_a_number')

    def test_validate_float_none(self):
        """Test float validation with None."""
        assert _validate_float(None, accept_none=True) is None
        with pytest.raises((ValueError, TypeError)):
            _validate_float(None, accept_none=False)

    def test_validate_float_or_none(self):
        """Test _validate_float_or_None wrapper."""
        assert _validate_float_or_None(None) is None
        assert _validate_float_or_None(3.14) == 3.14


class TestValidateInt:
    """Tests for _validate_int function."""

    def test_validate_int_basic(self):
        """Test basic int validation."""
        assert _validate_int(5) == 5
        assert _validate_int(0) == 0
        assert _validate_int(-10) == -10

    def test_validate_int_from_float(self):
        """Test int validation from float."""
        assert _validate_int(3.14) == 3
        assert _validate_int(5.9) == 5

    def test_validate_int_from_string(self):
        """Test int validation from string."""
        assert _validate_int('42') == 42
        assert _validate_int('-5') == -5

    def test_validate_int_invalid(self):
        """Test int validation with invalid input."""
        with pytest.raises(ValueError):
            _validate_int('not_a_number')

    def test_validate_int_none(self):
        """Test int validation with None."""
        assert _validate_int(None, accept_none=True) is None
        with pytest.raises((ValueError, TypeError)):
            _validate_int(None, accept_none=False)

    def test_validate_int_or_none(self):
        """Test _validate_int_or_None wrapper."""
        assert _validate_int_or_None(None) is None
        assert _validate_int_or_None(42) == 42


class TestValidateBool:
    """Tests for _validate_bool function."""

    def test_validate_bool_basic(self):
        """Test basic boolean validation."""
        assert _validate_bool(True) is True
        assert _validate_bool(False) is False

    def test_validate_bool_from_int(self):
        """Test boolean validation from integers."""
        assert _validate_bool(1) is True
        assert _validate_bool(0) is False
        assert _validate_bool(5) is True  # Any non-zero is True

    def test_validate_bool_from_string(self):
        """Test boolean validation from strings."""
        assert _validate_bool('test') is True
        assert _validate_bool('') is False  # Empty string is False


class TestValidateDict:
    """Tests for _validate_dict function."""

    def test_validate_dict_basic(self):
        """Test basic dict validation."""
        test_dict = {'key': 'value', 'num': 123}
        result = _validate_dict(test_dict)
        assert result == test_dict
        assert isinstance(result, dict)

    def test_validate_dict_empty(self):
        """Test dict validation with empty dict."""
        result = _validate_dict({})
        assert result == {}
        assert isinstance(result, dict)

    def test_validate_dict_invalid(self):
        """Test dict validation with non-dict input."""
        with pytest.raises(ValueError):
            _validate_dict('not a dict')
        with pytest.raises(ValueError):
            _validate_dict([1, 2, 3])

    def test_validate_dict_none(self):
        """Test dict validation with None."""
        assert _validate_dict(None, accept_none=True) is None
        with pytest.raises(ValueError):
            _validate_dict(None, accept_none=False)

    def test_validate_dict_or_none(self):
        """Test _validate_dict_or_None wrapper."""
        assert _validate_dict_or_None(None) is None
        assert _validate_dict_or_None({'key': 'value'}) == {'key': 'value'}


class TestValidateList:
    """Tests for _validate_list function."""

    def test_validate_list_basic(self):
        """Test basic list validation."""
        result = _validate_list([1, 2, 3])
        assert result == [1, 2, 3]
        assert isinstance(result, list)

    def test_validate_list_from_tuple(self):
        """Test list validation from tuple."""
        result = _validate_list((1, 2, 3))
        assert result == [1, 2, 3]
        assert isinstance(result, list)

    def test_validate_list_from_string(self):
        """Test list validation from string (converts to list of chars)."""
        result = _validate_list('abc')
        assert result == ['a', 'b', 'c']

    def test_validate_list_none(self):
        """Test list validation with None."""
        assert _validate_list(None, accept_none=True) is None
        with pytest.raises((ValueError, TypeError)):
            _validate_list(None, accept_none=False)

    def test_validate_list_or_none(self):
        """Test _validate_list_or_None wrapper."""
        assert _validate_list_or_None(None) is None
        assert _validate_list_or_None([1, 2, 3]) == [1, 2, 3]


class TestValidateDictOrString:
    """Tests for _validate_dict_or_string function."""

    def test_validate_dict_or_string_with_string(self):
        """Test with string input."""
        result = _validate_dict_or_string('test_string')
        assert result == 'test_string'
        assert isinstance(result, str)

    def test_validate_dict_or_string_with_dict(self):
        """Test with dict input."""
        test_dict = {'key': 'value'}
        result = _validate_dict_or_string(test_dict)
        assert result == test_dict
        assert isinstance(result, dict)

    def test_validate_dict_or_string_with_int(self):
        """Test with int (should convert to string)."""
        result = _validate_dict_or_string(123)
        assert result == '123'
        assert isinstance(result, str)


class TestValidateColor:
    """Tests for _validate_color function (uses six.string_types)."""

    def test_validate_color_named(self):
        """Test with named colors."""
        assert _validate_color('red') == 'red'
        assert _validate_color('blue') == 'blue'

    def test_validate_color_hex_with_hash(self):
        """Test with hex colors (with #)."""
        assert _validate_color('#FF0000') is not None
        assert _validate_color('#00FF00') is not None

    def test_validate_color_hex_without_hash(self):
        """Test with hex colors (without #) - uses isinstance(s, six.string_types)."""
        result = _validate_color('FF0000')
        assert result == '#FF0000'  # Should add # prefix

    def test_validate_color_none_string(self):
        """Test with 'none' string."""
        assert _validate_color('none') == 'none'
        assert _validate_color('None') == 'none'
        assert _validate_color('NONE') == 'none'

    def test_validate_color_rgb_tuple(self):
        """Test with RGB tuple."""
        result = _validate_color((1.0, 0.0, 0.0))
        assert result == (1.0, 0.0, 0.0)


class TestValidateKdeOrColor:
    """Tests for _validate_kde_or_color function (uses six.string_types)."""

    def test_validate_kde_or_color_kde(self):
        """Test with 'kde' string."""
        assert _validate_kde_or_color('kde') == 'kde'
        assert _validate_kde_or_color('KDE') == 'kde'

    def test_validate_kde_or_color_none(self):
        """Test with 'none' string."""
        assert _validate_kde_or_color('none') == 'none'
        assert _validate_kde_or_color('None') == 'none'

    def test_validate_kde_or_color_hex(self):
        """Test with hex color (uses isinstance(s, six.string_types))."""
        result = _validate_kde_or_color('FF0000')
        assert result == '#FF0000'

    def test_validate_kde_or_color_named(self):
        """Test with named color."""
        result = _validate_kde_or_color('red')
        assert result == 'red'


# Tests for Parameters dict operations (uses six.iteritems)

class TestParametersDictOperations:
    """Tests for Parameters dict operations (uses six.iteritems in iteration)."""

    def test_parameters_iteration(self):
        """Test iterating over parameters."""
        # This internally uses six.iteritems in Parameters class
        all_params = list(Parameters)
        assert len(all_params) > 0
        assert all(isinstance(key, str) for key in all_params)

    def test_parameters_items(self):
        """Test getting parameter items."""
        # Get a subset of parameters
        plot_params = Parameters.find_all('plotting')
        assert isinstance(plot_params, dict)
        assert len(plot_params) > 0

        # Verify all keys and values
        for key, value in plot_params.items():
            assert isinstance(key, str)
            assert 'plotting' in key

    def test_parameters_dict_iteration(self):
        """Test dict-style iteration over parameters (uses six.iteritems)."""
        # Iterate using dict() to convert Parameters to dict
        params_dict = dict(Parameters)
        assert isinstance(params_dict, dict)
        assert len(params_dict) > 0

        # Verify structure
        for key, value in params_dict.items():
            assert isinstance(key, str)

    def test_parameters_keys_values(self):
        """Test accessing keys and values."""
        keys = list(Parameters.keys())
        values = list(Parameters.values())

        assert len(keys) > 0
        assert len(values) > 0
        assert len(keys) == len(values)


# Tests for string type checking with Parameters

class TestParametersStringHandling:
    """Tests for Parameters handling of string types (six.text_type usage)."""

    def test_parameter_string_value(self):
        """Test setting and getting string parameter."""
        original = Parameters['plotting.cmap']

        try:
            # Test with regular string
            Parameters['plotting.cmap'] = 'viridis'
            assert Parameters['plotting.cmap'] == 'viridis'
            assert isinstance(Parameters['plotting.cmap'], str)

            # Test with different string
            Parameters['plotting.cmap'] = 'plasma'
            assert Parameters['plotting.cmap'] == 'plasma'
        finally:
            Parameters['plotting.cmap'] = original

    def test_parameter_update_with_dict(self):
        """Test updating multiple parameters with dict."""
        original_cmap = Parameters['plotting.cmap']
        original_verbose = Parameters['config.verbose']

        try:
            # Update multiple parameters
            Parameters.update({
                'plotting.cmap': 'plasma',
                'config.verbose': False
            })

            assert Parameters['plotting.cmap'] == 'plasma'
            assert Parameters['config.verbose'] is False
        finally:
            Parameters.update({
                'plotting.cmap': original_cmap,
                'config.verbose': original_verbose
            })


