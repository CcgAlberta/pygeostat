#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from pygeostat.pygeostat_parameters import Parameters

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


