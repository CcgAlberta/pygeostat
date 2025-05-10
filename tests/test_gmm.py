#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import pandas as pd

from pygeostat.plotting.gaussian_mv import GmmUtility

# Fixture component to be use on the tests
@pytest.fixture
def gmm_util():
    """Create a GMM utility fixture with standard normal test data."""
    x, y = np.random.multivariate_normal(
        mean=[0, 0], 
        cov=[[1, 0], [0, 1]], 
        size=5000
    ).T
    
    data = pd.DataFrame({
        'var1': x,
        'var2': y
    })
    
    return GmmUtility(
        mean_vector_list=[[0, 0], [0, 0]], 
        data=data, 
        variable_names=['var1', 'var2'],
        covariance_matrix_list=[[[1, 0], [0, 1]], [[1, 0], [0, 1]]], 
        contribution_list=[0.5, 0.5]
    )

def test_conditional_mean(gmm_util):
    """Test that conditional mean is correctly computed."""
    mean, _, _, _ = gmm_util.conditional_moments([0, None])
    assert abs(mean) < 1e-5, "Conditional mean should be approximately zero"

def test_conditional_variance(gmm_util):
    """Test that conditional variance is correctly computed."""
    _, var, _, _ = gmm_util.conditional_moments([0, None])
    assert abs(var - 1) < 1e-5, "Conditional variance should be approximately one"

def test_conditional_skewness(gmm_util):
    """Test that conditional skewness is correctly computed."""
    _, _, skew, _ = gmm_util.conditional_moments([0, None])
    assert abs(skew) < 1e-5, "Conditional skewness should be approximately zero"

def test_conditional_kurtosis(gmm_util):
    """Test that conditional kurtosis is correctly computed."""
    _, _, _, kurtosis = gmm_util.conditional_moments([0, None])
    assert abs(kurtosis - 3) < 1e-5, "Conditional kurtosis should be approximately three"

