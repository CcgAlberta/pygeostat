#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import pandas as pd


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




#class GmmUtilityTest(unittest.TestCase):
#
#    '''
#        This class is used to implement unittests for PyEDW.ServerConnect class
#    '''
#
#    def setUp(self):
#        '''
#        Method called to prepare the test fixture. This is called immediately before calling the test method;
#        other than AssertionError or SkipTest, any exception raised by this method will be considered
#        an error rather than a test failure. The default implementation does nothing.
#        '''
#        x, y = np.random.multivariate_normal([0, 0], [[1, 0], [0, 1]], 5000).T
#        x = x.reshape(5000, 1)
#        y = y.reshape(5000, 1)
#        data = pd.DataFrame(columns=['var1'], data=x)
#        data['var2'] = y
#        self.gmm_util = gs.GmmUtility(mean_vector_list=[[0, 0], [0, 0]], data=data, variable_names=['var1', 'var2'],
#                                   covariance_matrix_list=[[[1, 0], [0, 1]], [[1, 0], [0, 1]]], contribution_list=[0.5, 0.5])
#
#    def test_conditional_mean(self):
#        mean, _, _, _ = self.gmm_util.conditional_moments([0, None])
#
#        self.assertTrue(abs(mean) < 0.00001)
#
#    def test_conditional_variance(self):
#        _, var, _, _ = self.gmm_util.conditional_moments([0, None])
#
#        self.assertTrue(abs(var - 1) < 0.00001)
#
#    def test_conditional_skewness(self):
#        _, _, skew, _ = self.gmm_util.conditional_moments([0, None])
#
#        self.assertTrue(skew < 0.00001)
#
#    def test_conditional_kurtosis(self):
#        _, _, _, kurtosis = self.gmm_util.conditional_moments([0, None])
#
#        self.assertTrue(abs(kurtosis - 3) < 0.00001)
#
#
#if __name__ == '__main__':
#    unittest.main()
