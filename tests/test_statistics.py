#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the accuracy plot statistics in pygeostat"""

import pytest
import numpy as np
import pygeostat as gs


@pytest.fixture
def accuracy_data():
    """True values and realizations from the accuracy plot example data"""
    data = gs.ExampleData('accuracy_plot')
    return data[data.columns[0]].values, data[list(data.columns[1:])].values


def test_accuracy_summary_by_hand():
    """Summary statistics for fractions 0.25, 0.75 and 0.5 in the 25%, 50% and 75% intervals"""
    from pygeostat.statistics.utils import _accuracy_summary
    pints = np.array([0.25, 0.5, 0.75])
    propindic = {0.25: [1, 0, 0, 0], 0.5: [1, 1, 1, 0], 0.75: [1, 1, 0, 0]}

    propavg, stats = _accuracy_summary(pints, propindic, variances=[1.0, 3.0])

    np.testing.assert_allclose(propavg['FracIn'], [0.25, 0.75, 0.5])
    assert stats['avgvar'] == pytest.approx(2.0)
    assert stats['mse'] == pytest.approx((0 + 0.25 ** 2 + 0.25 ** 2) / 3)
    assert stats['acc'] == pytest.approx(2 / 3)  # on the line counts as accurate
    assert stats['pre'] == pytest.approx(1 - 2 * 0.25 / 3)
    assert stats['goo'] == pytest.approx(1 - (0.25 + 2 * 0.25) / 3)


def test_accmik_matches_accsim(accuracy_data):
    """MIK probabilities at many thresholds give nearly the same fractions as the realizations"""
    truth, reals = accuracy_data
    thresholds = np.quantile(reals, np.linspace(0.01, 0.99, 99))
    mikprobs = (reals[:, :, None] <= thresholds).mean(axis=1)

    mik, _ = gs.accmik(truth, thresholds, mikprobs)
    sim, _ = gs.accsim(truth, reals)

    np.testing.assert_allclose(mik['FracIn'], sim['FracIn'], atol=0.02)


def test_accmik_below_lowest_threshold():
    """A true value below the lowest threshold has a probability below F there"""
    thresholds = np.array([1.0, 2.0, 3.0])
    # The first CDF already starts at 0, the second needs its lower tail padded
    mikprobs = np.array([[0.0, 0.5, 1.0],
                         [0.2, 0.5, 0.9]])
    truth = np.array([2.0, 0.5])  # p = 0.5, and p = 0.1 on the padded tail

    propavg, _ = gs.accmik(truth, thresholds, mikprobs, pinc=0.05)
    frac = propavg.set_index(propavg['ProbInt'].round(2))['FracIn']

    assert (frac.loc[0.05:0.75] == 0.5).all()  # p = 0.1 is outside these intervals
    assert (frac.loc[0.85:0.95] == 1.0).all()  # and inside these


def test_accmik_above_highest_threshold():
    """A true value above the highest threshold has a probability above F there"""
    thresholds = np.array([1.0, 2.0, 3.0])
    # The first CDF already ends at 1, the second needs its upper tail padded
    mikprobs = np.array([[0.0, 0.5, 1.0],
                         [0.2, 0.5, 0.8]])
    truth = np.array([2.0, 3.5])  # p = 0.5, and p = 0.9 on the padded tail

    propavg, _ = gs.accmik(truth, thresholds, mikprobs, pinc=0.05)
    frac = propavg.set_index(propavg['ProbInt'].round(2))['FracIn']

    assert (frac.loc[0.05:0.75] == 0.5).all()  # p = 0.9 is outside these intervals
    assert (frac.loc[0.85:0.95] == 1.0).all()  # and inside these


def test_accmik_average_variance_is_exact():
    """U is the exact variance of each piecewise linear CDF, the same on every call"""
    truth = np.array([1.0])

    # Uniform on [0, 2]: variance 1/3
    _, uniform = gs.accmik(truth, np.array([0.0, 1.0, 2.0]), np.array([[0.0, 0.5, 1.0]]))
    # 0.8 spread over [0, 1] and 0.2 over [1, 3]: mean 0.8, E[X^2] = 0.8/3 + 0.2 * 13/3
    _, two_step = gs.accmik(truth, np.array([0.0, 1.0, 3.0]), np.array([[0.0, 0.8, 1.0]]))

    assert uniform['avgvar'] == pytest.approx(1 / 3)
    assert two_step['avgvar'] == pytest.approx(0.8 / 3 + 0.2 * 13 / 3 - 0.8 ** 2)


def test_accsim_truth_must_be_one_column(accuracy_data):
    """Passing several columns as the truth is an error, not extra realizations"""
    truth, reals = accuracy_data

    with pytest.raises(ValueError, match='single column'):
        gs.accsim(np.c_[truth, truth], reals)
