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
    """A, P and G for fractions 0.25, 0.75 and 0.5 in the 25%, 50% and 75% intervals, with local
    variances of 0.2 and 0.4 in Gaussian units (global variance 1)"""
    from pygeostat.statistics.utils import _accuracy_summary
    pints = np.array([0.25, 0.5, 0.75])
    propindic = {0.25: [1, 0, 0, 0], 0.5: [1, 1, 1, 0], 0.75: [1, 1, 0, 0]}

    propavg, stats = _accuracy_summary(pints, propindic, variances=[1.0, 3.0],
                                       gaussian_variances=[0.2, 0.4])

    np.testing.assert_allclose(propavg['FracIn'], [0.25, 0.75, 0.5])
    assert stats['avgvar'] == pytest.approx(2.0)
    assert stats['mse'] == pytest.approx((0 + 0.25 ** 2 + 0.25 ** 2) / 3)
    # Only the 75% interval is below the 1:1 line, by 0.25: A = (1 - 2 * 0.25 / 3) * 100
    assert stats['acc'] == pytest.approx(100 * (1 - 2 * 0.25 / 3))
    assert stats['pre'] == pytest.approx(70.0)
    assert stats['goo'] == pytest.approx(2 / (1 / stats['acc'] + 1 / 70.0))
    assert min(stats['acc'], stats['pre']) < stats['goo'] < max(stats['acc'], stats['pre'])


def test_accuracy_summary_limits():
    """A is 0 when no true value is in any interval and 100 when all intervals are on or above the
    line; P is negative (not clipped) when the local variance exceeds the global one, and G is then
    not defined; G is 0 when A is 0; P and G are NaN without Gaussian variances."""
    from pygeostat.statistics.utils import _accuracy_summary
    pints = np.array([0.25, 0.5, 0.75])
    none_in = {p: [0, 0] for p in pints}
    all_in = {p: [1, 1] for p in pints}

    _, nothing = _accuracy_summary(pints, none_in, [1.0], gaussian_variances=[0.5])
    assert nothing['acc'] == pytest.approx(100 * (1 - 2 * 0.5))      # mean of the intervals is 0.5
    assert nothing['acc'] == 0 and nothing['goo'] == 0 and nothing['pre'] == pytest.approx(50)
    _, perfect = _accuracy_summary(pints, all_in, [1.0], gaussian_variances=[0.0])
    assert perfect['acc'] == 100 and perfect['pre'] == 100 and perfect['goo'] == pytest.approx(100)
    _, wide = _accuracy_summary(pints, all_in, [1.0], gaussian_variances=[1.5, 1.5])
    assert wide['pre'] == pytest.approx(-50) and np.isnan(wide['goo'])
    _, global_twice = _accuracy_summary(pints, all_in, [1.0], gaussian_variances=[1.0], global_variance=2.0)
    assert global_twice['pre'] == pytest.approx(50)
    _, no_gaussian = _accuracy_summary(pints, all_in, [1.0])
    assert np.isnan(no_gaussian['pre']) and np.isnan(no_gaussian['goo']) and no_gaussian['acc'] == 100


def test_normal_score_table():
    """Distinct values get the normal quantile of the midpoint of their rank; weights move the
    midpoints; tied values share the mean score; values beyond the data keep the end scores."""
    from scipy.stats import norm
    from pygeostat.statistics.utils import _normal_score_table
    values, scores = _normal_score_table([3.0, 1.0, 2.0, 4.0])
    np.testing.assert_array_equal(values, [1, 2, 3, 4])
    np.testing.assert_allclose(scores, norm.ppf([0.125, 0.375, 0.625, 0.875]))
    values, scores = _normal_score_table([1.0, 2.0], weights=[3.0, 1.0])
    np.testing.assert_allclose(scores, norm.ppf([0.375, 0.875]))
    values, scores = _normal_score_table([1.0, 2.0, 2.0, 3.0])
    np.testing.assert_allclose(scores[1], norm.ppf([0.375, 0.625]).mean())
    np.testing.assert_allclose(np.interp([0.0, 10.0], *_normal_score_table([1.0, 2.0])),
                               norm.ppf([0.25, 0.75]))


def test_accsim_precision_in_gaussian_units(accuracy_data):
    """P comes from the local variance of Gaussian realizations: given directly, or made by the
    normal score transform of the data; the true values play no part in it."""
    truth, reals = accuracy_data
    rng = np.random.default_rng(1)
    gaussian = rng.normal(0, 0.5, reals.shape)               # local variance about 0.25
    _, direct = gs.accsim(truth, reals, reals_gaussian=gaussian)
    assert direct['pre'] == pytest.approx(100 * (1 - np.mean(np.var(gaussian, axis=1))))
    assert 70 < direct['pre'] < 80
    _, shuffled = gs.accsim(rng.permutation(truth), reals, reals_gaussian=gaussian)
    assert shuffled['pre'] == pytest.approx(direct['pre'])
    _, transformed = gs.accsim(truth, reals, data=truth)     # transform made from the data
    assert transformed['pre'] < 100 and np.isfinite(transformed['goo'])
    _, without = gs.accsim(truth, reals)
    assert np.isnan(without['pre']) and np.isnan(without['goo']) and 0 <= without['acc'] <= 100
    with pytest.raises(ValueError, match='shape'):
        gs.accsim(truth, reals, reals_gaussian=gaussian[:, :5])


def test_accsim_default_probability_intervals(accuracy_data):
    """9 intervals (10% to 90%) for fewer than 300 true values and 19 (5% to 95%) otherwise, unless
    the increment is given."""
    truth, reals = accuracy_data
    assert len(truth) < 300
    propavg, _ = gs.accsim(truth, reals)
    np.testing.assert_allclose(propavg['ProbInt'], np.arange(1, 10) / 10)
    propavg, _ = gs.accsim(np.tile(truth, 3), np.tile(reals, (3, 1)))
    assert len(propavg) == (9 if len(truth) * 3 < 300 else 19)
    many = np.random.default_rng(0).normal(size=(300, 20))
    propavg, _ = gs.accsim(many[:, 0], many)
    assert len(propavg) == 19
    propavg, _ = gs.accsim(many[:, 0], many, pinc=0.25)
    assert len(propavg) == 3


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


def test_accsim_precision_with_a_different_number_of_realizations_per_location():
    """Locations may have different numbers of realizations, and the global variance must be positive."""
    rng = np.random.default_rng(0)
    truth = np.zeros(20)
    reals = np.empty(20, dtype=object)
    for i in range(20):
        reals[i] = rng.normal(size=5 + i % 3)
    _, stats = gs.accsim(truth, reals, data=rng.normal(size=50))
    assert np.isfinite(stats['pre']) and stats['pre'] < 100
    _, stats = gs.accsim(truth, reals, reals_gaussian=[r * 0 for r in reals])
    assert stats['pre'] == 100
    with pytest.raises(ValueError, match='positive'):
        gs.accsim(truth, reals, data=[1.0, 2.0], global_variance=0)
