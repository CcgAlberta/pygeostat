
.. _statistics:
.. public

Statistics
##########

Collection of tools for calculating statistics.



CDF
===

Cumulative Distribution Function
--------------------------------
.. autofunction:: pygeostat.statistics.cdf.cdf

Percentile from CDF
-------------------------
.. autofunction:: pygeostat.statistics.cdf.percentile_from_cdf

z_percentile
------------
.. autofunction:: pygeostat.statistics.cdf.z_percentile

Build Indicator CDF
-------------------
.. autofunction:: pygeostat.statistics.cdf.build_indicator_cdf

Variance From CDF
-----------------
.. autofunction:: pygeostat.statistics.cdf.variance_from_cdf

Stdev From CDF
--------------
.. autofunction:: pygeostat.statistics.cdf.stdev_from_cdf


Kernel Density Estimation Functions
===================================

Univariate KDE with StatsModels
-------------------------------
.. autofunction:: pygeostat.statistics.kde.kde_statsmodels_u

Multivariate KDE with StatsModels
---------------------------------
.. autofunction:: pygeostat.statistics.kde.kde_statsmodels_m

KDE with Scikit-learn
---------------------
.. autofunction:: pygeostat.statistics.kde.kde_sklearn

KDE with Scipy
--------------
.. autofunction:: pygeostat.statistics.kde.kde_scipy




Weighted Statistics
===================

Weighted Mean
-------------
.. autofunction:: pygeostat.statistics.utils.weighted_mean

Weighted Variance
-----------------
.. autofunction:: pygeostat.statistics.utils.weighted_variance

Weighted Skewness
-----------------
.. autofunction:: pygeostat.statistics.utils.weighted_skew


Weighted Kurtosis
-----------------
.. autofunction:: pygeostat.statistics.utils.weighted_kurtosis

Weighted Correlation
--------------------
.. autofunction:: pygeostat.statistics.utils.weighted_correlation

Weighted Rank Correlation
-------------------------
.. autofunction:: pygeostat.statistics.utils.weighted_correlation_rank

Weighted Covariance
-------------------
.. autofunction:: pygeostat.statistics.utils.weighted_covariance



Assorted Stats Functions
========================

Nearest Positive Definite Correlation Matrix
--------------------------------------------
.. autofunction:: pygeostat.statistics.utils.near_positive_definite

Accuracy Plot Statistics - Simulation
-------------------------------------
.. autofunction:: pygeostat.statistics.utils.accsim

Accuracy Plot Statistics - CDF thresholds
-----------------------------------------
.. autofunction:: pygeostat.statistics.utils.accmik


PostSim
=======
.. autofunction:: pygeostat.statistics.postsim.postsim_multfiles
