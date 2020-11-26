#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

__author__ = 'pygeostat development team'
__date__ = '2020-01-04'
__version__ = '1.0.0'

import os, sys
try:
	import pygeostat as gs
except (ImportError, ModuleNotFoundError):
	sys.path.append(os.path.abspath(os.path.join(os.path.dirname( __file__ ), r'..')))
	import pygeostat as gs
import unittest
import warnings
import subprocess
import pandas as pd
import numpy as np
import tempfile



class BaseTest(unittest.TestCase):
	'''
	base class for testing plotting tools
	'''
	def setUp(self):
		warnings.simplefilter('ignore', category=ImportWarning)
		self.out_dir = tempfile.NamedTemporaryFile().name
		gs.mkdir(self.out_dir)

	def tearDown(self):
		gs.clrmplmem()
		gs.rmdir(self.out_dir)
		

class LocationPlotTest(BaseTest):
	'''
	 Test suite for pygeostat location_plot
	'''

	def setUp(self):
		super().setUp()
		self.data = gs.ExampleData('point3d_ind_mv')

	def test_simple_plot(self):
		'''
		Test the capability of plotting by passing only the data file
		'''
		
		output_file = os.path.join(self.out_dir, 'location_plot1.png')
		try:
			_ = gs.location_plot(self.data, output_file = output_file)
		except Exception as ex:
			self.fail('Unable to test the location plot with simple settings \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)

	def test_plot_cmap(self):
		'''
		Test the capability of plotting by passing variable name to get the color map
		'''
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'location_plot2.png')
		_, ax = plt.subplots(1,1)
		try:
			_ = gs.location_plot(self.data, var='Phi', output_file = output_file)
			_ = gs.location_plot(self.data, var='Phi', output_file = output_file, ax =ax)
		except Exception as ex:
			self.fail('Unable to test the location plot with passing the variable name \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_plot_lines(self):
		'''
		Test the capability of plotting along xz/yz orientation and use line plots based on drill hole id
		'''
		
		output_file = os.path.join(self.out_dir, 'location_plot3.png')
		try:
			_ = gs.location_plot(self.data, var='Phi', orient='yz', aspect =5, plot_collar = True, output_file = output_file)
		except Exception as ex:
			self.fail('Unable to test the location plot along yz/xz orientation \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_plot_categorical(self):
		'''
		Test the capability of plotting along xz/yz orientation and use line plots based on drill hole id for a categorical variable
		'''
		from matplotlib import pyplot as plt
		
		output_file = os.path.join(self.out_dir, 'location_plot3.png')
		_, ax = plt.subplots(1,1)
		try:
			_ = gs.location_plot(self.data, var='Lithofacies', orient='yz', aspect =5, plot_collar = True, output_file = output_file)
			_ = gs.location_plot(self.data, var='Lithofacies', orient='yz', aspect =5, plot_collar = True, output_file = output_file, ax = ax)
		except Exception as ex:
			self.fail('Unable to test the location plot along yz/xz orientation for a categorical variable \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class AccuracyPlotTest(BaseTest):
	'''
	 Test suite for pygeostat accuracy_plot
	'''

	def setUp(self):
		super().setUp()

	def test_simple_plot(self):
		'''
		Test the capability of plotting by passing only the data file
		'''
		from matplotlib import pyplot as plt

		data = gs.ExampleData('accuracy_plot')
		reals = data[list(data.columns[1:])].values
		truth = data[list(data.columns)[0]].values
		output_file = os.path.join(self.out_dir, 'accuracy_plot1.png')
		_, ax = plt.subplots(1,1)

		try:
			gs.accuracy_plot(truth=truth, reals=reals, output_file = output_file)
			gs.accuracy_plot(truth=truth, reals=reals, output_file = output_file, ax = ax)
		except Exception as ex:
			self.fail('Unable to test the accuracy plot with simple settings \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)
	

class SlicePlotTest(BaseTest):
	'''
	 Test suite for pygeostat slice_plot
	'''

	def setUp(self):
		super().setUp()
		self.griddef = gs.GridDef([40,0.5,1,40,0.5,1,40,0.5,1])
		self.data = gs.ExampleData('3d_grid', griddef=self.griddef)

	def test_simple_plot(self):
		'''
		Test the capability of plotting by passing only the data file
		'''
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'slice_plot1.png')
		_, ax = plt.subplots(1,1)

		try:
			gs.slice_plot(self.data, griddef=self.data.griddef, orient='xy', cmap='viridis', output_file=output_file)
			gs.slice_plot(self.data, griddef=self.data.griddef, orient='xy', cmap='viridis', output_file=output_file, ax =ax)
		except Exception as ex:
			self.fail('Unable to test the slice plot with simple settings \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_with_thickness_plot(self):
		'''
		Test the capability of plotting by passing only the data file
		'''
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'slice_plot2.png')
		_, ax = plt.subplots(1,1)

		try:
			gs.slice_plot(self.data, orient='xy', cmap='viridis', output_file=output_file, slice_thickness=20)
			gs.slice_plot(self.data, orient='xy', cmap='viridis', output_file=output_file, slice_thickness=20, ax =ax)
		except Exception as ex:
			self.fail('Unable to test the slice plot with simple settings \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class GridSlicePlotTest(SlicePlotTest):
	'''
	 Test suite for pygeostat gird_slice_plot
	'''

	def setUp(self):
		super().setUp()

	def test_simple_plot(self):
		'''
		Test the capability of plotting by passing only the data file
		'''
		
		output_file = os.path.join(self.out_dir, 'grid_slice_plot1.png')
		try:
			_ = gs.grid_slice_plot(self.data, output_file=output_file)
		except Exception as ex:
			self.fail('Unable to test the grid slice plot with simple settings \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_orient_slice_number_plot(self):
		'''
		Test the capability of specifying number of slices and orientation
		'''
		
		output_file = os.path.join(self.out_dir, 'grid_slice_plot2.png')
		try:
			_ = gs.grid_slice_plot(self.data, orient='xz', n_slice=5, output_file=output_file)
		except Exception as ex:
			self.fail('Unable to test the grid slice plot with number of slices and orientation \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)
	

	def test_nrow_ncol_plot(self):
		'''
		Test the capability of specifying number of rows or columns required for the slices
		'''
		
		output_file = os.path.join(self.out_dir, 'grid_slice_plot3.png')
		try:
			_ = gs.grid_slice_plot(self.data, orient='xz', n_slice=6, ncol=2, nrow=3, output_file=output_file)
		except Exception as ex:
			self.fail('Unable to test the grid slice plot number of rows or columns required for the slices \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_custom_plot(self):
		'''
		Test the ability to specify slice_plot kwargs using this function, so we can apply consistent
		custom formatting to all of the subplots:
		'''
		
		output_file = os.path.join(self.out_dir, 'grid_slice_plot4.png')
		try:
			_ = gs.grid_slice_plot(self.data, nrow=2, ncol=5, start_slice=10, end_slice=25, n_slice=10, cmap='hot', vlim=(-3,3), output_file=output_file)
		except Exception as ex:
			self.fail('Unable to test the grid slice plot to specify slice_plot kwargs \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class HistogramPlotTest(BaseTest):
	'''
	 Test suite for pygeostat Histogram Plot
	'''

	def setUp(self):
		super().setUp()
		self.data = gs.ExampleData("point3d_ind_mv")

	def test_simple_plot(self):
		'''
		Test the capability of plotting a simple histogram and saving the plot
		'''
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'histogram_plot1.png')
		_, ax = plt.subplots(1,1)

		try:
			gs.histogram_plot(self.data, var="Phi", bins=30, output_file = output_file)
			gs.histogram_plot(self.data, var="Phi", bins=30, output_file = output_file, ax = ax)
		except Exception as ex:
			self.fail('Unable to test the histogram plot with simple settings \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)

	def test_plot_density_with_color(self):

		'''
		Test the capability of plotting a histogram with density and custom colour and saving the plot
		'''
		from matplotlib import pyplot as plt
		output_file = os.path.join(self.out_dir, 'histogram_plot2.png')

		_, ax = plt.subplots(1,1)

		try:
			gs.histogram_plot(self.data, var="Phi", bins=30, color='#c2e1e5', sigfigs=5, log=True, density=True, output_file = output_file)
			gs.histogram_plot(self.data, var="Phi", bins=30, color='#c2e1e5', sigfigs=5, log=True, density=True, output_file = output_file, ax = ax)
		except Exception as ex:
			self.fail('Unable to test the histogram plot with density and custom color \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_plot_icdf(self):

		'''
		Test the capability of plotting a histogram fo idcf and custom stat location and saving the plot
		'''
		from matplotlib import pyplot as plt
		output_file = os.path.join(self.out_dir, 'histogram_plot3.png')
		_, ax = plt.subplots(1,1)

		try:
			gs.histogram_plot(self.data, var="Phi", icdf=True, color=3, lw=3.5, stat_xy=(1, 0.75), output_file = output_file)
			gs.histogram_plot(self.data, var="Phi", icdf=True, color=3, lw=3.5, stat_xy=(1, 0.75), output_file = output_file, ax = ax)
		except Exception as ex:
			self.fail('Unable to test the histogram plot for icdf and custom stat location \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_plot_var_per_category(self):

		'''
		Test the capability of plotting a histogram within different categorical variables
		'''

		output_file = os.path.join(self.out_dir, 'histogram_plot4.png')


		try:
			cats = [1, 2, 3, 4, 5]
			colors = gs.catcmapfromcontinuous("Spectral", 5).colors
			self.data.catdict = {c: "RT {:02d}".format(c) for c in cats}
			colordict =  {c: colors[i] for i, c in enumerate(cats)}

			gs.histogram_plot(self.data, var='Phi', cat=True, color=colordict, bins=40, figsize=(8, 4), xlabel=False, output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the histogram plot for a variable within different categories \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_plot_cdf_per_category(self):

		'''
		Test the capability of plotting a cdf within different categorical variables
		'''
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'histogram_plot5.png')


		try:
			cats = [1, 2, 3, 4, 5]
			colors = gs.catcmapfromcontinuous("Spectral", 5).colors
			self.data.catdict = {c: "RT {:02d}".format(c) for c in cats}
			colordict =  {c: colors[i] for i, c in enumerate(cats)}
			
			_, axs = plt.subplots(2, 2, figsize=(12, 9))
			axs=axs.flatten()
			
			for var, ax in zip(self.data.variables, axs):
				gs.histogram_plot(self.data, var=var, icdf=True, cat=True, color=colordict, ax=ax, output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the histogram plot for a variable within different categories \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_plot_proportions(self):

		'''
		Test the capability of plotting a histogram for categorical proportions
		'''
		from matplotlib import pyplot as plt
		output_file = os.path.join(self.out_dir, 'histogram_plot6.png')

		_, ax = plt.subplots(1,1)
		try:
			cats = [1, 2, 3, 4, 5]
			colors = gs.catcmapfromcontinuous("Spectral", 5).colors
			self.data.catdict = {c: "RT {:02d}".format(c) for c in cats}
			colordict =  {c: colors[i] for i, c in enumerate(cats)}

			gs.histogram_plot(self.data, cat=True, color=colordict, bins=40, figsize=(8, 4), xlabel=False, output_file=output_file)
			gs.histogram_plot(self.data, cat=True, color=colordict, bins=40, figsize=(8, 4), xlabel=False, output_file=output_file, ax = ax)

		except Exception as ex:
			self.fail('Unable to test the histogram plot for categorical proportions\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class HistogramPlotTestSimulation(BaseTest):
	
	def setUp(self):
		super().setUp()

		self.griddef = gs.GridDef([10,1,0.5, 10,1,0.5, 1,1,0.5])
		size = self.griddef.count()
		mean = 2.5
		scale = 3
		self.nreal = 100
		self.ref_data = pd.DataFrame({'value':np.random.normal(mean, scale, size = size)})

		for i, offset in enumerate(np.random.rand(self.nreal)*0.04 - 0.02):
			data = pd.DataFrame({'value':np.random.normal(mean+offset, scale, size = size)})
			gs.write_gslib(data, os.path.join(self.out_dir, 'real{}.out'.format(i+1)))

	def test_histogram_reproduction_plot(self):

		'''
		Test the capability of plotting histogram reproduction
		'''

		output_file = os.path.join(self.out_dir, 'histplt_sim_plot1.png')


		try:
			gs.histogram_plot_simulation(os.path.join(self.out_dir,'real*.out'), self.ref_data, reference_variable='value', 
								 		reference_color='C3', simulated_column=1, griddef = self.griddef, simulation_color='cyan',
										title='var', grid=True, nreal=self.nreal, output_file=output_file)
		except Exception as ex:
			self.fail('Unable to test the histogram reproduction plot \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)
			

class CorrelationMatrixPlotTest(BaseTest):
	'''
	 Test suite for pygeostat Correlation Matrix Plot
	'''

	def setUp(self):
		super().setUp()
		dfl = gs.ExampleData("point3d_ind_mv")
		data = dfl[dfl.variables]
		self.data_cor = data.corr() 


	def test_simple_plot(self):

		'''
		Test the capability of plotting a simple correlation matrix
		'''
		from matplotlib import pyplot as plt
		output_file = os.path.join(self.out_dir, 'correlation_matrix_plot1.png')

		_, ax = plt.subplots(1,1)
		try:
			gs.correlation_matrix_plot(self.data_cor, output_file=output_file)
			gs.correlation_matrix_plot(self.data_cor, output_file=output_file, ax =ax)
		except Exception as ex:
			self.fail('Unable to test the correlation matrix plot \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_simple_lower_matrix(self):

		'''
		Test the capability of plotting lower part of a simple correlation matrix 
		'''
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'correlation_matrix_plot2.png')

		_, ax = plt.subplots(1,1)

		try:
			gs.correlation_matrix_plot(self.data_cor, lower_matrix=True, annotation=True, output_file=output_file)
			gs.correlation_matrix_plot(self.data_cor, lower_matrix=True, annotation=True, output_file=output_file, ax = ax)
		except Exception as ex:
			self.fail('Unable to test the correlation matrix plot with only lower matrix option\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class ScatterPlotTest(BaseTest):

	def setUp(self):
		super().setUp()
		self.data = gs.ExampleData("point3d_ind_mv")

	def test_simple_plot(self):

		'''
		Test the capability of plotting a scatter plot with probability density estimate
		'''
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'scatter_plot1.png')

		_, ax = plt.subplots(1,1)

		try:
			gs.scatter_plot(self.data[self.data.variables[0]], self.data[self.data.variables[1]], cmap='hot', cbar=True, output_file = output_file)
			gs.scatter_plot(self.data[self.data.variables[0]], self.data[self.data.variables[1]], cmap='hot', cbar=True, output_file = output_file, ax =ax)

		except Exception as ex:
			self.fail('Unable to test the scatter plot \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class ScatterPlotsTest(ScatterPlotTest):

	def setUp(self):
		super().setUp()

	def test_scatter_plots(self):

		'''
		Test the capability of plotting a scatter plot with probability density estimate for multivariate data set
		'''

		output_file = os.path.join(self.out_dir, 'scatter_plots1.png')


		try:
			_ = gs.scatter_plots(self.data, nmax=1000, stat_xy=(0.95, 0.95), pad=(-1, -1), s=10, figsize=(10, 10), output_file = output_file)

		except Exception as ex:
			self.fail('Unable to test the scatter plots \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_scatter_plots_lu(self):

		'''
		Test the capability of plotting a scatter plot with probability density estimate for multivariate data set, creating an upper/lower matrix triangle of
		scatterplots for comparing the scatter of multiple variables in two data sets
		'''

		output_file = os.path.join(self.out_dir, 'scatter_plots2.png')


		try:
			_ = gs.scatter_plots_lu(self.data, self.data, titles=('Data', 'Realization'), s=10, nmax=1000,
							 stat_xy=(0.95, 0.95), pad=(-1, -1), figsize=(10, 10),
							 align_orient=True, output_file = output_file)

		except Exception as ex:
			self.fail('Unable to test the scatter plots creating an upper/lower matrix triangle \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class QqPlotTest(BaseTest):

	def setUp(self):
		super().setUp()

	def test_simple_plot(self):

		'''
		Test the capability of plotting a simple quantile by quantile comparison between two distributions
		'''
		import numpy as np
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'qq_plot1.png')
		_, ax = plt.subplots(1,1)

		try:
			_ = gs.qq_plot(np.random.randn(1000),np.random.randn(1000), npoints=500, output_file=output_file)
			_ = gs.qq_plot(np.random.randn(1000),np.random.randn(1000), npoints=500, output_file=output_file, ax =ax)

		except Exception as ex:
			self.fail('Unable to test the qq_plot \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)
		

class ProbabilityPlotTest(BaseTest):

	def setUp(self):
		super().setUp()

	def test_simple_plot(self):

		'''
		Test the capability of plotting a CDF with probability scale
		'''
		from matplotlib import pyplot as plt

		output_file = os.path.join(self.out_dir, 'probability_plot1.png')

		_, ax = plt.subplots(1,1)
		try:
			data = gs.ExampleData('oilsands')
			_ = gs.probability_plot(data.data['Fines'], logscale=False, output_file=output_file)
			_ = gs.probability_plot(data.data['Fines'], logscale=False, output_file=output_file, ax =ax)
		except Exception as ex:
			self.fail('Unable to test the probability_plot \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class ValidationPlotTest(BaseTest):

	def setUp(self):
		super().setUp()


	def test_simple_plot(self):
		'''
		Test the capability of plotting a simple scatter plot to compare estimate vs truth
		'''
		from matplotlib import pyplot as plt 

		output_file = os.path.join(self.out_dir, 'validation_plot1.png')
		_, ax = plt.subplots(1,1)

		try:
			data = gs.ExampleData('3d_estimate')
			_ = gs.validation_plot(data.data['Estimate'], data.data['True'], stat_blk='minimal', output_file = output_file)
			_ = gs.validation_plot(data.data['Estimate'], data.data['True'], stat_blk='minimal', output_file = output_file, ax =ax)

		except Exception as ex:
			self.fail('Unable to test the validation_plot with simple settings\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)

	def test_detailed_plot(self):
		'''
		Test the capability of plotting a simple scatter plot to compare estimate vs truth
		'''
		from matplotlib import pyplot as plt 

		output_file = os.path.join(self.out_dir, 'validation_plot2.png')
		_, ax =plt.subplots(1,1)

		try:
			mean = [0, 0]
			cov = [[1, 0.95], [0.95, 1]]  # diagonal covariance
			x, y = np.random.multivariate_normal(mean, cov, 5000).T
			gs.validation_plot(x,y,vlim=(-3.5, 3.5) ,grid=True, stat_xy=(1, 0.68), output_file=output_file)
			gs.validation_plot(x,y,vlim=(-3.5, 3.5) ,grid=True, stat_xy=(1, 0.68), output_file=output_file, ax =ax)

		except Exception as ex:
			self.fail('Unable to test the validation_plot \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class LoadingsPlotTest(BaseTest):

	def setUp(self):
		super().setUp()
		refdata = gs.ExampleData('3d_correlation')
		self.loadmat = refdata.data.corr().iloc[3:6,6:9]

	def test_simple_plot1(self):
		'''
		Test the capability of plotting a simple loadings plot for original variables and the principal components (using pandas dataframe)
		'''

		output_file = os.path.join(self.out_dir, 'loadings_plot1.png')

		try:
			_ = gs.loadings_plot(self.loadmat, output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the loadings_plot with simple settings (using pandas dataframe)\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_simple_plot2(self):
		'''
		Test the capability of plotting a simple loadings plot for original variables and the principal components (using a numpy array)
		'''

		output_file = os.path.join(self.out_dir, 'loadings_plot2.png')

		try:
			_ = gs.loadings_plot(self.loadmat.values, output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the loadings_plot with simple settings (using a numpy array)\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_with_labels(self):
		'''
		Test the capability of plotting a loadings plot with detailed settings for original variables and the principal components
		'''

		output_file = os.path.join(self.out_dir, 'loadings_plot3.png')

		try:
			_ = gs.loadings_plot(self.loadmat.values, figsize=(4,4), xticklabels=['PC1', 'PC2', 'PC3'], yticklabels=['InputVariable1', 'InputVariable2', 'InputVariable3'], output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the loadings_plot with  tick labels and figsize\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_with_tickrotation(self):
		'''
		Test the capability of plotting a loadings plot with detailed settings for original variables and the principal components
		'''

		output_file = os.path.join(self.out_dir, 'loadings_plot4.png')

		try:
			_ = gs.loadings_plot(self.loadmat.values, figsize=(5,5), xticklabels=['PC1', 'PC2', 'PC3'],
				 yticklabels=['InputVariable1', 'InputVariable2', 'InputVariable3'],
				 rotateticks=(False, True), output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the loadings_plot with tick rotation settings\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class DrillPlotTest(BaseTest):

	def setUp(self):
		super().setUp()
		dat = gs.ExampleData('point3d_ind_mv')
		self.data = dat.data[dat.data['HoleID'] == 3]

	def test_continuous_plot(self):
		'''
		Test the capability of plotting a continuous variable
		'''

		output_file = os.path.join(self.out_dir, 'drill_plot1.png')

		try:
			_ = gs.drill_plot(self.data['Elevation'], self.data['Sw'], grid = True, output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the drill_plot for a continuous variable\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_categorical_plot(self):
		'''
		Test the capability of plotting a categorical variable
		'''

		output_file = os.path.join(self.out_dir, 'drill_plot2.png')

		try:
			_ = gs.drill_plot(self.data['Elevation'], self.data['Lithofacies'], output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the drill_plot for a categorical variable\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


	def test_categorical_custom_plot(self):
		'''
		Test the capability of plotting a categorical variable with custom categorical dictionary
		'''

		output_file = os.path.join(self.out_dir, 'drill_plot3.png')

		try:
			cat_dict = {1: {'name': 'Sand', 'color': 'gold'},
			3: {'name': 'SHIS','color': 'orange'},
			4: {'name': 'MHIS','color': 'green'},
			5: {'name': 'MUD','color': 'gray'}}

			_ = gs.drill_plot(self.data['Elevation'], self.data['Lithofacies'], categorical_dictionary=cat_dict, output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the drill_plot for a categorical variable with custom categorical dictionary\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class VariogramPlotTest(BaseTest):

	def setUp(self):
		super().setUp()

	def test_experimental_plot(self):
		'''
		Test the capability of plotting experimental variograms
		'''

		output_file = os.path.join(self.out_dir, 'variogram_plot1.png')

		try:
			dfl = gs.ExampleData('experimental_variogram')
			_ = gs.variogram_plot(dfl, experimental=True, figsize=(8,4), output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the variogram_plot for experimental variogram\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)

	def test_model_plot(self):
		'''
		Test the capability of plotting variogram models
		'''

		output_file = os.path.join(self.out_dir, 'variogram_plot2.png')

		try:
			dfl = gs.ExampleData('variogram_model')
			_ = gs.variogram_plot(dfl, experimental=False, lw = 3, figsize=(8,4), output_file=output_file)

		except Exception as ex:
			self.fail('Unable to test the variogram_plot for experimental variogram\n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)


class ContourPlotTest(BaseTest):
	'''
	 Test suite for pygeostat countour_plot
	'''

	def setUp(self):
		super().setUp()
		grid_str = '''120 5.00  10.00 -nx, xmn, xsiz
					110 1205.00  10.00 -ny, ymn, ysiz
					1 0.5 1.0  -nz, zmn, zsiz'''
		
		griddef = gs.GridDef(grid_str=grid_str)
		self.data = gs.ExampleData("grid2d_surf", griddef=griddef)

	def test_simple_plot(self):
		'''
		Test the capability of plotting by passing only the essential parameters
		'''
		

		output_file = os.path.join(self.out_dir, 'countour_plot1.png')
		try:
			_ = gs.contour_plot(self.data, var="Thickness", clabel=True, output_file=output_file)
		except Exception as ex:
			self.fail('Unable to test the countour plot with simple settings \n{}'.format(str(ex)))

		self.assertEqual(os.path.isfile(output_file), True)
	
	def test_gmm_plot(self):
		'''
		Test the GMM plotting
		'''
		try:
			gmm_util = gs.GmmUtility(gmm_file=os.path.join(os.path.normpath(os.getcwd() + os.sep + os.pardir), r'pygeostat\pygeostat\data\example_data\gmm_fit.out'), 
                      data=gs.ExampleData('point2d_mv').data, variable_names=['Var1', 'Var2','Var3'])
		except Exception as ex:
			self.fail('Unable to test GmmUtility \n{}'.format(str(ex)))
		try:
			gmm_util.bivariate_plot(var_index=[1,2], cmap='viridis',title='Bivariate Plot')
		except Exception as ex:
			self.fail('Unable to test bivariate_plot \n{}'.format(str(ex)))
		try:
			gmm_util.summary_plot(pad=0.1)
		except Exception as ex:
			self.fail('Unable to test summary_plot \n{}'.format(str(ex)))
		try:
			gmm_util.univariate_conditional_plot(conditioning_data=[0, 0,None])
		except Exception as ex:
			self.fail('Unable to test univariate_conditional_plot \n{}'.format(str(ex)))

if __name__ == '__main__':
	subprocess.call([sys.executable, '-m', 'unittest', str(__file__), '-v'])
