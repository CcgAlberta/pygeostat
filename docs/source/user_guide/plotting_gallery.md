(plotting_gallery)=
# Plotting Gallery

PyGeoStat provides a comprehensive library of plotting functions for geostatistical data analysis and visualization. This gallery showcases the main plotting capabilities.

:::{seealso}
For detailed documentation of each function, see the {ref}`Plotting API Reference <plotting>`.
:::

---

## Statistical Plots

::::{grid} 3
:gutter: 3

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.histogram_plot
:img-top: ../figures/PlottingGallery/histogram_plot1.png

**histogram_plot**

Standard histogram with optional CDF overlay and summary statistics.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.histogram_plot
:img-top: ../figures/PlottingGallery/histogram_plot2.png

**histogram_plot (weighted)**

Histogram with weighted data and custom styling.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.histogram_plot
:img-top: ../figures/PlottingGallery/histogram_plot3.png

**histogram_plot (log scale)**

Logarithmic histogram for skewed distributions.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.histogram_plot_simulation
:img-top: ../figures/PlottingGallery/histogram_plot_simulation.png

**histogram_plot_simulation**

Compare multiple realizations against reference data.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.probability_plot
:img-top: ../figures/PlottingGallery/probability_plot.png

**probability_plot**

Probability plot for distribution analysis.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.qq_plot
:img-top: ../figures/PlottingGallery/qq_plot.png

**qq_plot**

Quantile-quantile plot for comparing distributions.
:::

::::

---

## Scatter & Correlation Plots

::::{grid} 3
:gutter: 3

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.scatter_plot
:img-top: ../figures/PlottingGallery/scatter_plot.png

**scatter_plot**

Basic scatter plot with correlation statistics.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.scatter_plots
:img-top: ../figures/PlottingGallery/scatter_plots.png

**scatter_plots**

Matrix of scatter plots for multiple variables.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.scatter_plots_lu
:img-top: ../figures/PlottingGallery/scatter_plots_lu.png

**scatter_plots_lu**

Lower-upper triangle scatter plot matrix.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.correlation_matrix_plot
:img-top: ../figures/PlottingGallery/Correlation_plot.png

**correlation_matrix_plot**

Correlation matrix heatmap visualization.
:::

::::

---

## Spatial Plots

::::{grid} 3
:gutter: 3

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.location_plot
:img-top: ../figures/PlottingGallery/location_plot.png

**location_plot**

Plot spatial data at sample locations.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.slice_plot
:img-top: ../figures/PlottingGallery/slice_plot.png

**slice_plot**

2D slice through 3D spatial data.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.grid_slice_plot
:img-top: ../figures/PlottingGallery/grid_slice_plot.png

**grid_slice_plot**

Gridded data slice visualization.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.contour_plot
:img-top: ../figures/PlottingGallery/contour_plot.png

**contour_plot**

Contour plot for spatial data.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.drill_plot
:img-top: ../figures/PlottingGallery/drill_plot.png

**drill_plot**

Drillhole log visualization.
:::

::::

---

## Geostatistical Plots

::::{grid} 3
:gutter: 3

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.variogram_plot
:img-top: ../figures/PlottingGallery/variogram_plot.png

**variogram_plot**

Experimental and model variogram display.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.accuracy_plot
:img-top: ../figures/PlottingGallery/accuracy_plot.png

**accuracy_plot**

Accuracy assessment for estimates.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.validation_plot
:img-top: ../figures/PlottingGallery/validation_plot.png

**validation_plot**

Cross-validation results visualization.
:::

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.pit_plot
:img-top: ../figures/PlottingGallery/pitplot_mr.png

**pit_plot**

Probability integral transform (PIT) plot.
:::

::::

---

## Multivariate Analysis

::::{grid} 3
:gutter: 3

:::{grid-item-card}
:link: plotting.html#pygeostat.plotting.loadings_plot
:img-top: ../figures/PlottingGallery/loadings_plot.png

**loadings_plot**

Factor loadings visualization for PCA/FA.
:::

::::

---

## Getting Started with Plotting

### Basic Example

```python
import pygeostat as gs

# Load example data
data = gs.ExampleData('oilsands')

# Create a histogram
gs.histogram_plot(data, var='Bitumen')

# Create a location plot
gs.location_plot(data, var='Bitumen',
                 x='Easting', y='Northing')

# Create a scatter plot
gs.scatter_plot(data, var1='Bitumen', var2='Fines')
```

### Customizing Plots

Configure default plotting styles using Parameters:

```python
import pygeostat as gs

# Set global plotting parameters
gs.Parameters['plotting.cmap'] = 'viridis'
gs.Parameters['plotting.unit'] = 'm'
gs.Parameters['plotting.xname'] = 'Easting'
gs.Parameters['plotting.yname'] = 'Northing'

# All plots will now use these defaults
gs.location_plot(data, var='Grade')
```

### Saving Figures

```python
import pygeostat as gs

# Create plot and get figure handle
fig, ax = gs.histogram_plot(data, var='Grade')

# Save to file
fig.savefig('histogram.png', dpi=300, bbox_inches='tight')
```

---

## Next Steps

::::{grid} 2
:gutter: 3

:::{grid-item-card} 📚 Plotting API Reference
:link: plotting
:link-type: ref

Complete documentation of all plotting functions with parameters and examples.
:::

:::{grid-item-card} ⚙️ Configuration Guide
:link: configuration
:link-type: ref

Learn how to customize default plotting styles and parameters.
:::

:::{grid-item-card} 🚀 Getting Started
:link: getting_started
:link-type: ref

Introduction to PyGeoStat workflows and basic usage.
:::

:::{grid-item-card} 🧪 Examples
:link: examples
:link-type: ref

Complete tutorials demonstrating PyGeoStat in real-world workflows.
:::

::::
