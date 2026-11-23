# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- 3D viewer (optional, `pip install pygeostat[view3d]`): `gs.view3d()`, `gs.Viewer`, `gs.drillhole_mesh()`
  and `gs.block_mesh()` for drill holes, block models and sub-blocks, with a user guide page (#147)
- Interactive layer panel for the 3D viewer (optional, `pip install "pygeostat[view3d-panel]"`):
  `Viewer.panel()`, with a slicer (a plane through the scene), and `Viewer.set_layer()`, `set_slicer()`,
  `set_axis_origin()`, `get_view()` and `set_view()` to do the same from a script (#149)
- Display options and picking for the 3D layer panel: colormap, colour by another variable, show only some
  categories, cell edges, tube radius, vertical exaggeration, background, scene boxes and ruler ticks, camera
  views and fit view, picking a drill hole, and screenshot and HTML export; from a script with
  `Viewer.set_background()`, `set_scene()`, `set_axis_ticks()`, `fit_view()`, `set_parallel_projection()`,
  `set_camera_view()`, `interval_info()`, `save_screenshot()` and `export_html()`, new arguments of
  `Viewer.set_layer()` and the `other_vars` argument of `gs.view3d()` (#151)
- Shared plot options in `pygeostat.plotting.utils` (`_format_ticks`, `_format_fontsize`, `_format_axes`) and the parameter `gs.Parameters['plotting.ticks']` (default `True`), so that plots offer `ticks`, `fontsize`, `xlim`, `ylim`, `grid` and `background` with the same meaning
- `histogram_plot()` takes `ticks`, `fontsize` and `background`, next to its existing `xlim`, `ylim` and `grid`
- `validation_plot()` takes `outliers` and `outlier_color` to draw the points furthest from the 45 degree line in colour with their index, `outlier_labels` and `reg_line` to leave out those numbers and the regression line, and the shared options `ticks`, `fontsize`, `background`
- `histogram_plot_simulation()` shows the number of values per realization (`'realn'` in `stat_blk`) and takes `ylim`, `grid`, `ticks`, `fontsize` and `background`; `variogram_plot_simulation()` takes `grid`, `ticks`, `fontsize` and `background`, so the reproduction plots have the same options as the other plots
- `histogram_plot_simulation()` takes `lw_real`, `ylabel` and `stat_fontsize`, and `variogram_plot_simulation()` takes `sill_color`, `sill_lw` and `real_label`, for values that were fixed in the code. The defaults are unchanged
- `accuracy_plot()` takes `color` and `label` to compare models on one axis, `reals_gaussian`, `data`, `weights` and `global_variance` for the precision and goodness

### Fixed
- `variogram_plot_simulation()` applies `ls`, which was accepted but ignored
- Accuracy statistics for MIK (`accmik`): it no longer crashes (missing imports, pandas 3 indexing), the lower tail is padded with 0 instead of 1, the upper tail is padded so true values above the highest threshold are counted, and the average variance is exact instead of a random estimate. A `truth` with several columns is now an error

### Changed
- `accuracy_plot()` follows the updated definitions of accuracy, precision and goodness of Deutsch (2026), Paper 40, CCG Annual Report 28, and the statistics block labels them A26, P26 and G26. It also applies `stat_fontsize` and `axis_xy` (accepted but ignored before), shows the statistics as decimals, and draws coloured deviation lines from each point to the 1:1 line
- The statistics block of `validation_plot()` shows the mean squared error (`MSE`) where it showed the root mean squared error (`RMSE`)
- Plots that use the shared options draw tick marks by default (4 points, pointing outwards). The plot styles had no major tick marks before; set `gs.Parameters['plotting.ticks'] = False` or pass `ticks=False` to get the old look
- The realizations of the histogram and variogram reproduction plots are drawn in a lighter grey by default (`simclr` `'0.5'`, `alpha` `0.3`)

## [1.2.0] - 2025-06-25

### Major Updates
- **Complete documentation overhaul** with modern Sphinx theme
- **Python version support** expanded to include 3.10, 3.11, 3.12, and 3.13
- **Dependency modernization** with updated versions of all major dependencies

### Added
- Documentation with pydata-sphinx-theme
- Modern badge-based README with clear installation instructions
- Support for latest Python versions (3.10-3.13)
- Modern build system using pyproject.toml

### Testing
- **Migrated test suite to pytest** from unittest for better testing experience
- Modern test fixtures and parametrized testing
- Automated testing across multiple Python versions in CI

### Changed
- **Breaking**: Minimum Python version raised to 3.10
- Documentation completely rebuilt with modern structure and examples
- Build system migrated from setup.py to pyproject.toml
- Updated package structure following modern Python standards
- Enhanced example data and tutorials

### Updated Dependencies
- matplotlib >= 3.10.0 (was < 3.0)
- numpy >= 2.1.3 (major version update)
- pandas >= 2.2.3 (major version update)
- scipy >= 1.15.2 (major version update)
- scikit-learn >= 1.6.1 (major version update)
- All other dependencies updated to latest stable versions

## [1.1.1] - 2021-11-13

### Previous Release
- Last stable release before major modernization
- Support for Python 3.6-3.9
- Legacy documentation and build system
- Original CCG geostatistics functionality

---

## Migration Guide

### Upgrading from 1.1.1 to 1.2.1

**Python Version Requirements:**
- **Minimum Python version is now 3.10** (was 3.6)
- Ensure you're using Python 3.10 or newer before upgrading

**Installation:**
```bash
# Uninstall old version
pip uninstall pygeostat

# Install new version
pip install pygeostat

# Verify installation
python -c "import pygeostat; print(pygeostat.__version__)"
```

**Code Compatibility:**
- Most existing code should work without changes
- Some deprecated matplotlib/pandas patterns may need updating
- Check the new documentation for updated examples

