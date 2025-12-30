(configuration)=
# Configuration & Parameters

PyGeoStat provides a powerful configuration system through the `Parameters` class that allows you to set default values for plotting styles, data handling, and workflow parameters. This eliminates the need to repeatedly specify the same arguments across function calls.

---

## Overview

The `Parameters` object is a dictionary-like object that stores default settings for:

- **Data parameters**: Coordinate columns, grid definitions, variable names
- **Plotting parameters**: Color maps, axis labels, figure sizes, fonts
- **Configuration options**: Verbosity, parallel processing, executable paths

All PyGeoStat functions check the `Parameters` object for default values. If you don't specify a parameter in a function call, PyGeoStat uses the value from `Parameters`.

:::{tip}
Function keyword arguments always override `Parameters` defaults for that specific call, without modifying the stored defaults.
:::

---

## Basic Usage

### Viewing Parameters

View all available parameters:

```python
import pygeostat as gs

# View all parameters
print(gs.Parameters)

# View specific parameter
print(gs.Parameters['plotting.cmap'])

# Search for parameters by pattern
plot_params = gs.Parameters.find_all('plotting')
print(plot_params)
```

### Setting Parameters

Modify parameters using dictionary-style assignment:

```python
import pygeostat as gs

# Set individual parameters
gs.Parameters['plotting.cmap'] = 'viridis'
gs.Parameters['plotting.unit'] = 'm'
gs.Parameters['plotting.xname'] = 'Easting'
gs.Parameters['plotting.yname'] = 'Northing'

# Update multiple parameters at once
gs.Parameters.update({
    'plotting.cmap': 'plasma',
    'plotting.unit': 'ft',
    'config.verbose': False
})
```

### Parameter Descriptions

Get help on any parameter:

```python
# Describe a specific parameter
gs.Parameters.describe('plotting.cmap')

# Describe all parameters (warning: very long output!)
gs.Parameters.describe()
```

---

## Common Configuration Scenarios

### Setting Up Spatial Coordinates

Configure default coordinate column names for your data:

```python
import pygeostat as gs

# Set coordinate column names
gs.Parameters['data.x'] = 'Easting'
gs.Parameters['data.y'] = 'Northing'
gs.Parameters['data.z'] = 'Elevation'

# Set drillhole ID column
gs.Parameters['data.dh'] = 'WellID'

# Now load data - PyGeoStat knows which columns are coordinates
data = gs.DataFile('mydata.csv')
```

### Configuring Grid Definitions

Set a default grid for your project:

```python
import pygeostat as gs

# Create a grid definition
griddef = gs.GridDef(
    '100 0.5 1.0, 100 0.5 1.0, 50 0.5 0.5'
)

# Set as default
gs.Parameters['data.griddef'] = griddef

# All functions that need a grid will now use this default
```

### Plotting Style Configuration

Customize the default appearance of all plots:

```python
import pygeostat as gs

# Color scheme
gs.Parameters['plotting.cmap'] = 'viridis'
gs.Parameters['plotting.cmap_cat'] = 'tab20'  # For categorical data

# Axis labels and units
gs.Parameters['plotting.unit'] = 'm'
gs.Parameters['plotting.xname'] = 'Easting'
gs.Parameters['plotting.yname'] = 'Northing'
gs.Parameters['plotting.zname'] = 'Elevation'
gs.Parameters['plotting.xabbrev'] = 'E'
gs.Parameters['plotting.yabbrev'] = 'N'
gs.Parameters['plotting.zabbrev'] = 'Elev'

# Histogram settings
gs.Parameters['plotting.histogram_plot.facecolor'] = '0.9'
gs.Parameters['plotting.histogram_plot.edgecolor'] = 'k'
gs.Parameters['plotting.histogram_plot.histbins'] = 20

# Location plot settings
gs.Parameters['plotting.location_plot.s'] = 30  # Point size
gs.Parameters['plotting.location_plot.c'] = '0.4'  # Default color
```

### Performance and System Configuration

```python
import pygeostat as gs

# Set number of parallel processes
gs.Parameters['config.nprocess'] = 8

# Control verbosity
gs.Parameters['config.verbose'] = True

# Use shipped executables
gs.Parameters['config.use_shipped_executable'] = True
```

---

## Saving and Loading Configurations

### Save Configuration to File

Save your current parameters to a file:

```python
import pygeostat as gs

# Configure your parameters
gs.Parameters['plotting.cmap'] = 'plasma'
gs.Parameters['plotting.unit'] = 'm'
gs.Parameters['data.x'] = 'Easting'

# Save to file
gs.Parameters.save('my_project_config.json')
```

### Load Configuration from File

Load parameters from a saved configuration:

```python
import pygeostat as gs

# Load parameters from file
gs.Parameters.load('my_project_config.json')

# Now all your saved settings are active
print(gs.Parameters['plotting.cmap'])  # 'plasma'
```

---

## System Defaults

PyGeoStat can store your preferred defaults in your home directory so they're automatically loaded every session.

### Set System Defaults

Save your current configuration as the system default:

```python
import pygeostat as gs

# Configure your preferences
gs.Parameters['plotting.cmap'] = 'viridis'
gs.Parameters['plotting.unit'] = 'm'
gs.Parameters['config.verbose'] = False

# Save as system defaults (stored in ~/.Pygeostat/Parameters.json)
gs.Parameters.set_systemdefault()
```

### Enable Auto-Loading

Enable automatic loading of system defaults:

```python
# Enable auto-load
gs.Parameters['config.autoload.parameters'] = True
gs.Parameters.set_systemdefault()

# Next time you import pygeostat, your defaults will load automatically
```

### Load System Defaults Manually

```python
import pygeostat as gs

# Load system defaults without auto-load
gs.Parameters.get_systemdefault()
```

### Reset to Factory Defaults

Restore all parameters to their original values:

```python
import pygeostat as gs

# Reset to pygeostat defaults
gs.Parameters.restore_defaults()

# Or reset system defaults file to factory settings
gs.Parameters.reset_systemdefault()
```

---

## Plot Style Management

PyGeoStat includes a `PlotStyle` class for managing matplotlib plotting styles. This is separate from `Parameters` and controls the visual appearance of matplotlib plots.

### Available Plot Styles

PyGeoStat provides several built-in styles:

- `'mpldefault'`: Matplotlib default style
- `'ccgpaper'`: Style for academic papers
- `'darkcontent'`: Dark theme for presentations
- `'presentation'`: Optimized for presentation slides

### Setting Plot Style

```python
import pygeostat as gs

# Set a plot style
gs.PlotStyle.set_style('ccgpaper')

# Use custom matplotlib rcParams
custom_style = {
    'figure.figsize': (10, 6),
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.grid': True,
    'grid.alpha': 0.3
}
gs.PlotStyle.set_style(custom_dict=custom_style)
```

### Saving and Loading Plot Styles

```python
import pygeostat as gs

# Save current plot style
gs.PlotStyle.save('my_plot_style.json')

# Load plot style from file
gs.PlotStyle.load('my_plot_style.json')

# Enable auto-loading of plot style
gs.Parameters['config.autoload.plot_style'] = True
gs.Parameters.set_systemdefault()
```

---

## Complete Example: Project Setup

Here's a complete example showing how to configure PyGeoStat for a mining project:

```python
import pygeostat as gs

# ===== Data Configuration =====
# Set coordinate columns
gs.Parameters.update({
    'data.x': 'East',
    'data.y': 'North',
    'data.z': 'RL',  # Relative Level
    'data.dh': 'HoleID'
})

# Define project grid
project_grid = gs.GridDef(
    '200 5.0 10.0, 200 5.0 10.0, 100 2.0 5.0'
)
gs.Parameters['data.griddef'] = project_grid

# ===== Plotting Configuration =====
# Spatial labels and units
gs.Parameters.update({
    'plotting.unit': 'm',
    'plotting.xname': 'Easting',
    'plotting.yname': 'Northing',
    'plotting.zname': 'RL',
    'plotting.xabbrev': 'E',
    'plotting.yabbrev': 'N',
    'plotting.zabbrev': 'RL'
})

# Color and style
gs.Parameters.update({
    'plotting.cmap': 'viridis',
    'plotting.cmap_cat': 'Set3',
    'plotting.grid': True,
    'plotting.sigfigs': 3
})

# Histogram defaults
gs.Parameters.update({
    'plotting.histogram_plot.histbins': 25,
    'plotting.histogram_plot.facecolor': '#e8e8e8',
    'plotting.histogram_plot.edgecolor': '#333333'
})

# ===== System Configuration =====
gs.Parameters.update({
    'config.nprocess': 12,  # Use 12 cores for parallel processing
    'config.verbose': True,
    'config.use_shipped_executable': True
})

# ===== Save Configuration =====
# Save project-specific config
gs.Parameters.save('mining_project_config.json')

# Optionally set as system default
# gs.Parameters['config.autoload.parameters'] = True
# gs.Parameters.set_systemdefault()

# ===== Set Plot Style =====
gs.PlotStyle.set_style('ccgpaper')

print("✓ Project configuration complete!")
```

Now you can use PyGeoStat without specifying these parameters in every function call:

```python
# Load data - uses configured coordinate columns
data = gs.DataFile('assay_data.csv')

# Create plots - uses configured colors, labels, units
gs.histogram_plot(data, var='Au')
gs.location_plot(data, var='Au')

# Uses configured grid
gs.slice_plot(data, var='Au', orient='xy', sliceno=50)
```

---

## Parameter Validation

PyGeoStat validates all parameter values:

```python
import pygeostat as gs

# Valid assignment
gs.Parameters['plotting.cmap'] = 'viridis'  # ✓ Valid colormap

# Invalid assignment - raises error
try:
    gs.Parameters['data.tmin'] = 'not_a_number'  # ✗ Expected float
except ValueError as e:
    print(f"Error: {e}")

# Invalid parameter name - shows warning
gs.Parameters['invalid.parameter'] = 'value'  # ⚠ Warning
```

---

## Available Parameter Categories

### Configuration Parameters (`config.*`)

- `config.verbose`: Print status messages
- `config.autoload.parameters`: Auto-load system defaults
- `config.autoload.plot_style`: Auto-load plot style
- `config.nprocess`: Number of parallel processes
- `config.use_shipped_executable`: Use bundled executables

### Data Parameters (`data.*`)

- `data.x`, `data.y`, `data.z`: Coordinate column names
- `data.dh`: Drillhole ID column
- `data.griddef`: Default grid definition
- `data.tmin`, `data.tmax`: Trimming limits

### Plotting Parameters (`plotting.*`)

- `plotting.unit`: Spatial unit label
- `plotting.xname`, `plotting.yname`, `plotting.zname`: Axis labels
- `plotting.cmap`: Default colormap
- `plotting.grid`: Show grid lines
- `plotting.sigfigs`: Significant figures for statistics

### Function-Specific Parameters

- `plotting.histogram_plot.*`: Histogram plot defaults
- `plotting.location_plot.*`: Location plot defaults
- `plotting.scatter_plot.*`: Scatter plot defaults

For a complete list of all parameters, use:

```python
import pygeostat as gs
for param in sorted(gs.Parameters.keys()):
    print(param)
```

---

## Tips and Best Practices

:::{tip}
**Project Configuration File**

Create a configuration file for each project and load it at the beginning of your scripts:

```python
import pygeostat as gs
gs.Parameters.load('project_config.json')
```
:::

:::{tip}
**System Defaults for Personal Preferences**

Use system defaults for your personal preferences (colors, fonts, etc.) and project-specific files for data-related settings (grid definitions, column names).
:::

:::{tip}
**Function Arguments Override Defaults**

Remember that explicit function arguments always take precedence over defaults:

```python
gs.Parameters['plotting.cmap'] = 'viridis'

# Uses 'plasma' for this call only, doesn't change the default
gs.histogram_plot(data, var='Au', cmap='plasma')

# Next call uses 'viridis' again
gs.histogram_plot(data, var='Cu')
```
:::

:::{warning}
**Parameter File Location**

System defaults are stored in `~/.Pygeostat/Parameters.json`. Don't edit this file directly - use `Parameters.set_systemdefault()` instead.
:::

---

## See Also

- {ref}`API Reference: Parameters <parameters>` - Complete API documentation
- {ref}`Getting Started <getting_started>` - Introduction to PyGeoStat
- {ref}`Plotting Functions <plotting>` - Using plotting parameters
