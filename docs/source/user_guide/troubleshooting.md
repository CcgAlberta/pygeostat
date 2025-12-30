(troubleshooting)=
# Troubleshooting & FAQ

This guide helps you resolve common issues when installing and using PyGeoStat.

---

## Installation Issues

### Python Version Error

**Problem:** Installation fails with Python version error

```
ERROR: Package 'pygeostat' requires Python >=3.10 but you are using Python 3.9
```

**Solution:**

Upgrade to Python 3.10 or higher. PyGeoStat v1.2.0+ requires Python 3.10 minimum.

```bash
# Check your Python version
python --version

# If using conda, create environment with Python 3.10+
conda create -n pygeostat_env python=3.10
conda activate pygeostat_env
pip install pygeostat

# If using system Python, install Python 3.10+
# Then use pip install pygeostat
```

:::{seealso}
See the {ref}`Installation Guide <installation>` for detailed setup instructions.
:::

---

### Dependency Conflicts

**Problem:** pip reports dependency conflicts during installation

```
ERROR: pip's dependency resolver does not currently take into account all the packages that are installed...
```

**Solution:**

Use a fresh virtual environment:

```bash
# Create a new virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install pygeostat
pip install pygeostat
```

For conda users:

```bash
# Create fresh conda environment
conda create -n pygeostat_clean python=3.10
conda activate pygeostat_clean
pip install pygeostat
```

---

### Import Errors After Installation

**Problem:** PyGeoStat imports fail even though installation succeeded

```python
>>> import pygeostat
ModuleNotFoundError: No module named 'pygeostat'
```

**Solution:**

1. **Verify you're in the correct environment:**

```bash
# Check which Python is active
which python  # On Windows: where python

# Check if pygeostat is installed
pip list | grep pygeostat
```

2. **Reinstall in correct environment:**

```bash
pip uninstall pygeostat
pip install pygeostat
```

3. **Check for multiple Python installations:**

```bash
# Use python -m to ensure correct environment
python -m pip install pygeostat
```

---

## Data Loading Issues

### File Not Found Errors

**Problem:** DataFile cannot find your file

```python
data = gs.DataFile('mydata.csv')
# FileNotFoundError: [Errno 2] No such file or directory: 'mydata.csv'
```

**Solution:**

1. **Use absolute paths:**

```python
import os
import pygeostat as gs

# Get absolute path
file_path = os.path.abspath('mydata.csv')
print(f"Looking for file at: {file_path}")

data = gs.DataFile(file_path)
```

2. **Check current working directory:**

```python
import os
print(f"Current directory: {os.getcwd()}")
print(f"Files here: {os.listdir('.')}")
```

---

### Column Not Found Errors

**Problem:** PyGeoStat can't find coordinate columns

```python
gs.location_plot(data, var='Grade')
# KeyError: 'x'
```

**Solution:**

Set coordinate column names in Parameters or specify explicitly:

```python
import pygeostat as gs

# Option 1: Set default coordinate columns
gs.Parameters['data.x'] = 'Easting'
gs.Parameters['data.y'] = 'Northing'
gs.Parameters['data.z'] = 'Elevation'

# Now plotting works
gs.location_plot(data, var='Grade')

# Option 2: Specify columns directly
gs.location_plot(data, var='Grade', x='Easting', y='Northing')
```

:::{seealso}
See the {ref}`Configuration Guide <configuration>` for managing default parameters.
:::

---

### Grid Definition Errors

**Problem:** Grid definition format is invalid

```python
grid = gs.GridDef('invalid')
# AssertionError: Grid definition string is invalid
```

**Solution:**

Use correct GSLIB grid format: `nx dx xmn, ny dy ymn, nz dz zmn`

```python
import pygeostat as gs

# Correct format: count spacing origin
grid = gs.GridDef('100 1.0 0.0, 100 1.0 0.0, 50 1.0 0.0')

# Or use keyword arguments
grid = gs.GridDef(nx=100, dx=1.0, xmn=0.0,
                  ny=100, dy=1.0, ymn=0.0,
                  nz=50, dz=1.0, zmn=0.0)
```

---

## Plotting Issues

### Figures Not Displaying

**Problem:** Plotting functions run but no figure appears

**Solution:**

1. **Enable interactive backend:**

```python
import matplotlib.pyplot as plt
plt.ion()  # Turn on interactive mode

import pygeostat as gs
gs.histogram_plot(data, var='Grade')
plt.show()  # Force display
```

2. **Use Jupyter/IPython magic:**

```python
# In Jupyter notebooks
%matplotlib inline

import pygeostat as gs
gs.histogram_plot(data, var='Grade')
```

3. **Save figure explicitly:**

```python
import pygeostat as gs

fig, ax = gs.histogram_plot(data, var='Grade')
fig.savefig('histogram.png', dpi=300, bbox_inches='tight')
```

---

### Colormap Errors

**Problem:** Invalid colormap name

```python
gs.Parameters['plotting.cmap'] = 'invalid_cmap'
# ValueError: 'invalid_cmap' is not a valid colormap name
```

**Solution:**

Use valid matplotlib colormaps:

```python
import pygeostat as gs

# Valid sequential colormaps
gs.Parameters['plotting.cmap'] = 'viridis'  # Default
# Other options: 'plasma', 'inferno', 'magma', 'cividis'

# Diverging colormaps
gs.Parameters['plotting.cmap'] = 'RdBu'
# Other options: 'seismic', 'coolwarm', 'bwr'

# Categorical colormaps
gs.Parameters['plotting.cmap_cat'] = 'tab20'
# Other options: 'Set1', 'Set2', 'Set3', 'Paired'

# List all available colormaps
import matplotlib.pyplot as plt
print(plt.colormaps())
```

---

### Plot Size Issues

**Problem:** Plots are too small or too large

**Solution:**

Control figure size through PlotStyle or function arguments:

```python
import pygeostat as gs

# Method 1: Set default figure size via PlotStyle
gs.PlotStyle.set_style(custom_dict={'figure.figsize': (12, 8)})

# Method 2: Specify figsize per plot
fig, ax = gs.histogram_plot(data, var='Grade', figsize=(10, 6))

# Method 3: Adjust existing figure
import matplotlib.pyplot as plt
fig = plt.gcf()
fig.set_size_inches(12, 8)
```

---

## CCG Software Issues

### Executable Not Found

**Problem:** PyGeoStat can't find CCG software executables

```python
gs.Program(...)
# FileNotFoundError: CCG executable not found
```

**Solution:**

1. **Check if executables are in PATH:**

```bash
# On Linux/Mac
which varmap
which nscore

# On Windows
where varmap.exe
where nscore.exe
```

2. **Add CCG software to PATH:**

```bash
# Linux/Mac (add to ~/.bashrc or ~/.zshrc)
export PATH="/path/to/ccg/bin:$PATH"

# Windows: Add to System Environment Variables
# Control Panel → System → Advanced → Environment Variables
```

3. **Configure executable path in Parameters:**

```python
import pygeostat as gs

# Specify executable directory
gs.Parameters['config.executable_dir'] = '/path/to/ccg/bin'

# Or use shipped executables (if available)
gs.Parameters['config.use_shipped_executable'] = True
```

:::{note}
CCG software executables are **not included** with PyGeoStat. They must be obtained separately from the [CCG](http://www.ccgalberta.com).
:::

---

## Performance Issues

### Slow Data Loading

**Problem:** Loading large CSV files is slow

**Solution:**

1. **Use efficient data formats:**

```python
import pygeostat as gs

# Instead of CSV, use Parquet (much faster for large files)
data.to_parquet('mydata.parquet')
data = gs.DataFile('mydata.parquet')

# Or HDF5
data.to_hdf('mydata.h5', key='data')
data = gs.DataFile('mydata.h5', key='data')
```

2. **Load only needed columns:**

```python
import pandas as pd
import pygeostat as gs

# Load only specific columns
cols = ['X', 'Y', 'Z', 'Grade', 'Density']
df = pd.read_csv('large_file.csv', usecols=cols)
data = gs.DataFile(data=df)
```

---

### Memory Errors with Large Datasets

**Problem:** Out of memory errors when working with large grids

```python
# MemoryError: Unable to allocate array
```

**Solution:**

1. **Process data in chunks:**

```python
import pygeostat as gs

# Process subsets of data
for i, chunk in enumerate(pd.read_csv('large_file.csv', chunksize=10000)):
    data_chunk = gs.DataFile(data=chunk)
    # Process chunk
    result = process_data(data_chunk)
    # Save partial results
    result.to_csv(f'result_chunk_{i}.csv')
```

2. **Use data types efficiently:**

```python
import numpy as np
import pygeostat as gs

# Use appropriate dtypes to reduce memory
data['Grade'] = data['Grade'].astype(np.float32)  # Instead of float64
data['Category'] = data['Category'].astype('category')
```

3. **Increase parallel processing:**

```python
import pygeostat as gs

# Use more cores (if available)
gs.Parameters['config.nprocess'] = 16
```

---

## Common Error Messages

### "Key 'x' not found in Parameters"

**Cause:** Trying to use a function that requires coordinate columns, but they're not set.

**Fix:**

```python
import pygeostat as gs

gs.Parameters['data.x'] = 'Easting'
gs.Parameters['data.y'] = 'Northing'
gs.Parameters['data.z'] = 'Elevation'
```

---

### "GridDef object required"

**Cause:** Function needs a grid definition but none is provided.

**Fix:**

```python
import pygeostat as gs

# Option 1: Set default grid
grid = gs.GridDef('100 1.0 0.0, 100 1.0 0.0, 50 1.0 0.0')
gs.Parameters['data.griddef'] = grid

# Option 2: Pass grid to function
my_function(data, griddef=grid)
```

---

### "ValueError: Could not convert to float"

**Cause:** Trying to set a numeric parameter with a string value.

**Fix:**

```python
import pygeostat as gs

# Wrong
gs.Parameters['plotting.gammasize'] = 'large'  # ✗

# Correct
gs.Parameters['plotting.gammasize'] = 3.0  # ✓
```

---

## Frequently Asked Questions

### How do I update PyGeoStat?

```bash
pip install --upgrade pygeostat
```

Check current version:

```python
import pygeostat as gs
print(gs.__version__)
```

---

### Can I use PyGeoStat with Python 3.9?

No. PyGeoStat v1.2.0+ requires **Python 3.10 or higher**. This is due to dependencies on modern NumPy, pandas, and matplotlib versions.

If you must use Python 3.9, install the legacy version:

```bash
pip install pygeostat==1.1.1
```

:::{warning}
Version 1.1.1 is **no longer maintained** and lacks recent features and bug fixes.
:::

---

### Where are my system defaults stored?

System defaults are stored in:

- **Linux/Mac:** `~/.Pygeostat/Parameters.json`
- **Windows:** `C:\Users\<username>\.Pygeostat\Parameters.json`

View the location:

```python
import os
config_file = os.path.join(os.path.expanduser('~'), '.Pygeostat', 'Parameters.json')
print(f"Config file: {config_file}")
```

---

### How do I reset to factory defaults?

```python
import pygeostat as gs

# Reset current session
gs.Parameters.restore_defaults()

# Reset system defaults file
gs.Parameters.reset_systemdefault()
```

---

### Can I use PyGeoStat without CCG software?

**Yes!** PyGeoStat provides extensive functionality without CCG executables:

- Data loading and manipulation
- Statistical analysis
- Visualization and plotting
- Grid operations
- Data transformations

CCG executables are **only required** for:
- Running GSLIB-style programs
- Specific geostatistical modeling workflows

---

### How do I use PyGeoStat in a Jupyter notebook?

```python
# Enable inline plotting
%matplotlib inline

import pygeostat as gs
import matplotlib.pyplot as plt

# Load and plot data
data = gs.ExampleData('oilsands')
gs.histogram_plot(data, var='Bitumen')
plt.show()
```

---

### Why are my plots not showing colors?

Check your colormap settings:

```python
import pygeostat as gs

# View current colormap
print(gs.Parameters['plotting.cmap'])

# Set a valid colormap
gs.Parameters['plotting.cmap'] = 'viridis'

# For categorical data
gs.Parameters['plotting.cmap_cat'] = 'tab20'
```

---

### How do I save my configuration?

```python
import pygeostat as gs

# Configure parameters
gs.Parameters['plotting.cmap'] = 'plasma'
gs.Parameters['data.x'] = 'Easting'

# Save to project file
gs.Parameters.save('my_project.json')

# Load later
gs.Parameters.load('my_project.json')
```

:::{seealso}
See the {ref}`Configuration Guide <configuration>` for more details.
:::

---

## Getting Help

If you encounter issues not covered here:

1. **Check the documentation:**
   - {ref}`Installation Guide <installation>`
   - {ref}`Getting Started <getting_started>`
   - {ref}`Configuration Guide <configuration>`
   - {ref}`API Reference <api>`

2. **Search existing issues:**
   Visit the [GitHub Issues](https://github.com/CcgAlberta/pygeostat/issues) page

3. **Report a bug:**
   Create a new issue with:
   - Python version (`python --version`)
   - PyGeoStat version (`import pygeostat; print(pygeostat.__version__)`)
   - Operating system
   - Minimal reproducible example
   - Full error message and traceback

4. **Get community help:**
   - Visit [CCG Website](http://www.ccgalberta.com)
   - Check [Geostatistics Lessons](http://geostatisticslessons.com/)

---

## Debugging Tips

### Enable Verbose Mode

See more diagnostic information:

```python
import pygeostat as gs

gs.Parameters['config.verbose'] = True
```

### Check Loaded Data

Inspect your data after loading:

```python
import pygeostat as gs

data = gs.DataFile('mydata.csv')

# Basic info
print(f"Shape: {data.shape}")
print(f"Columns: {list(data.columns)}")
print(f"Data types:\n{data.dtypes}")

# Check for missing values
print(f"Missing values:\n{data.isnull().sum()}")

# Summary statistics
print(data.describe())
```

### Test with Example Data

Verify PyGeoStat works correctly:

```python
import pygeostat as gs

# Load example data
data = gs.ExampleData('oilsands')
print(f"Example data loaded: {data.shape}")

# Test plotting
gs.histogram_plot(data, var='Bitumen')
print("✓ Plotting works")

# Test statistics
print(data.describe())
print("✓ Statistics work")
```

If example data works but your data doesn't, the issue is likely with your data format or configuration.
