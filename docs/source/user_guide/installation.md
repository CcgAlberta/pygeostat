(installation)=
# Installation

pygeostat is compatible with **Python 3.10 or newer** and follows modern Python packaging standards.

:::{important}
**Minimum Python version:** Python 3.10

For historical releases supporting older Python versions, see pygeostat `v1.1.1`.
:::

---

## Requirements

Before installing pygeostat, ensure that you have:

- Python **3.10+**
- `pip` (or a compatible Python package installer)

Generic Python installation instructions are intentionally omitted here.
If you need help installing Python, see the [official Python documentation](https://www.python.org/downloads/).

---

## Quick install (recommended)

pygeostat is published on the Python Package Index (PyPI).

Install it with:

```shell
pip install pygeostat
```

This installs the core library and all required Python dependencies.

---

## Virtual environments (recommended)

Using a virtual environment is strongly recommended to isolate dependencies.

### Using `venv` (standard Python)

```shell
python -m venv pygeostat-env
source pygeostat-env/bin/activate   # Linux / macOS
# pygeostat-env\Scripts\activate    # Windows

pip install pygeostat
```

### Using `conda` (alternative)

1. Create the environment

```shell
conda create -n <environment_name> python=3.10
```

2. Activate the environment

```shell
conda activate <environment_name>
```

3. Install pygeostat in the activated environment

```shell
pip install pygeostat
```

---

## Additional Software ([CCG](https://CCGAlberta.com/) Members)

[CCG](https://CCGAlberta.com/) members have the option to install CCG/GSLIB software (executable files) to enable pygeostat {ref}`scripting features <programs>`. This can be done using the following function. While GSLIB executable files are available through a public repository, for [CCG](https://CCGAlberta.com/) software a valid access token is required to download executables files from a private repository. The access token is available for CCG members at CCG knowledge base website.

```{eval-rst}
.. autofunction:: pygeostat.utility.get_executable
   :no-index:
```
