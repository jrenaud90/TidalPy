# TidalPy
<p align="center">
  <img src="https://github.com/jrenaud90/TidalPy/tree/main/Documentation/_static/images/2025-11-28_Logo_2-4.svg" width="350" alt="TidalPy Logo">
</p>
<div style="text-align: center;">
    <a href="https://app.readthedocs.org/projects/tidalpy/builds/?version__slug=latest"><img src="https://app.readthedocs.org/projects/tidalpy/badge/?version=latest&style=flat" alt="TidalPy Documentation" /></a>
    <a href="https://doi.org/10.5281/zenodo.7017475"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7017475.svg" alt="DOI: 10.5281/zenodo.7017475"></a><br />
    <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/Python-3.9|3.10|3.11|3.12|3.13-blue" alt="Python Version 3.9-3.13" /></a>
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_win.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_win.yml/badge.svg?branch=main" alt="Windows Tests" /></a>
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_mac.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_mac.yml/badge.svg?branch=main" alt="MacOS Tests" /></a>
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_ubun.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_ubun.yml/badge.svg?branch=main" alt="Ubuntu Tests" /></a><br />
    <a href="https://pypi.org/project/TidalPy/"><img alt="PyPI - Downloads" src="https://img.shields.io/pypi/dm/TidalPy?label=PyPI%20Downloads" /></a>
    <a href="https://anaconda.org/conda-forge/tidalpy"> <img alt="Conda Downloads" src="https://img.shields.io/conda/d/conda-forge/tidalpy" /></a>
    <a href="https://emac.gsfc.nasa.gov?cid=2207-034"><img src="https://img.shields.io/badge/EMAC-2207--034-blue" alt="EMAC 2207-034"></a>
</div>

---

<a href="https://github.com/jrenaud90/TidalPy/releases"><img src="https://img.shields.io/badge/TidalPy-0.7.0a0.dev7-orange" alt="TidalPy Version 0.7.0a0.dev7" /></a>

**Tidal Dynamics and Thermal-Orbital Evolution Toolkit Implemented in Cython and Python**

TidalPy is an open source software suite that utilizes a semi-analytic approach to estimate tidal heating,
spin-orbit resonances, tidal & loading Love numbers, and thermal evolution for rocky and icy worlds. It has been used
to simulate the thermal-orbital evolution of moons within our Solar System as well as exoplanets beyond. TidalPy's
`RadialSolver` package can accurately estimate the viscoelastic Love and Shida numbers for a multi-layered, rocky or
icy world, including the effects of liquid layers, compressibility, dynamic tides, and advanced rheological models.
This module has been used to study the tidal response of Mercury, Venus, Earth, our Moon, Mars, and much more.

Have any questions or suggestions? Feel free to leave an [issue](https://github.com/jrenaud90/TidalPy/issues) or send
a message to [TidalPy@gmail.com](mailto:tidalpy@gmail.com).

## How to Install

### Compatibility

TidalPy has been developed to work on most modern operating systems. We specifically test it on the latest versions of
Ubuntu, Windows, and MacOS. We also pre-build binaries for these operating systems and provide them via
[PyPI](https://pypi.org/project/TidalPy/) or [Conda-Forge](https://anaconda.org/conda-forge/tidalpy). If a pre-built
binary is not available for your operating system version then see details about
[building TidalPy from source](https://tidalpy.readthedocs.io/en/latest/Overview/Readme.html#building-tidalpy-from-source).

* **Windows-Latest**: [![Windows Tests](https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_win.yml/badge.svg?branch=main)](https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_win.yml)
* **MacOS-Latest**: [![MacOS Tests](https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_mac.yml/badge.svg?branch=main)](https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_mac.yml)
* **Ubuntu-Latest**: [![Ubuntu Tests](https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_ubun.yml/badge.svg?branch=main)](https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_ubun.yml)

### Basic Installation

Installing TidalPy is as simple as ensuring 64-bit [Python 3.9+](https://www.python.org/) is installed on your 
system and running from a terminal:

`pip install TidalPy`

Alternatively you can use conda to install TidalPy:

`conda install -c conda-forge TidalPy`

or 

`mamba install TidalPy`

We recommend using a virtual environments (via a manager like [uv](https://docs.astral.sh/uv/pip/environments/) or
[miniforge](https://conda-forge.org/download/)) when installing TidalPy.

### Accessing Jupyter Notebooks
There are several demos provided with TidalPy that make use of [Jupyter notebooks](https://jupyter.org/) found in
the /Demos/ folder of TidalPy's [repository](https://github.com/jrenaud90/TidalPy). In order to access these you
will need to make sure you install Jupyter and a few related packages:

`pip install ipympl ipython ipywidgets jupyter`

or 

`conda install ipympl ipython ipywidgets jupyter`

You can then clone TidalPy's GitHub repository,

`git clone https://github.com/jrenaud90/TidalPy`

to a local directory. Navigate to this directory and the `Demos` sub-directory then access the notebooks by using the
command, `jupyter notebook`.

### Cartopy

TidalPy offers the ability to make 2D surface projection plots using the
[cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html) package. In turn, cartopy relies on 
[GEOS](https://trac.osgeo.org/geos/) which is not a python package and must be installed outside of pip.

#### Installing Cartopy using `conda`
The easiest way to install cartopy is using a conda environment by,

`conda install -c conda-forge cartopy`

#### Installing Cartopy using `pip`
If you are not using a conda environment then you will need to find and install the GEOS binaries manually:

**Windows:** [Follow instructions here](https://trac.osgeo.org/osgeo4w/)
**On Ubuntu:** `sudo apt-get install libgeos-dev`
**On MacOS:** `brew install geos`

After GEOS is installed you can pip install the rest,

`pip install pyproj shapely pyshp cartopy`

### Installation Troubleshooting

_If you ran into a problem that is not listed below please [submit an issue](https://github.com/jrenaud90/TidalPy/issues) and we will work on addressing it!_

**Known Problems:**
* The `setuptools` package is required before TidalPy can be installed. Usually it is automatically installed, but if
  you are starting with a clean virtual environment it may not have been.
  * For conda: `conda install setuptools`
  * Or pip: `pip install setuptools`

### Building TidalPy from Source
We automatically provide pre-built binaries for the latest version of MacOS, Ubuntu, and Windows via 
[PyPI](https://pypi.org/project/TidalPy/) or [Conda-Forge](https://anaconda.org/conda-forge/tidalpy). If your OS
version does not have pre-built binaries or if you are running into problems with the pre-builds, then you can build
TidalPy from its source code.

To do so, you will need to make sure that your environment has access to a C and C++ compiler that supports
C++20 standards, a recent version of Python, and has Cython 3.0+ installed. 

#### PyPI Build from Source
Using the source code uploaded to PyPI by running,

```bash
python -m pip install TidalPy -v --no-binary TidalPy
```

#### GitHub Repo
Alternatively you can clone the latest version of the GitHub repo and build locally,
```bash
git clone https://www.GitHub.com/jrenaud90/TidalPy.git
python -m pip install . -v --no-binary TidalPy  # The . assumes you have navigated to the directory with `pyproject.toml`
```

This is also the approach you would take to build TidalPy if you plan to edit its code. See more details about 
developing TidalPy [here](https://tidalpy.readthedocs.io/en/latest/Overview/Contributing.html).

#### Special consideration for MacOS
On MacOS, If you run into problems installing TidalPy then reinstall using the verbose flag (`pip install -v .`) to
look at the installation log. If you see an error that looks like "clang: error: unsupported option '-fopenmp'" then
you are likely using the default compiler or other compiler that does not support OpenMP. Read more about this issue
[here](https://github.com/facebookresearch/xformers/issues/157) and the steps taken
[here](https://github.com/jrenaud90/CyRK/blob/main/.github/workflows/push_tests_mac.yml). A fix for this issue is to
use `llvm`'s clang compiler. This can be done by doing the following in your terminal before installing TidalPy.

_Note this error can also occur when installing "CyRK" a critical dependency of TidalPy. The fix is the same as below,
just swap out "TidalPy" for "CyRK"._

```bash
brew install llvm
brew install libomp

# If on ARM64 (Apple Silicon) then do:
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
export CC=/opt/homebrew/opt/llvm/bin/clang
export CXX=/opt/homebrew/opt/llvm/bin/clang++

# Otherwise change these directories to:
export LDFLAGS="-L/usr/local/opt/llvm/lib"
export CPPFLAGS="-I/usr/local/opt/llvm/include"
export LDFLAGS="-L/usr/local/opt/libomp/lib"
export CPPFLAGS="-I/usr/local/opt/libomp/include"
export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++

pip install CyRK --no-binary="CyRK"
```

### TidalPy Versioning
TidalPy uses the major.minor.bugfix versioning scheme. In TidalPy's current state, we only promise to provide support
for the latest version (as found on the [GitHub release page](https://github.com/jrenaud90/TidalPy/releases)).
Therefore, older minor versions may not get critical patches (_e.g._, TidalPy 0.7.x will get patches where 0.6.x will
not even if the patches are applicable to that earlier version). If you would like to see a bug fix back ported to an
older version please open a [issue](https://github.com/jrenaud90/TidalPy/issues). 

In the future we hope to support multiple minor versions of TidalPy. At that time, all supported versions will be
listed in this section.


## Using TidalPy

Check out the [Getting Started](https://tidalpy.readthedocs.io/en/latest/Overview/1_Getting_Started.html) guide to
learn about TidalPy's features. The `Demos` directory is another good resource to learn by looking at
[Jupyter notebooks](https://jupyter.org/) that can teach you how to use TidalPy's features.

### Doing Science with TidalPy

TidalPy has been used in several studies already, and we encourage you to use it in yours. We would appreciate you
including a link back to this [page](https://github.com/jrenaud90/TidalPy) and cite one of the papers mentioned [here](https://tidalpy.readthedocs.io/en/latest/Readme.html#citing-tidalpy). We also would love to hear where TidalPy is being used! Please feel free to send us an
email: [TidalPy@gmail.com](mailto:TidalPy@gmail.com) when a paper or presentation utilized TidalPy. Anyone is welcome to
create forks or copies of TidalPy as long as their work references back to this page. License information can be found [here](https://tidalpy.readthedocs.io/en/latest/License.html).

### Citing TidalPy

If you use TidalPy for your research please cite the package by using the preferred citation found in the
[citation.cff](https://github.com/jrenaud90/TidalPy/blob/main/citation.cff). Currently, that is its Zenodo
[doi: 10.5281/zenodo.7017474](https://doi.org/10.5281/zenodo.7017474).

```bibtex
@software{2022zndo...7017475R,
       author = {{Renaud}, Joe P.},
        title = "{TidalPy}",
         year = 2022,
        month = jan,
          eid = {10.5281/zenodo.7017474},
          doi = {10.5281/zenodo.7017474},
    publisher = {Zenodo},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022zndo...7017475R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

It is good practice to cite the specific version of TidalPy you are using. Major versions have their own DOI on
[Zenodo](https://zenodo.org/records/16883555).

The science used in TidalPy is described in the following papers and software (with additional references therein):

* Rheology Module:
  * [Tidally Heated Terrestrial Exoplanets: Viscoelastic Response Models](https://ui.adsabs.harvard.edu/abs/2009ApJ...707.1000H/abstract)
  * [Increased Tidal Dissipation Using Advanced Rheological Models](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...98R/abstract)
* Dynamics and Tides Module:
  * [Tidal Dissipation in Dual-Body, Highly Eccentric, and Non-synchronously Rotating System](https://ui.adsabs.harvard.edu/abs/2021PSJ.....2....4R/abstract)
  * [Tidal Evolution of the Keplerian Elements](https://ui.adsabs.harvard.edu/abs/2019CeMDA.131...30B/abstract)
* RadialSolver (Love Number Calculator) Module:
  * [Constraining the Venus Interior Structure](https://ui.adsabs.harvard.edu/abs/2023PSJ.....4...65C).
  * [Seismic Surface Waves](https://www.doi.org/10.1016/B978-0-12-460811-5.50010-6)
* Third Party Software:
  * *Interior Model*: [BurnMan](https://github.com/geodynamics/burnman)
  * *Integration Routines*: [CyRK](https://zenodo.org/records/8329446)
  * *Graphics*: [Scientific Color Maps](https://doi.org/10.5281/zenodo.1243862), [Cartopy](https://scitools.org.uk/cartopy/docs/latest/), [Matplotlib](https://www.doi.org/10.1109/MCSE.2007.55)
  * *Exoplanet data*: [Astroquery](https://github.com/astropy/astroquery/blob/main/astroquery/CITATION), [AstroPy](https://www.astropy.org/acknowledging.html)
  * *Scientific Python*: [NumPy](https://doi.org/10.1038/s41586-020-2649-2), [SciPy](https://doi.org/10.1038/s41592-019-0686-2)
  * *Performance*: [Numba](https://doi.org/10.1145/2833157.2833162), [Cython](https://www.doi.org/10.1109/MCSE.2010.118)

## Contribute to TidalPy

TidalPy is in early alpha and there are lots of areas where it can improve! If you are interested in helping out, please
check out the information in `Documentation\Contribute.md`.

**Found a bug or have an idea for a new feature?**

* Go to TidalPy's [Github page](https://github.com/jrenaud90/TidalPy) and click the "Issues" tab then make a new report.
  * If you ran into a bug please include a code snippet (in markdown: code is designated by Grave accents surrounding
    the text) that reproduces the error (please keep this snippet as concise as possible).
  * It is helpful to triage issues when they are made. If you think you know the severity of a bug or can provide any
    other *at-a-glance* context, consider adding a "label" (right-hand side of the github issue form) to the issue.

# Acknowledgements
TidalPy was partially developed with support from NASA Goddard Space Flight Center's 
Sellers Exoplanet Environments Collaboration (SEEC) and Geodesy ISFMs. 
TidalPy is partially based upon work supported by NASA under award number 80GSFC21M0002 and the
Center for Research and Exploration in Space Science & Technology II (CRESST II) administered at the University of
Maryland, College Park.

TidalPy has been improved by numerous contributors some of which you can find [here](https://github.com/jrenaud90/TidalPy/graphs/contributors).

TidalPy has benefited from work and conversations with the following:
- Wade G. Henning (U. of Maryland, College Park / NASA GSFC)
- Michael Efroimsky (U.S. Naval Observatory)
- Michaela Walterová (Charles University)
- Sander Goossens (NASA GSFC)
- Marc Neveu (U. of Maryland, College Park / NASA GSFC)
- Gael Cascioli (U. of Maryland, Baltimore County / NASA GSFC)
- Nick Wagner (Brown University)

# License and Copyright
Copyright 2025 by [Joe P. Renaud](https://github.com/jrenaud90).
TidalPy is licensed under the Apache License, Version 2.0 (the "License"); you may not use this code except in
compliance with the License. You may obtain a copy of the License at
[www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0) or in this repository's LICENSE.md
file. Unless required by applicable law or agreed to in writing, software distributed under the License is distributed
on an "AS IS" basis, without warranties or conditions of any kind, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

You are welcome to copy/fork TidalPy and make modifications assuming the following conditions are met:
* Code repositories link back to TidalPy's original [repository](https://github.com/jrenaud90/TidalPy).
* Any published research cites this code using the preferred citation found in the
  [citation.cff file](https://github.com/jrenaud90/TidalPy/blob/main/citation.cff).

TidalPy's logo was originally designed by Ruhul Amin and modified by Joe P. Renaud.