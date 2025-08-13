# TidalPy

<div style="text-align: center;">
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_win.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_win.yml/badge.svg?branch=main" alt="Windows Tests" /></a>
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_mac.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_mac.yml/badge.svg?branch=main" alt="MacOS Tests" /></a>
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_ubun.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests_ubun.yml/badge.svg?branch=main" alt="Ubuntu Tests" /></a><br />
    <!--<a href="https://codecov.io/github/jrenaud90/TidalPy" ><img src="https://codecov.io/github/jrenaud90/TidalPy/branch/main/graph/badge.svg?token=35OY4ZLOA5"/></a><br />-->
    <a href="https://doi.org/10.5281/zenodo.7017475"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7017475.svg" alt="DOI"></a>
    <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/Python-3.9|3.10|3.11|3.12|3.13-blue" alt="Python Version 3.9-3.13" /></a>
</div>

---

<a href="https://github.com/jrenaud90/TidalPy/releases"><img src="https://img.shields.io/badge/TidalPy-0.6.4 Alpha-orange" alt="TidalPy Version 0.6.5 Alpha" /></a>

**Tidal Dynamics and Thermal-Orbital Evolution Software Suite Implemented in Cython and Python**

TidalPy is an open-source software suite that utilizes a semi-analytic approach to estimate tidal heating,
orbit-rotation evolution, and thermal changes for rocky and icy worlds. It has been used to simulate the thermal-orbital
evolution of moons within our Solar System as well as exoplanets beyond. TidalPy's `RadialSolver` package can accurately
estimate the viscoelastic Love and Shida numbers for a multi-layered, rocky or icy world, including the effects of liquid layers,
compressibility, dynamic tides, and advanced rheological models. This module has been used to study the tidal response
of Mercury, Venus, Earth, our Moon, Mars, and much more.

Have any questions? Feel free to leave an [issue](https://github.com/jrenaud90/TidalPy/issues) or send a message to
[TidalPy@gmail.com](mailto:tidalpy@gmail.com).

# How to Install

## Compatibility

* **Windows-Latest**: *Installation & tests passed.*
* **MacOS-Latest**: *Installation & tests passed.*
* **Ubuntu-Latest**: *Installation & tests passed.*

## Basic Installation

Installing TidalPy is ass simple as ensuring 64-bit [Python 3.9+](https://www.python.org/) is installed on your 
system and performing the following in a terminal:

`pip install TidalPy`

Alternatively you can use conda to install TidalPy:

`conda install -c conda-forge TidalPy`

or 

`mamba install TidalPy`

## Accessing Jupyter Notebooks
There are several jupyter notebooks with TidalPy demos found in the /Demos/ folder of this repository.
In order to access these you will need to make sure you install Jupyter and a few related packages:

`pip install ipympl ipython ipywidgets jupyter`

or 

`conda install ipympl ipython ipywidgets jupyter`

You can then clone this GitHub repository,

`git clone https://github.com/jrenaud90/TidalPy`

to a local directory. Navigate to this directory and the `Demos` sub-directory then access the notebooks by using the command,
`jupyter notebook`.

## Cartopy

TidalPy offers the ability to make nice 2D plots using the [cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html) package for some of 
3d projection map plotting. In turn, cartopy relies on [GEOS](https://trac.osgeo.org/geos/) which is not a python
package and must be installed outside of pip.

The easiest way to install cartopy is using an Anaconda environment by,

`conda install -c conda-forge cartopy`

If you are not using a conda environment then you will need to find and install the GEOS binaries manually:

**Windows:** [Follow instructions here](https://trac.osgeo.org/osgeo4w/)
**On Ubuntu:** `sudo apt-get install libgeos-dev`
**On MacOS:** `brew install geos`

After GEOS is installed you can pip install the rest,

`pip install pyproj shapely pyshp cartopy`

## Installation Troubleshooting

_If you ran into a problem that is not listed below please [submit an issue](https://github.com/jrenaud90/TidalPy/issues) and we will work on addressing it!_

**Known Problems:**
* The `setuptools` package is required before TidalPy can be installed. Usually it is automatically installed, but if
  you are starting with a clean virtual environment it may not have been.
  * For Anaconda: `conda install setuptools`
  * Or for regular Python: `pip install setuptools`

# How to Use TidalPy

Check out the `Documentation\1 - Getting Started.md` file. This is pretty bare bones at the moment but offers some basic
info about TidalPy. For now the best way to learn how to use TidalPy is by checking out the `Demos` directory. There
are "beginner" [Jupyter notebooks](https://jupyter.org/) that are a great starting point.

## Using TidalPy for Science

TidalPy has been used in several studies already, and we encourage you to use it in yours. We would appreciate you
include a link back to this [page](https://github.com/jrenaud90/TidalPy) and cite one of the papers discussed in 
the next section. We also would love to see where TidalPy is being used! Please feel free to send us an
email: [TidalPy@gmail.com](mailto:TidalPy@gmail.com) when a paper or presentation utilized TidalPy. Anyone is welcome to
make forks or copies of TidalPy as long as their work references back to this page. License information can be found at
the end of this file.

## Citing TidalPy

If you use TidalPy for your research please cite its Zenodo [doi: 10.5281/zenodo.7017474](https://doi.org/10.5281/zenodo.7017474).

The science used in TidalPy is described in the following papers and software (with additional references therein):

* Rheological Modeling Package:
  * [Tidally Heated Terrestrial Exoplanets: Viscoelastic Response Models](https://ui.adsabs.harvard.edu/abs/2009ApJ...707.1000H/abstract)
  * [Increased Tidal Dissipation Using Advanced Rheological Models](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...98R/abstract)
* Non-synchronous Rotation Evolution and High Eccentricity Truncation Packages:
  * [Tidal Dissipation in Dual-Body, Highly Eccentric, and Non-synchronously Rotating System](https://ui.adsabs.harvard.edu/abs/2021PSJ.....2....4R/abstract)
  * [Tidal Evolution of the Keplerian Elements](https://ui.adsabs.harvard.edu/abs/2019CeMDA.131...30B/abstract)
* Third Party Software:
  * *Interior Model*: [BurnMan](https://github.com/geodynamics/burnman)
  * *Integration Routines*: [CyRK](https://zenodo.org/records/8329446)
  * *CVD Conscious Color Maps*: [Geodynamic Color Maps](http://doi.org/10.5281/zenodo.5501399)
  * *Projection Maps*: [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)
  * *Exoplanet data*: [Astroquery](https://github.com/astropy/astroquery/blob/main/astroquery/CITATION), [AstroPy](https://www.astropy.org/acknowledging.html)

## Contribute to TidalPy

TidalPy is in early alpha and there are lots of areas where it can improve! If you are interested in helping out, please
check out the information in `Documentation\Contribute.md`.

**Found a bug or have an idea for a new feature?**

* Go to TidalPy's [Github page](https://github.com/jrenaud90/TidalPy) and click the "Issues" tab then make a new report.
  * If you ran into a bug please include a code snippet (in markdown: code is designated by Grave accents surrounding
    the text) that reproduces the error (please keep this snippet as concise as possible).
  * It is helpful to triage issues when they are made. If you think you know the severity of a bug or can provide any
    other *at-a-glance* context, consider adding a "label" (right-hand side of the github issue form) to the issue.

# License Information
You are welcome to copy/fork TidalPy and make modifications assuming the following conditions are met:
* Links are included that point back to this [page](https://github.com/jrenaud90/TidalPy).
* Any software derived from TidalPy must remain open-source and non-commercial.

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view
a copy of this license,
visit [http://creativecommons.org/licenses/by-nc-sa/4.0/](http://creativecommons.org/licenses/by-nc-sa/4.0/) or send a
letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

# Acknowledgements
TidalPy was partially developed with support from NASA Goddard Space Flight Center's 
Sellers Exoplanet Environments Collaboration (SEEC) and Geodesy ISFM. 
TidalPy is partially based upon work supported by NASA under award number 80GSFC21M0002 and the
Center for Research and Exploration in Space Science & Technology II (CRESST II) administered at the University of
Maryland, College Park.

TidalPy has been improved by numerous contributors some of which you can find [here](https://github.com/jrenaud90/TidalPy/graphs/contributors).

TidalPy has benefited from work and conversations with the following:
- Wade G. Henning (U. of Maryland, College Park / NASA GSFC)
- Michael Efroimsky (U.S. Naval Observatory)
- Sander Goossens (NASA GSFC)
- Marc Neveu (U. of Maryland, College Park / NASA GSFC)
- Gael Cascioli (U. of Maryland, Baltimore County / NASA GSFC)
- Nick Wagner (Brown U.)
