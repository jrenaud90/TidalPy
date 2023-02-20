<h1 align="center">TidalPy</h1>
<h3 align="center">v0.4.0 Beta</h3>

<p align="center">
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/pr_tests_win.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/pr_tests_win.yml/badge.svg?branch=main" alt="Windows Tests" /></a>
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/pr_tests_mac.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/pr_tests_mac.yml/badge.svg?branch=main" alt="MacOS Tests" /></a>
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/pr_tests_ubun.yml"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/pr_tests_ubun.yml/badge.svg?branch=main" alt="Ubuntu Tests" /></a>
    <a href="https://codecov.io/github/jrenaud90/TidalPy" ><img src="https://codecov.io/github/jrenaud90/TidalPy/branch/main/graph/badge.svg?token=35OY4ZLOA5"/></a><br />
    <a href="https://mybinder.org/v2/gh/jrenaud90/TidalPy/master?filepath=%2FDemos%2F"><img src="https://mybinder.org/badge_logo.svg" alt="Binder" /></a>
    <a href="https://doi.org/10.5281/zenodo.7017475"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7017475.svg" alt="DOI"></a>
</p>

## Purpose

TidalPy is an open-source software suite designed to assist researchers in the semi-analytic calculation of tidal
dissipation and subsequent orbit-spin evolution for rocky and icy worlds.

**TidalPy is intended to be a...**

* Black Box (in the documentation this is referred to as the "*OOP* scheme" for Object-Oriented Programming)
    * TidalPy serves as simple to install (cross-platform) and, hopefully, simple to use package that users can pick up
      and hit the ground running.
    * The OOP scheme performs many calculations with very little input from the user. The major drawbacks are
      performance (in some situations) and that many assumptions are opaque to the user without some digging.
* Tool Box (referred to as the "*Functional* scheme")
    * TidalPy also contains many efficient functions to perform calculations relevant to tides and thermal-orbital
      coupling. These can be quickly imported and used in a custom scripts.
        * In general, the functional scheme will have much higher performance, flexibility, and extensibility than OOP.
          It also generally makes assumptions more visible to the user. The downside is the user may need to be more
          familiar with the underlying physics.

*Once you are comfortable with TidalPy, it is usually a good idea to mix the two schemes: take some aspects of OOP that
you don't want to deal with and build on them with some of TidalPy's or your own functions.*

### Limitations

The major limitations of the current version of TidalPy are...

* A multilayer model has now been implemented, but it is not currently part of the OOP scheme.
* Chemical and phase changes within a planet's layers have not been implemented.

### Related Software

Below is a non-exhaustive list of publicly available software that performs similar or parallel calculations as TidalPy.

* Are you interested in the habitability of a planet? With considerations of tides, atmospheres, water content, solar
  interactions? Check out...
    * [VPLanet](https://github.com/VirtualPlanetaryLaboratory/vplanet)
* Are you interested in the orbital evolution of multiple planets with each planet influencing one another? Consider an
  N-body approach like...
    * [Posidonius (formerly Mercury-T)](https://github.com/marblestation/posidonius)
    * [ReboundX](https://github.com/dtamayo/reboundx)
* Don't care about tides or orbital dynamics but are instead interested in interior structure and composition of
  planets?
    * [BurnMan](https://github.com/geodynamics/burnman)
    * [PerpleX](http://www.perplex.ethz.ch/)
* Are you interested in tides, interiors, *and* the chemical evolution of small worlds but don't care about
  non-synchronous rotation or compressibility of planets?
    * [IcyDwarf](https://github.com/MarcNeveu/IcyDwarf)

However, if you want high fidelity tidal, orbital, spin, and interior models --- then you have come to the right place!
Read below for instructions on how to install and use TidalPy.

## How to Install

### Compatibility

*As of TidalPy v0.4.0*:

* **Windows-Latest**: *Installation & tests passed.*
* **MacOS-Latest**: *Installation & tests passed.*
* **Ubuntu-Latest**: *Installation & tests passed.*

### Simple Installation

As simple as ensuring 64-bit [Python 3.8+](https://www.python.org/) is installed on your system and performing the
following in a terminal:

`pip install TidalPy`

However, there can be several gotchas that come with this simple installation process. It is recommended to use the
advanced installation described in the next section.

_Note: As of TidalPy v0.3.4, if you do not install via anaconda then you will need to manually install 
[proj v8.0.0+](https://proj.org/install.html) and [geos 3.7.2+](https://anaconda.org/conda-forge/geos) 
(before installing TidalPy. Otherwise, the project map graphic utility may not work. 
[Read more here](https://scitools.org.uk/cartopy/docs/latest/installing.html)_ 

### Accessing Jupyter Notebooks
There are several jupyter notebooks with TidalPy demos found in the /Demos/ and /Derivative/ directories.
In order to access these you will need to make sure you install Jupyter and a few related packages:

`pip install ipympl ipython ipywidgets jupyter`

or 

`conda install ipympl ipython ipywidgets jupyter`

Then you can navigate to these directories in a terminal and access the notebooks by using the command,
`jupyter notebook`.

### Advanced Installation

It is highly recommended that you use the [Anaconda](https://www.anaconda.com/distribution/) distribution of Python.
This has pre-compiled binaries for several packages that TidalPy uses and will generally negate a lot of potential
headaches. It is also recommended that you use a virtual environment. Using Anaconda, a new virtual environment can be
made with `conda create -n <name> python=3.9` and switched to with `conda activate <name>`.

* Get the latest version of TidalPy from Github.
    * Ensure you have the latest version of [git](https://git-scm.com/downloads)
      or [github](https://desktop.github.com/). Clone the TidalPy git
      using `git clone https://github.com/jrenaud90/TidalPy.git`.
        * Whenever you want to update TidalPy simply navigate to this directory and use `git pull`. Since TidalPy is in
          early development, it is recommended you check for updates regularly. Updates will **not**
          download automatically.

* Using a terminal, navigate to the TidalPy directory that contains `setup.py` and then:
    * For Anaconda Python:
        * Run `conda install --file conda_requirements.txt -c defaults -c conda-forge; pip install -e .` *(That trailing
          period is important, don't leave it out!)*
    * For non-Anaconda Python:
        * Run `pip install -e .` *(That trailing period is important, don't leave it out!)*

* Test your installation:
    * Navigate to the TidalPy directory that contains `setup.py` in a terminal.
    * Ensure you have `pytest` package installed (`conda install pytest` or `pip install pytest`).
    * Run pytest by simply using the command `pytest` from your terminal:
        * Running all the tests can take a while (currently around 10 minutes), if all you are interested in is checking that
          TidalPy installed correctly then you can let pytest check the first dozen tests if they are passing then you
          can quit the test suite early.
        * If no errors show up (warnings are okay and expected) then you should hopefully be good to go.
    * Open a new terminal *not in the TidalPy directory* (e.g., your desktop).
        * Run `python` and then try to `import TidalPy`; if that works try the command `TidalPy.version` if you do not
          get any import errors, and the version number is as you expect, then TidalPy was successfully installed.

### Additional Dependencies

TidalPy will only install 3rd-party packages that are absolutely needed to run the majority of 
its functionality. However, there are some functions throughout TidalPy that take advantage of
additional packages. You can install these additional packages by calling (note you need to install
cartopy and julia's dependencies first --- see the next few sections),

`pip install -r additional_dependencies.txt`

This command should be run from an elevated terminal to avoid permission issues.

#### Cartopy

TidalPy utilizes the [cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html) package for some of 
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

#### Diffeqpy / Julia

TidalPy provides the option to use the [Julia](https://julialang.org/) programing language's differential equation 
solver for python: [diffeqpy](https://github.com/SciML/diffeqpy). To utilize this package you first need to ensure
that Julia is installed on your machine and available via the system's environment path.

* Install the Julia language from [https://julialang.org/downloads/](https://julialang.org/downloads/)
* Add Julia's directory and its `bin` subdirectory to your system's path.
* Open an elevated ("as administrator") terminal, command prompt, or powershell.
    * If you are using a virtual Python environment make sure it is active.
* Install `julia` and `diffeqpy` for python using pip
    * Run `pip install julia diffeqpy`
* Open Python on your elevated terminal (the following steps may take a while to compile). 
  * Run `import julia; julia.install(); import diffeqpy; diffeqpy.install()`

### Installation Troubleshooting

* The `setuptools` package is required before TidalPy can be installed. Usually it is automatically installed, but if
  you are starting with a clean virtual environment it may not have been.
    * For Anaconda: `conda install setuptools`
    * Or for regular Python: `pip install setuptools`
* The current version of TidalPy is in Alpha and will receive many updates on a relatively fast schedule. So, it is
  recommended that you run it from an IDE and/or install it as
  an [editable package](https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs). If you do not wish to
  install as an editable package then please remove all `-e` flags.

## How to Use

Check out the `Documentation\Getting Started.md` file. This is pretty bare bones at the moment but offers some basic
info about TidalPy. For now the best way to learn how to use TidalPy is by checking out the `Demos` directory. There
are "beginner" [Jupyter notebooks](https://jupyter.org/) that are a great starting point.

## Using TidalPy for Science

TidalPy has been used in several studies already, and we encourage you to use it in yours. We would appreciate you
include a link back to this [page](https://github.com/jrenaud90/TidalPy) and cite one of the papers below (if you
utilized a specific package). We also would love to see where TidalPy is being used! Please feel free to send us an
email: [TidalPy@gmail.com](mailto:TidalPy@gmail.com) when a paper or presentation utilized TidalPy. Anyone is welcome to
make forks or copies of TidalPy as long as their work references back to this page. License information can be found at
the end of this file.

### Referencing TidalPy

The science used in TidalPy is described in the following papers (and references therein):

* Rheological Modeling Package:
    * [Tidally Heated Terrestrial Exoplanets: Viscoelastic Response Models](https://ui.adsabs.harvard.edu/abs/2009ApJ...707.1000H/abstract)
    * [Increased Tidal Dissipation Using Advanced Rheological Models](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...98R/abstract)
* Non-synchronous Rotation Evolution and High Eccentricity Truncation Packages:
    * [Tidal Dissipation in Dual-Body, Highly Eccentric, and Non-synchronously Rotating System](https://ui.adsabs.harvard.edu/abs/2021PSJ.....2....4R/abstract)
    * [Tidal Evolution of the Keplerian Elements](https://ui.adsabs.harvard.edu/abs/2019CeMDA.131...30B/abstract)
* Third Party Software:
    * *Interior Model*: [BurnMan](https://github.com/geodynamics/burnman)
    * *Integration Routines*: [diffeqpy](https://github.com/SciML/diffeqpy), [Julia DiffEq](https://diffeq.sciml.ai/v2.0/)
    * *CVD Conscious Color Maps*: [Geodynamic Color Maps](http://doi.org/10.5281/zenodo.5501399)
    * *Projection Maps*: [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)

## Contribute

TidalPy is in early alpha and there are lots of areas where it can improve! If you are interested in helping out, please
check out the information in `Documentation\Contribute.md`.

**Found a bug or have an idea for a new feature?**

* Go to TidalPy's [Github page](https://github.com/jrenaud90/TidalPy) and click the "Issues" tab then make a new report.
    * If you ran into a bug please include a code snippet (in markdown: code is designated by Grave accents surrounding
      the text) that reproduces the error (please keep this snippet as concise as possible).
    * It is helpful to triage issues when they are made. If you think you know the severity of a bug or can provide any
      other *at-a-glance* context, consider adding a "label" (right-hand side of the github issue form) to the issue.

## License Information
You are welcome to copy/fork TidalPy and make modifications assuming the following conditions are met:
* Links are included that point back to this [page](https://github.com/jrenaud90/TidalPy).
* Any software derived from TidalPy must remain open-source and non-commercial.

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view
a copy of this license,
visit [http://creativecommons.org/licenses/by-nc-sa/4.0/](http://creativecommons.org/licenses/by-nc-sa/4.0/) or send a
letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
