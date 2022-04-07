<h1 style="text-align: center">TidalPy</h1>
<h3 style="text-align: center">v0.3.5 Beta</h3>

<p style="text-align: center">
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests.yml?query=branch%3Amaster"><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/push_tests.yml/badge.svg?branch=master" alt="Push Test Pass/Fail" /></a>
    <a href="https://codecov.io/gh/jrenaud90/TidalPy"><img src="https://codecov.io/gh/jrenaud90/TidalPy/branch/master/graph/badge.svg?token=35OY4ZLOA5" alt="Code Coverage"/></a><br />
    <a href="https://mybinder.org/v2/gh/jrenaud90/TidalPy/master?filepath=%2FCookbooks%2F"><img src="https://mybinder.org/badge_logo.svg" alt="Binder" /></a>
    <a href="https://gitter.im/TidalPy/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge"><img src="https://badges.gitter.im/TidalPy/community.svg" alt="Gitter" /></a>
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

*As of TidalPy v0.3.0a*:

* **Win10**: *Installation & tests passed.*
* **MacOS (Catalina)**: *Installation & tests passed.*
* **CentOS7**: *TBD*

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

### Installation Troubleshooting

* The `setuptools` package is required before TidalPy can be installed. Usually it is automatically installed, but if
  you are starting with a clean virtual environment it may not have been.
    * For Anaconda: `conda install setuptools`
    * Or for regular Python: `pip install setuptools`
* The current version of TidalPy is in Alpha and will receive many updates on a relatively fast schedule. So, it is
  recommended that you run it from an IDE and/or install it as
  an [editable package](https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs). If you do not wish to
  install as an editable package then please remove all `-e` flags.

#### Installing `Julia` and `diffeqpy` (time integration suite)

By default, TidalPy will utilize the `SciPy.integrate` package for solving differential equations. However, it may be
more optimal to use the `Julia` language which has ODE integrators that can be called directly from Python. In order to
use this functionality you will need to install `Julia` and the `diffeqpy` package.

* Install the Julia language from [https://julialang.org/downloads/](https://julialang.org/downloads/)
* Add Julia's directory and its `bin` subdirectory to your system's path.
* Open an elevated ("as administrator") terminal, command prompt, or powershell.
    * If you are using a virtual Python environment make sure it is active.

*As of TidalPy v0.3.0, the `diffeqpy` that is available from pypi is not up-to-date with the version found on the
project's github page. TidalPy uses functionality that is only available from this new version. If you would like to use
Julia with TidalPy you must clone the github version. Keep in mind this is an unreleased version so more bugs are
likely.*

* Create a new directory to clone the `diffeqpy` repository.
    * Run `git clone https://github.com/SciML/diffeqpy`
    * With your browser navigated to the directory that contains `setup.py`, run `pip install .` (the trailing period is
      important).
* Open Python on your elevated terminal (the following steps may take a while to compile).
    * Run `import diffeqpy; diffeqpy.install()`
    * Run `import julia; julia.install()`

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
