<h1 align="center">TidalPy</h1>
<h3 align="center">v0.3 Alpha</h3>

<p align="center">
    <a href="https://github.com/jrenaud90/TidalPy/actions/workflows/main_tests.yml?query=branch%3Amaster" /><img src="https://github.com/jrenaud90/TidalPy/actions/workflows/main_tests.yml/badge.svg?branch=master" /></a> <a href="https://codecov.io/gh/jrenaud90/TidalPy"><img src="https://codecov.io/gh/jrenaud90/TidalPy/branch/master/graph/badge.svg?token=35OY4ZLOA5"/></a><br />
    <a href="https://mybinder.org/v2/gh/jrenaud90/TidalPy/master?filepath=%2FCookbooks%2F"><img src="https://mybinder.org/badge_logo.svg" /></a>
    <a href="https://gitter.im/TidalPy/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge"><img src="https://badges.gitter.im/TidalPy/community.svg" /></a>
</p>

## Purpose
TidalPy is an open-source software suite designed to assist researchers in the semi-analytic calculation of tidal dissipation and subsequent orbit-spin evolution for rocky and icy worlds.

**TidalPy is intended to be a...**
* Black Box (referred to as "*OOP* scheme" for `Object-Oriented Programming` throughout the documentation)
    * TidalPy serves as simple to install (cross-platform) and, hopefully, simple to use package that users can pick up and hit the ground running.
    * The OOP scheme performs many calculations with very little input from the user. The major drawbacks are performance (in some situations) and that assumptions have been made that are opaque to the user without some digging. 
* Tool Box (referred to as "*Functional* scheme")
    * TidalPy also contains many efficient functions to perform calculations relevant to tides and thermal-orbital coupling. These can be quickly imported and used in a custom script by a more experienced user.
        * In general, the functional scheme will have much higher performance, flexibility, and extensibility than OOP. It also generally makes assumptions more visible to the user. 
          
*Once you are comfortable with TidalPy, it is usually a good idea to mix the two schemes: take some aspects of OOP that you don't want to deal with and build on them with some of TidalPy's or your own functions and code.*

### Limitations

The major limitations of the current version of TidalPy are...
* A multilayer model has now been implemented, but it is not currently part of the OOP scheme.
* Chemical and phase changes within a planet's layers have not yet been implemented.

### Related Software

Below is a list (non-exhaustive) of publicly available software that performance similar or parallel calculations as TidalPy.

* Are you interested in the habitability of a planet? With considerations of tides, atmospheres, water content, solar interactions? Check out...
    * [VPLanet](https://github.com/VirtualPlanetaryLaboratory/vplanet)
* Are you interested in the orbital evolution of multiple planets with each planet influencing one another? Consider an N-body approach like...
    * [Posidonius (formerly Mercury-T)](https://github.com/marblestation/posidonius)
    * [ReboundX](https://github.com/dtamayo/reboundx)
* Don't care about tides or orbital dynamics but are instead interested in interior structure and composition of planets?
    * [BurnMan](https://github.com/geodynamics/burnman)
    * [PerpleX](http://www.perplex.ethz.ch/)
* Are you interested in tides, interiors, *and* the chemical evolution of small worlds but don't care about non-synchronous rotation or compressibility of planets?
    * [IcyDwarf](https://github.com/MarcNeveu/IcyDwarf)

However, if you want high fidelity tidal, orbital, spin, and interior models --- then you have come to the right place! Read below for instructions on how to install and use TidalPy.

## How to Install

### Compatibility
*As of TidalPy v0.3.0a*:
* **Win10**: *Installation & tests passed.*
* **MacOS (Catalina)**: *TBD*
* **CentOS7**: *TBD*

### Pre-Install

Pre-Install Requirements:
* Python version 3.8+ must be installed on your machine.
    * It is highly recommended that you use the [Anaconda](https://www.anaconda.com/distribution/) distribution of Python. This has pre-compiled binaries for several packages that TidalPy uses and will generally negate a lot of potential headaches. If you don't want to use Anaconda you can find the regular Python distribution [here](https://www.python.org/).
    * Make sure that your Python (Anaconda or regular) is 64-bit if you are on a 64-bit machine.
* Working internet connection (for the initial installation only).
* The `setuptools` package is required before TidalPy can be installed. Usually it is automatically installed, but if you are starting with a clean virtual environment it may not have been.
    * For Anaconda: `conda install setuptools`
    * Or for regular Python: `pip install setuptools`
* Unless you plan to download the source code from github directly, make sure you have [git](https://git-scm.com/downloads) or [github](https://desktop.github.com/) installed on your machine.

### Install
The current version of TidalPy is in Alpha and will receive many updates on a relatively fast schedule. It is, therefore, recommended that you run it from an IDE (more on that below) and/or install it as an [editable package](https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs). If you do not wish to install as an editable package then please remove all `-e` flags used below.

* Get the latest version of TidalPy from Github
    * Open a terminal and navigate to an easy-to-access directory where you would like to install TidalPy.
    * Clone the TidalPy git using `git clone https://github.com/jrenaud90/TidalPy.git`.
        * Whenever you want to update TidalPy simply navigate to this directory and use `git pull` (to pull from the master branch; other branches are not recommended).
        * Since TidalPy is in early development, it is recommended you check for updates regularly. Updates will **not** download automatically. 
        * Always make a backup of the TidalPy installation directory in case new versions break whatever you were working on.

**Before continuing:** You might consider using a new virtual environment so that these new packages do not overwrite packages that you may be using for different projects on your machine.
* Install Burnman:
    * The latest released version of [Burnman](https://github.com/geodynamics/burnman) has several issues that will prevent TidalPy from running. However, the latest update on its github is working.
        * Install from github: `python -m pip install git+https://github.com/geodynamics/burnman.git`
    
* Install TidalPy:
    * Using a terminal, navigate to the TidalPy directory that contains `setup.py` and then:
        * For Anaconda Python:
            * Run `conda install --file conda_requirements.txt; pip install -e .` *(That trailing period is important, don't leave it out!)*
        * For non-Anaconda Python:
            * Run `pip install -e .` *(That trailing period is important, don't leave it out!)*
    * This will automatically ensure that your python installation (Anaconda or regular) has the required third party packages.
* Test your installation:
    * Navigate to the TidalPy directory that contains `setup.py` in a terminal.
    * Ensure you have `pytest` package installed (`conda install pytest` or `pip install pytest`).
    * Run pytest by simply using the command `pytest` from your terminal:
        * Running all the tests can take a while (currently 3-10 minutes), if all you are interested in is checking that TidalPy installed correctly then you can let pytest check first handful or so if they are passing then you can quit the test suite early.
        * If no errors show up (warnings are okay and expected) then the first check is good.
    * Open a new terminal *not in the TidalPy directory* (e.g., your desktop).
        * Run `python` and then try to `import TidalPy`; if that works try the command `TidalPy.version` if you do not get any import errors, and the version number is as you expect, then TidalPy was successfully installed.
    
#### Using TidalPy from an IDE
A good Integrated Development Environment can automatically set paths to TidalPy and allows you to use TidalPy without actually "installing" it. If you are comfortable with IDEs then this may be an easier way to use TidalPy, especially during its alpha phase.

#### Installing `Julia` and `diffeqpy` (integration suite)
By default TidalPy will utilize the `SciPy.integrate` package for solving differential equations. However, it may be more optimal to use the `Julia` language which has ODE integrators that can be called directly from Python. In order to use this functionality you will need to install `Julia` and the `diffeqpy` package.

* Install the Julia language from [https://julialang.org/downloads/](https://julialang.org/downloads/)
* Add Julia's directory and its `bin` subdirectory to your system's path.
* Open an elevated ("as administrator") terminal, command prompt, or powershell.
    * If you are using a virtual Python environment make sure it is active.
    
*As of TidalPy v0.3.0, the `diffeqpy` that is available from pypi is not up to date with the version found on the project's github page. TidalPy uses functionality that is only available from this new version. If you would like to use Julia with TidalPy you must close the github version. Keep in mind this is an unreleased version so more bugs are likely.*
      
* Create a new directory to clone the `diffeqpy` repository.
    * Run `git clone https://github.com/SciML/diffeqpy`
    * With your browser navigated to the directory that contains `setup.py`, run `pip install .` (the trailing period is important).
* Open Python on your elevated terminal (the following steps may take a while to compile).
    * Run `import diffeqpy; diffeqpy.install()`
    * Run `import julia; julia.install()`

## How to Use
Check out the `Documentation\Getting Started.md` file. This is pretty bare bones at the moment but offers some basic info about TidalPy.
For now the best way to learn how to use TidalPy is by checking out the `Demos` directory. There are "beginner" [Jupyter notebooks](https://jupyter.org/) that are a great starting point.

## Using TidalPy for Science
TidalPy has been used in several studies already, and we encourage you to use it in yours. We would appreciate you include a link back to this [page](https://github.com/jrenaud90/TidalPy) and cite one of the papers below (if you utilized a specific package). We also would love to see where TidalPy is being used! Please feel free to send us an email: [TidalPy@gmail.com](mailto:TidalPy@gmail.com) when a paper or presentation utilized TidalPy.
Anyone is welcome to make forks or copies of TidalPy as long as their work references back to this page. License information can be found at the end of this file.

### TidalPy's Science
The science used in TidalPy is described in the following papers (and references therein):
* Rheological Modeling Package:
    * [Tidally Heated Terrestrial Exoplanets: Viscoelastic Response Models](https://ui.adsabs.harvard.edu/abs/2009ApJ...707.1000H/abstract)
    * [Increased Tidal Dissipation Using Advanced Rheological Models](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...98R/abstract)
* Non-synchronous Rotation Evolution and High Eccentricity Truncation Packages:
    * [Tidal Dissipation in Dual-Body, Highly Eccentric, and Non-synchronously Rotating System](https://ui.adsabs.harvard.edu/abs/2020arXiv201011801R/abstract)
    * [Tidal Evolution of the Keplerian Elements](https://ui.adsabs.harvard.edu/abs/2019CeMDA.131...30B/abstract)
* Third Party Software:
    * *Interior Model*: [BurnMan](https://github.com/geodynamics/burnman)
    * *Integration Routines*: [diffeqpy](https://github.com/SciML/diffeqpy), [Julia DiffEq](https://diffeq.sciml.ai/v2.0/)
    * *CVD Conscious Color Maps*: [Geodynamic Color Maps](http://doi.org/10.5281/zenodo.1243862)

## Contribute
TidalPy is in early alpha and there are lots of areas where it can improve! If you are interested in helping out, please check out the information in `Documentation\Contribute.md`.

**Found a bug or have an idea for a new feature?**
* Go to TidalPy's [Github page](https://github.com/jrenaud90/TidalPy) and click the "Issues" tab then make a new report.
    * If you ran into a bug please include a code snippet (in markdown: code is designated by Grave accents surrounding the text) that reproduces the error (please keep this snippet as concise as possible).
    * It is helpful to triage issues when they are made. If you think you know the severity of a bug or can provide any other *at-a-glance* context, consider adding a "label" (right-hand side of the github issue form) to the issue.

## License Information
You are welcome to make a copy/fork of TidalPy and make modifications assuming the following conditions are met:
* Links are included that point back to this [page](https://github.com/jrenaud90/TidalPy).
* Any software derived from TidalPy must remain open-source and non-commercial.

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this license, visit [http://creativecommons.org/licenses/by-nc-sa/4.0/](http://creativecommons.org/licenses/by-nc-sa/4.0/) or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
