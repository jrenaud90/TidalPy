# TidalPy
**v0.2 Alpha**  ![Travis-ci](https://travis-ci.com/jrenaud90/TidalPy.svg?token=hTmV5nwCsy8qF9GmqKXP&branch=master) 

## Purpose
TidalPy is an open-source software suite designed to assist researchers in the semi-analytic calculation of tidal dissipation and subsequent orbit-spin evolution for rocky and icy worlds.

** TidalPy is intended to be a...**
* Black Box (referred to as "*OOP* scheme" for `Object-Oriented Programming` throughout the documentation)
    * TidalPy serves as simple to install (cross-platform) and, hopefully, simple to use package that non-experts (to tides and planetary interiors) can pick up and use quickly.
    * The OOP scheme performs many calculations with very little input from the user. The major drawbacks are peformance (in some situations) and that assumptions have been made that are opaque to the user without some digging. 
* Tool Box (referred to as "*Functional* scheme")
    * TidalPy also contains many efficient functions to perform calculations relevant to tides and thermal-orbital coupling. These can be quickly imported and used in a custom script by a more experienced user.
        * In general, the functional scheme will have much higher performance, flexibility, and extensibility than OOP. It also ensures that assumptions in the various models are made clear to the use. 
          
*Once you are comfortable with TidalPy, it is usually a good idea to mix the two schemes: take some aspects of OOP that you don't want to deal with and build on them with some of TidalPy's or your own functions and code.*

### Limitations

The major limitations of the current version of TidalPy are...
* The interiors of planets are treated as 0-dimensional layers. Multiple layers are allowed, but a multi-shell dissipation approach has not been implemented.
* Effects caused by compressibility are not considered.
* Chemical and phase changes within a planet's layers have not yet been implemented.

### Related Software

Below is a list (non-exhaustive) of publicly available software that performance similar or parallel calculations as TidalPy.

* Are you interested in the orbital evolution of multiple planets with each planet influencing one another? Consider an N-body approach like...
    * [VPLanet](https://github.com/VirtualPlanetaryLaboratory/vplanet)
    * [Posidonius (formerly Mercury-T)](https://github.com/marblestation/posidonius)
    * [ReboundX](https://github.com/dtamayo/reboundx)
* Don't care about tides or orbital dynamics but are instead interested in interior structure and composition of planets?
    * [BurnMan](https://github.com/geodynamics/burnman)
    * [PerpleX](http://www.perplex.ethz.ch/)
* Are you interested in tides, interiors, *and* the chemical evolution of small worlds but don't care about non-synchronous rotation?
    * [IcyDwarf](https://github.com/MarcNeveu/IcyDwarf)

Instead, if you are okay with a less accurate interior model and no multi-planet support, then you have come to the right place! Look below on how to install and use TidalPy.

## How to Install

### Compatibility
*As of TidalPy v0.2.1.a3*:
* **Win10**: *Installation & tests passed.*
* **MacOS (Catalina)**: *Installation & tests passed.*
* **CentOS7**: *Conflicts during some parts of conda install, however after a while installation completes and all tests pass.*

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
The current version of TidalPy is in Alpha and will receive many updates on a relatively fast schedule. It is, therefore, recommended that you run it from an IDE (more on that below) and/or install it as an [editable package](https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs).

#### Install as an editable package
* Get the latest version of TidalPy from Github
    * Open a terminal and navigate to an easy-to-access directory where you would like to install TidalPy.
    * Clone the TidalPy git using `git clone https://github.com/jrenaud90/TidalPy.git`.
        * Whenever you want to update TidalPy simply navigate to this directory and use `git pull` (to pull from the master branch; other branches are not recommended).
        * Since TidalPy is in early development, it is recommended you check for updates regularly. Updates will **not** download automatically. 
        * Always make a backup of the TidalPy installation directory in case new versions break whatever you were working on.
* TidalPy source code should now be in your directory but Python does not know this, so using `import TidalPy` inside of Python will only work if performed from a terminal that is navigated to this directory.
* To install TidalPy so that it can be accessed from the terminal anywhere:
    * **Before continuing:** You might consider using a new virtual environment so that these new packages do not overwrite packages that you may be using for different projects on your machine.
    * Using a terminal, navigate to the TidalPy directory that contains `setup.py` and then run `pip install -e .` *(That trailing period is important, don't leave it out!)*
    * This will automatically ensure that your python installation (Anaconda or regular) has the required third party packages that TidalPy requires.
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

## How to Use
TBA!

For now, check out the `Cookbooks` and `Documentation` directories. There are "beginner" [Jupyter notebooks](https://jupyter.org/) that are a great starting point.

Eventually there will be more "getting started information" in the `Documentation/Getting Started.md`.

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
