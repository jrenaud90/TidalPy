# TidalPy - v0.2 Alpha
![Travis-ci](https://travis-ci.com/jrenaud90/TidalPy.svg?token=hTmV5nwCsy8qF9GmqKXP&branch=master) 
## Purpose
TidalPy is an open-source software suite designed to assist researchers in the semi-analytic calculation of tidal dissipation and evolution for rocky and icy worlds.

**TidalPy is intended to be good as a...**
* Black Box (referred to as "*OOP* scheme" for `Object-Oriented Programming` throughout the documentation)
    * `TidalPy` serves as simple to install (cross-platform) and, hopefully, simple to use package that non-experts (to tides and planetary interiors) can pick up and use quickly.
    * The OOP scheme performs many calculations with very little input from the user.
* Tool Box (referred to as "*Functional* scheme")
    * `TidalPy` also contains many efficient functions to perform calculations relevant to tides and thermal-orbital coupling. These can be quickly imported and used in a custom script by a more experienced user.
        * In general, the functional scheme will have much higher performance, flexibility, and extensibility than OOP.
        * Once you are comfortable with `TidalPy`, it is usually a good idea to mix the two schemes: take some aspects of OOP that you don't want to deal with and build off of them with some of `TidalPy`'s or your own functions and code.

**Related Software**

Below is a list of (non-exhaustive) publicly available software that is similar to `TidalPy`.

* Are you more interested in the orbital evolution of multiple planets with each planet influencing one another? Consider a N-body approach like...
    * [VPLanet](https://github.com/VirtualPlanetaryLaboratory/vplanet)
    * [Posidonius (formerly Mercury-T)](https://github.com/marblestation/posidonius)
    * [ReboundX](https://github.com/dtamayo/reboundx)
* Don't care about tides or orbital dynamics but are instead interested in interior structure and composition of planets?
    * [BurnMan](https://github.com/geodynamics/burnman)
* Are you interested in tides, interiors, *and* the chemical evolution of icy worlds but don't care too much about non-synchronous rotation?
    * [IcyDwarf](https://github.com/MarcNeveu/IcyDwarf)
* Instead, if you are okay with a less accurate interior model and no multi-body support, then you have come to the right place! Look below on how to download and use `TidalPy`.

## How to Install
### Pre-Install
Note: Installation has only been tested on Windows 10 and Ubuntu operating systems.

Pre-Install Requirements:
* Python version 3.7+ must be installed on your machine.
    * It is highly recommended that you use the [Anaconda](https://www.anaconda.com/distribution/) distribution of Python. This has pre-compiled binaries for several packages that TidalPy uses and will generally negate a lot of potential installation headaches. If you don't want to use Anaconda you can find the regular Python distribution [here](https://www.python.org/).
    * Make sure that your Python (Anaconda or regular) is 64-bit if you are on a 64-bit machine.
* Working internet connection (for the initial installation only).

### Install
The current version of TidalPy is in Alpha and will receive many updates on a relatively fast schedule. It is, therefore, recommended that you run it from an IDE (more on that below) and/or install it as an [editable package](https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs).

#### Install as an editable package
* Get the the latest version from Github
    * Ensure you have [git](https://git-scm.com/downloads) or [github](https://desktop.github.com/) installed on your machine.
    * Open a terminal and navigate to an easy-to-access directory where you would like to install TidalPy.
    * Clone the TidalPy git using `git clone https://github.com/jrenaud90/TidalPy.git`.
        * Whenever you want to update TidalPy simply navigate to this directory and use `git pull` (to pull from the master branch; other branches are not recommended).
        * Since TidalPy is in early development, it is recommended you check for updates regularly. Updates will **not** download automatically. 
        * Always make a backup of the TidalPy installation directory in case new versions break whatever you were working on.
* TidalPy source code should now be in your directory but Python does not know this, so using `import TidalPy` will only work if performed from a terminal that is navigated to this directory.
* To install TidalPy so it can be accessed from the terminal anywhere:
    * Using a terminal, navigate to the TidalPy directory that contains `setup.py` and then run `pip install -e .`.
        * That trailing period is important, don't leave it out!
        * Optionally, you may add the `-v` flag (before `-e`) to see more installation information --- but this tends to be too much info.
        * This will automatically ensure that your python installation (Anaconda or regular) has the required third party packages that TidalPy requires.
        * **Before you run this:** You might consider using a new virtual environment so that these new packages do not overwrite packages that you may be using for different projects on your machine.
* Test your install:
    * Navigate to the TidalPy directory that contains `setup.py` in a terminal.
    * Ensure you have `pytest` package installed (`conda install pytest` or `pip install pytest`).
    * Run pytest by simply using the command `pytest` from your terminal:
        * Running all the tests can take some time, if all you are interested in is checking that TidalPy installed correctly then you can let pytest check the first dozen or so if they are passing then you can quit the test suite early.
        * If no errors show up (warnings are okay and expected) then the first check is good.
    * Open a new terminal *not in the TidalPy directory* (you can use the desktop for instance).
        * Run `python` and then try to `import TidalPy`; if that works try the command `TidalPy.__version__` if you do not get any import errors, and the version number is as you expect, then TidalPy was successfully installed.

#### Using TidalPy from an IDE
A good Integrated Development Environment can automatically set paths to TidalPy and allows you to use TidalPy without actually "installing" it. If you are comfortable with IDEs then this may be an easier way to use TidalPy, especially during its alpha phase.

## How to Use
Coming Soon! For now, check out the `Cookbooks and Scripts` and `Documentation` directories. There are "beginner" [Jupyter notebooks](https://jupyter.org/) that are a great starting point.

Eventually there will be more "getting started information" in the `Documentation/Getting Started.md`.

## Using TidalPy for Science
TidalPy has been used in several studies already and we encourage you to use it in yours. We would appreciate you include a link back to this [page](https://github.com/jrenaud90/TidalPy) and cite one of the papers below (if you utilized a specific package). We also would love to see where TidalPy is being used! Please feel free to send us an email: [TidalPy@gmail.com](mailto:TidalPy@gmail.com) when a paper or presentation utilized TidalPy.
Anyone is welcome to make forks or copies of TidalPy as long as their work references back to this page. License information can be found at the end of this file.

### TidalPy's Science
The science used in TidalPy is described in the following papers (and references therein):
* Rheological Modeling Package:
    * [Tidally Heated Terrestrial Exoplanets: Viscoelastic Response Models](https://ui.adsabs.harvard.edu/abs/2009ApJ...707.1000H/abstract)
    * [Increased Tidal Dissipation Using Advanced Rheological Models](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...98R/abstract)
* Non-synchronous Rotation Evolution and High Eccentricity Truncation Packages:
    * Tidal Dissipation in Dual-Body, Highly Eccentric, and Non-synchronously Rotating System (*in review*).
* Third Party Software:
    * [BurnMan](https://github.com/geodynamics/burnman)
    * [diffeqpy](https://github.com/SciML/diffeqpy)

## Contribute
TidalPy is in early alpha and needs lots of improvements and help! If you are interested in helping out, please check out the following.

**Found a bug or have an idea for a new feature?**
* Go to TidalPy's [Github page](https://github.com/jrenaud90/TidalPy) and click the "Issues" tab then make a new report.
    * If you ran into a bug please include a code snippet (in markdown: code is designated by Grave accents surrounding the text) that reproduces the error (please keep this snippet as concise as possible).
    * It is helpful to triage issues when they are made. If you think you know the severity of a bug or can provide any other *at-a-glance* context, consider adding a "label" (right-hand side of the github issue form) to the issue.

**Want to contribute directly?**

TidalPy is an open-source project and improves with each new contribution from the community. If you would like to contribute, please follow these steps:
* Find an issue or new feature you want to tackle (a good place to start is to look in the **Issues tab** for the label "Beginner").
* Fork a new branch off of the latest version of the `Master` branch on your Github account or local machine.
* Work on making the code corrections that fix the issue or implement a new feature.
    * If you are implementing a new feature, please try to build test cases in the /Tests/ directory so that future bugs can be quickly squashed.
* [Rebase](https://www.atlassian.com/git/tutorials/merging-vs-rebasing) from the master branch (to ensure your code works with any recent changes to TidalPy since you preformed your initial fork).
* Run TidalPy's tests (see below) to ensure your changes did not negatively impact other parts of the code.
* Assuming all tests pass, make a new [Pull Request](https://help.github.com/en/articles/creating-a-pull-request-from-a-fork) at github.com/jrenaud90/TidalPy

Do not hesitate to make a pull request *before* you have a working solution to the issue. A pull request provides a forum where you can get feedback or help on your changes from other contributors before things are finalized.

### How to run TidalPy tests
After you have installed TidalPy and have made some changes, you should run tests (also build new ones for whatever changes that were made!). 
* Open a terminal/shell and navigate to the TidalPy directory.
* Ensure you have pytest installed (`pip install pytest`).
* Simply type `pytest` and hit enter. Pytest will automatically look for all test cases in the `/Tests/` directory and run them.
    * Note that multiple warnings may show while you are running tests. These are likely normal warnings and are expected. TidalPy will try to raise an `Exception` (which `pytest` should catch automatically) when there is a serious problem.

## License Information
You are welcome to make a copy/fork of TidalPy and make modifications assuming the following conditions are met:
* Links are included that point back to this page.
* Any software derived from TidalPy must remain open-source and non-commercial.

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this license, visit [http://creativecommons.org/licenses/by-nc-sa/4.0/](http://creativecommons.org/licenses/by-nc-sa/4.0/) or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
