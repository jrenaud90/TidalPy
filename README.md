# TidalPy - v0.1 Beta
## Purpose
TidalPy is an open-source software suite designed to assist researchers in the semi-analytic calculation of tidal dissipation and evolution for rocky and icy worlds.

**TidalPy is intended to be good as a...**
* Black Box (referred to as "*OOP* scheme" for `Object-Oriented Programming` throughout the documentation)
    * `TidalPy` serves as simple to install (cross-platform) and, hopefully, simple to use package that non-experts (to tides and planetary interiors) can pick up and use quickly.
    * The OOP scheme performs many calculations with very little input from the user.
* Tool Box (referred to as "*Functional* scheme")
    * `TidalPy` also contains many efficient functions to perform calculations relevant to tides and thermal-orbital coupling. These can be quickly imported and used in a custom script by a more experienced user.
        * In general, the functional scheme will have much higher performance than OOP.
        * Once you are comfortable with `TidalPy`, it is usually a good idea to mix the two schemes: take some aspects of OOP that you don't want to deal with and build off of them with some of `TidalPy`'s or your own functions and code.

**Related Software**

Below publicly available software that is similar to `TidalPy`.

* Are you more interested in the orbital evolution than the interior evolution of a planet or a system of multiple planets?
    * [VPLanet](https://github.com/VirtualPlanetaryLaboratory/vplanet)
    * [Posidonius (formerly Mercury-T)](https://github.com/marblestation/posidonius)
* Do you not care about tides but are instead interested in interior structure and composition of planets?
    * [Burnman](https://github.com/geodynamics/burnman)
* Are you interested in tides, interiors, *and* the chemical evolution of icy worlds but don't care too much about non-synchronous rotation?
    * [IcyDwarf](https://github.com/MarcNeveu/IcyDwarf)
* Do you want all of the above (minus: chemical evolution and N-body integration)? Then you have come to the right place! Look below on how to download `TidalPy`.

## How to Install
### Pre-Install
Note: Installation has only been tested on Windows 10 and Ubuntu operating systems.

Pre-Install Requirements:
* Python version 3.7+ must be installed on your machine
    * I usually recommend the [Anaconda](https://www.anaconda.com/distribution/) distribution of Python, but lately I have been having lots of issues with it and Windows 10. If you don't want to use Anaconda you can find the regular Python distribution [here](https://www.python.org/).
    * Make sure that your Python (Anaconda or regular) is 64-bit if you are on a 64-bit machine.
* Working internet connection (for the initial installation only)

### Install
The current version of TidalPy is going to receive many updates on a relatively fast schedule. It is therefore recommended that you run it from an IDE (more on that below) and/or install it as an editable package.

#### Install as an editable package
* Git clone the latest version from Github
    * Ensure you have [git](https://git-scm.com/downloads) or [github](https://desktop.github.com/) installed on your machine
    * Open a terminal and navigate to an easy-to-access directory that you want to store TidalPy at
    * Clone the tidalpy git using `git clone https://github.com/jrenaud90/TidalPy.git`
        * Whenever you want to update TidalPy simply navigate to this directory and use `git pull`
        * Since TidalPy is in early development, it is recommended you check for updates regularly. Updates will **not** download automatically.
* TidalPy source code should now be in your directory but Python does not know this, so using `import TidalPy` will only work if performed from a terminal that has been navigated to this directory.
* To install TidalPy so it can be accessed from any terminal location:
    * Navigate to the TidalPy directory that contains `setup.py` in a terminal
    * Run `pip install -r requirements.txt`
        * This will ensure that your python installation has the required packages to run TidalPy.
        * **Before you run this:** You might consider starting a new virtual environment so that these new packages do not overwrite packages that you may be using for different projects on your machine
    * After that completes successfully you can now install TidalPy as an editable package
    * Ensure you are still navigated to the same directory and then run `pip install -e .`
        * That trailing period is important, don't leave it out!
* Test your install
    * Navigate to the TidalPy directory that contains `setup.py` in a terminal
    * Ensure you have `pytest` package installed (`pip install pytest`)
    * Run pytest by simply using the command `pytest` from your terminal
        * If no errors show up (warnings are okay and expected) then the first check is good.
    * Open a new terminal *not in the TidalPy directory* (you can use the desktop for instance).
        * Run `python` and then try to `import TidalPy` if you do not get any import errors then TidalPy was successfully installed.

#### Using TidalPy from an IDE
A good IDE can automatically set paths to TidalPy and allows you to use TidalPy without actually installing it. If you are comfortable with IDEs then this may be an easier way to get and use TidalPy (we recommend the PyCharm IDE).

## How to Use
Coming Soon! For now, check out the `Cookbooks` and `Documentation` directories.

Eventually there will be getting started information in the `Documentation/Getting Started.md`.

## Contribute
TidalPy is in very early alpha and needs lots of work and help! Please look at the following to see how you can help.

**Found a bug or have an idea for a new feature?**
* Go to TidalPy's [Github page](https://github.com/jrenaud90/TidalPy) and click the "Issues" tab then make a new report.
    * If you ran into a bug please include a code snippet (in markdown: code is designated by Grave accents surrounding the text) that reproduces the error (please keep this snippet as concise as possible).
    * It is helpful to triage issues when they are made. If you think you know the severity of a bug or can provide any other *at-a-glance* context, consider adding a "label" (right-hand side of the github issue form) to the issue.

**Want to contribute directly?**

TidalPy is an open-source project and depends upon new contributions from the community. If you would like to contribute, please follow these steps:
* Find an issue or new feature you want to tackle (a good place to start is to look in the **Issues tab** for the label "Beginner")
* Fork the latest version of the Master branch into a new branch on your Github account.
* Work on making the code corrections that fix the issue or implement the new feature.
   * If you are implementing a new feature please try to build test cases in the /Tests/ directory so that future bugs can be quickly spotted.
* [Rebase](https://www.atlassian.com/git/tutorials/merging-vs-rebasing) from the master branch (to ensure your code works with any recent changes to TidalPy since you forked).
* Run TidalPy's tests (see below) to ensure your changes did not negatively impact other parts of the code.
* Assuming all tests pass, make a new [Pull Request](https://help.github.com/en/articles/creating-a-pull-request-from-a-fork) at github.com/jrenaud90/TidalPy

Do not hesitate to make a pull request *before* you have a working solution to the issue. A pull request provides a forum where you can get feedback or help on your changes from other contributors before things are finalized.

### How to run TidalPy tests
After you have installed TidalPy and have made some changes, you should run tests (also build new ones for whatever changes that were made!). 
* Open a terminal/shell and navigate to the TidalPy directory
* Ensure you have pytest installed (`pip install pytest`)
* Simply type `pytest` and hit enter. Pytest will automatically look for all test cases in the `/Tests/` directory and run them.
    * Note that multiple warnings for `invalid escape sequence` for the `burnman` package may show up. You can ignore these or tell pytest to ignore them using the flag: `pytest -W ignore::DeprecationWarning`

## Current Build Info

![Travis-ci](https://travis-ci.com/jrenaud90/TidalPy.svg?token=hTmV5nwCsy8qF9GmqKXP&branch=master)
