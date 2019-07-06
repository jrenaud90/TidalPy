# TidalPy
TidalPy is a software suite designed to assist researchers in the calculation of tidal dissipation for rocky and icy worlds. 

## How to Install
Note: Installation has only been tested on Windows 10 and Ubuntu operating systems.

Pre-Install Requirements:
* Python version 3.6+ must be installed on your machine
* Working internet connection (for the initial installation only)

There are two methods for installing TidalPy. The first method directly installs the package to your python library. This is easy to do but will require uninstalling TidalPy each time a new update is pushed to github. This won't be an issue once TidalPy makes its way to PyPI, however, in the mean time, it is recommended to use Method 2. 
#### Method 1 - Direct Install:
_Easy to do, but harder to keep updated_
* Download the latest release at https://github.com/jrenaud90/TidalPy.git
* Unzip the package into an easy to a directory that is easy to access via a terminal
* Open a terminal and navigate to the outer most TidalPy directory (setup.py should be located there)
    * On Windows you need make sure your terminal is opened 'as administrator'
* Install the python package with the command 'python setup.py install'
* Python will then automatically install any required packages that your system is missing.
* Can now use 'import TidalPy' in your python environments.

#### Method 2 - Git Clone:
_Little harder to setup and use, but easy to update_
* Ensure you have the latest version of `git` or github installed
* Create an easy to access directory named `tidalpy`
* Open a terminal and navigate to this directory
* Clone the tidalpy git using `git clone https://github.com/jrenaud90/TidalPy.git`
    * Whenever you want to update TidalPy simply navigate to this directory and use `git pull`
* TidalPy source code will now be in your directory but Python does not know this, so using `import TidalPy` will only work if you use it from a terminal opened in this directory
* Python will also not check to see if you have all the dependencies to install those:
    * Open setup.py in a text editor like notepad and look for the list of strings after `install_requires=[`. These are the packages you must have.
    * In a terminal use `pip install <package 1> <package 2> ...` for each package in that list (do not type the version number or the >= symbols)
    * **Do not install** `burnman` this way, instead see the instructions in the next section.
* Much of these steps can be done automatically if you use an IDE that has version control (we recommend PyCharm)

#### Installing Burnman
As of 7-5-2019, the Burnman package has had several critical updates that have not made their way on to its official pyPI release. Therefore, in order to use TidalPy, you must download the latest BurnMan version via git.
* Ensure you have the latest version of `git` or github installed
* Create an easy to access directory named `burnman`
* Open a terminal and navigate to this directory
* Clone the lastest master branch of burnman using `git clone https://github.com/geodynamics/burnman`
* In a terminal navigate to the burnman directory that contains `setup.py`
* Install Burnman using `python setup.py install`

#### Test the Install
Once all packages have been installed, including Burnman, and you have completed install method 1 or 2, then you can test the installation by running `import TidalPy` in a python terminal.

Note: For method 2, the testing terminal must be navigated to the directory containing TidalPy.

## Basic Usage

## More Usage
