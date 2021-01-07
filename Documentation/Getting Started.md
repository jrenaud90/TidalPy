# Getting Started with TidalPy

## Installation
Please see the package's `readme.md` file for the most up to date installation information.
You can find this file, along with the project's source code, at its [GitHub page](https://github.com/jrenaud90/TidalPy).

## After Installation
At the current stage of `TidalPy` development things are changing too rapidly to make a comprehensive getting started guide, sorry!

For now look to the `Cookbooks` examples on how to use `TidalPy`. If you find any issues or have an idea about a new feature then feel free to leave a new github issue [here](https://github.com/jrenaud90/TidalPy/issues). Stuck on a problem or have a question? Feel free to contact joe.p.renaud@gmail.com try to include a minimal working snippet of code (if relevant to your question).

## TidalPy OOP structure
In addition to basic function calls, TidalPy contains an object-oriented system to quickly build planets, calculate physics, and return results. Under the hood this system is fairly complicated, making changes to one object may end up breaking functionality with other objects. This is particularly problemsome in the current version of TidalPy as the `Tests` suite does not fully capture all OOP operations. We have created a flow chart to help explain the OOP structure and hopefully mitigate issues. This chart can be found at `Documentation\OOP Structure (<version #>).pdf`.

This flow chart should be updated (via changing the `.svg`) whenever there is a change to TidalPy's OOP structure.

## Functionality still under development
In this early iteration of `TidalPy` you will occasionally come across functions or operations that are not fully implemented or are in a very early alpha state. These will generally be denoted with a `_dev` suffix. Be very careful when using any functionality that is not fully implemented and tested.