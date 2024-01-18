# Contribute
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

## How to run TidalPy tests
After you have installed TidalPy and have made some changes, you should run tests (also build new ones for whatever changes that were made!). 
* Open a terminal/shell and navigate to the TidalPy directory.
* Ensure you have pytest installed (`pip install pytest`).
* Simply type `pytest` and hit enter. Pytest will automatically look for all test cases in the `/Tests/` directory and run them.
    * Note that multiple warnings may show while you are running tests. These are likely normal warnings and are expected. TidalPy will try to raise an `Exception` (which `pytest` should catch automatically) when there is a serious problem.

## Marking Issues or ToDos Inside the Source Code
Ideally, all bugs or future features should be added to TidalPy's [Github page](https://github.com/jrenaud90/TidalPy) (in the issues tab). However, if you would like to mark a location inside the source code where there is an issue, or a spot where a model could be improved in the future, use the following comment *tag* syntax:
```
# FIXME: Here is a function that has a bug that should be fixed ASAP
<Broken Code>

# TODO: Here is a bit of code related to a model that could be improved, perhaps there is a new research that has a better method.
<Working Code, but better model out there>

# OPT: This bit of code is working okay, but it could be coded up to be computationally faster. Use this tag and feel free to offer optimization advice.
<Working Code that is unoptimized>
```
These tags are particularly helpful when making a pull request that has issues you need help on.