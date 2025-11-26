---
title: 'TidalPy: Software Suite for Solving Problems in Tidal Dynamics'
tags:
  - Python
  - C++
  - Cython
  - Astronomy
  - Material physics
  - Orbital Mechanics
  - Dynamics
  - Tides
  - Solid body tides
  - Tidal disruption
  - Tidal heating
  - Rheology
  - Love numbers
authors:
  - name: Joe P. Renaud
    orcid: 0000-0002-8619-8542
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
    corresponding: true
affiliations:
  - name: University of Maryland, College Park, Maryland USA
    index: 1
    ror: 047s2c258
  - name: NASA Goddard Space Flight Center, Greenbelt, Maryland USA
    index: 2
    ror: 0171mag52
date: 22 November 2025
bibliography: paper.bib

---



# Summary

`TidalPy` is an open-source Python package for modeling tidal deformation, internal dissipation,
and rotational–orbital evolution in planetary systems. The software combines a user-friendly Python
interface with performance-critical routines implemented in C++ and Cython. Key capabilities include
dynamic Love number computation, advanced rheological models, and coupled spin–orbit evolution.
TidalPy addresses a gap in existing tools by providing these capabilities within a single accessible
framework that easily interfaces with other popular Python packages used in the field. It has been 
extensively tested in a variety of studies examining tidal dissipation in Solar System planets and moons,
as well as exoplanets.

# Statement of need

`TidalPy` is a open-source Python package developed at NASA Goddard Space Flight
Center and University of Maryland College Park to support tidal research in Planetary Science.
The package provides a flexible, accessible, and performant toolkit for solving problems in
tides and tidal dynamics. The same tides that cause Earth's ocean to rise twice each day
can churn the interiors of other planets and moons to the point that significant fractions
of their bulk can melt, greatly altering the long-term thermal evolution of these worlds.
The energy that drives this heat originates in the orbits and rotations of the planet and
its hosts. TidalPy provides functions and frameworks to apply the latest tidal modeling
theories and methods to a wide variety of Solar System and exoplanetary worlds.

# Overview

TidalPy is written primarily in Python, with performance-critical components implemented
in C++ and Cython [REF]. Its API is designed to be intuitive and consistent with modern
conventions, enabling both early career and experienced researchers to quickly learn its
syntax and incorporate it in their scientific projects. TidalPy joins a robust community
of other packages that perform similar calculations [REF] and expands on this prior
work in three major areas described in this section.

## Love Number Solver (RadialSolver Module)

_Learn more about TidalPy's RadialSolver Module [here](REF)_

**Love Numbers** quantify a planet or moon's ability to respond to tidal forces [@love1911; @shida1912]. 
They are dynamic and depend on many physical factors such as a world's thermal state, physical
structures (e.g., a presence of a solid or liquid core), past stress events, and orbit/spin state.
These numbers can be measured, albeit with difficulty (particularly challenging
if we are unable to send a fly-by or orbiting satellite). Therefore, it is useful to perform forward modeling
utilizing our best estimates of a world's structure and composition to provide a range of tidal efficiency.

TidalPy provides a Love number solver that uses information about a planet's interior structure
and thermal state to estimate these numbers. A user can turn on or off a variety of assumptions
to determine their impact. This solver can be used in other routines to, for example, determine
the effect of long-term heating on a world's tidal dissipation or in a Markov Chain Monte-Carlo scheme to
predict its most likely interior structure. 

![The percent difference between Love number using the dynamic and static assumption is shown for a icy moon with a significant ocean layer as a function of tidal forcing period. Several reference periods are shown to give a sense of when dynamic tides may be important to consider. \label{fig:dynamic}](Papers\2025-JOSS\figures\dynamic_vs_static_tides.png){ width=50% }

TidalPy's solver uses a shooting method [@TakeuchiSaito1972] to find tidal and loading
Love numbers. This approach is advantageous as it enables more advanced physics, providing 
a more accurate description of a world. Specifically, TidalPy's solver allows for: liquid layers
and oceans, bulk compressibility (See \autoref{fig:compressibility}), and dynamic tides (See \autoref{fig:dynamic}).
This additional physics has been shown to be important for certain worlds during certain epochs. TidalPy's
Love number solver has been benchmarked against others tools that provide some of the same functionality including
`ALMA3` [@alma3; @pyalma] and `LoadDef` [@loaddef]. Other tools exist that, unlike the current version of TidalPy, 
can calculate multi-dimensional Love numbers [@Qin+2014nov; @Rovira-Navarro+2024maya; Berne+2023nov].

![Bulk dissipation can lead to significant differences in both the Tidal (left) and Loading (right) love numbers in this simplified Venus model. \label{fig:compressibility}](Papers\2025-JOSS\figures\compressibility_effects_venus.png){ width=50% }

## Advanced Rheological Modeling (Rheology Module)

_Learn more about TidalPy's Rheology Module [here](REF)_

The calculation of tidal Love numbers requires knowing the viscoelastic state of a planet. 
This can be described through the shear and bulk modulus as well as the shear and bulk viscosity. 
The former describe how sound waves travel through a planet's bulk, the later describe how material flows
on long timescales. Linking these properties to tides requires making assumptions about the dominant
mechanism driving dissipation in the rocks and ices [@Bagheri+2022aug; RenaudHenning2018apr]. For example,
microscopic grains of ice will tend to move more freely then larger solid crystalline chunks. Likewise, 
rock that has experienced significant fracturing or is porous tends to have more opportunity to create frictional
heat. The choice of **Rheology** determines which dissipation mechanism is dominant within a world. 

![Tidal heating is shown for four different rheology models for a simplified model of Jupiter's moon Io. Heating in certain viscoelastic phase spaces can be orders of magnitude different depending on your choice in rheology.\label{fig:rheology}](io_rheology_comparison.png){ width=50% }

TidalPy provides several different rheological models in its Rheology Module. Most rheologies have empirical 
parameters which are relatively unknown for rocks and ices at planetary temperatures and pressures. TidalPy 
suggests typical values used in the literature but allows you to vary them. These efficient rheological functions can be
used with other TidalPy methods, like the Love number solver, or in your own scripts alongside other tools. 

- Spin-Orbit Dynamics & Coupling - Tidal dissipation transforms the energy from the
rotation and orbit of planets and moons into tidal heat. Complex dynamics can occur when
a planet has a non-synchronous rotation, high eccentricity, and/or has a non-zero obliquity.
TidalPy uses a sophisticated Spin-Orbit Resonance coupling system to account for potential
resonance capture and other long-term dynamics [REF; REF].

![Spin-Orbit Resonance "ledges" calculated with TidalPy. A planet can become trapped on a ledge (stuck at a certain spin rate) for millions of years depending on its interior structure. [adapted from REF].\label{fig:sor}](SOR.png){ width=20% }

TidalPy has proven to be a powerful tool in investigating tides within the Solar
System [REFs] and beyond [REFs]. Future releases will focus on increasing performance,
improving usability, and incorporating more physics. Get started using TidalPy by installing
the package and checking out the documentation over at [https://tidalpy.info](https://tidalpy.info)

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Availability

TidalPy's source code is available and kept up to date on its [GitHub Repository](https://github.com/jrenaud90/TidalPy).
All versions are released on GitHub as well as [PyPI](https://pypi.org/project/TidalPy/) and [Conda-Forge](https://anaconda.org/channels/conda-forge/packages/tidalpy).
Major versions are also released with dedicated DOI's on TidalPy's [Zenodo page](https://zenodo.org/records/10656488).
Anyone is welcome to open pull requests, create forks, or issue bug reports, suggestions, and questions. The latter can be made on the [GitHub issue tracker](https://github.com/jrenaud90/TidalPy/issues).
TidalPy can also be found on NASA's [Exoplanet Modeling and Analysis Center](https://emac.gsfc.nasa.gov/?cid=2207-034) [@emac].

## License

TidalPy is licensed under Creative Commons Attribution-ShareAlike (CC BY-SA). This allows any user to share and reproduce
TidalPy in whole or in part as long as attribution is made back to the original repository and cites this paper. All adapted versions
must carry a similar license (share alike). Full details about the license can be found in the [repository's license file](https://github.com/jrenaud90/TidalPy/blob/main/LICENSE.md).

# Acknowledgements

TidalPy was greatly improved by conversations, contributions, and testing performed by many in the community. We would
like to specifically thank Wade G. Henning, Michael Efroimsky, Michaela Walterová, Sander Goossens, Marc Neveu, Nick Wagner, and Gael Cascioli.
The development of TidalPy was supported by NASA Sellers' Exoplanet Environments Collaboration and Geodesy ISFMs.
J. Renaud was additionally supported during its development by the CRESST-II cooperative agreement (NASA award 80GSFC24M0006).
TidalPy makes extensive use of the following software: CyRK [@cyrk], BurnMan [@burnman], Numpy [@numpy], SciPy [@scipy], Numba [@numba], Matplotlib [@matplotlib], and cmcrameri [@cmcrameri; @scicmap].

# References