# TidalPy Roadmap

* Version 0.1.0 - Basic Usage - June 2019
    * Be able to load a planet and calculate tidal heating for both rocky and icy layers (no ice EOS)
* Version 0.2.0 - Testing Scheme - June 2019
    * Implement a system to quickly run tests so that future patches can be quickly checked to minimize bugs
    * Should work with both Github pull requests and python unittest
* Version 0.3.0 - Orbital Upgrade - July 2019
    * Have both functional and OOP support for orbital dynamics including duel-body dissipation (semi-analytical time domain studies)
* Version 0.4.0 - Ice Upgrade - July 2019
    * Incorporate a realistic Ice EOS (for at least low-pressure and high-pressure cases)
* Version 0.5.0 - Quality of Life Update - July/August 2019
    * Finalize file IO
    * Implement command line argument usage
    * Expand number of cookbooks
    * Mark areas that need documentation for 1.0.0 release
    * Clean-up github issues and add / modify issue categories
    * Expand graphics module
* Version 0.6.0 - Multi-layer Upgrade - Fall 2019
    * Add the S&V multi-layer method
* Version 0.7.0 - Multi-Layer and Orbit Coupling - Fall 2019
    * Couple the multi-layer method to orbital integration (time domain)
* Version 1.0.0 - Release Version
    * Finish Documentation
    * Finalize setup.py and add to pyPI

## Beyond 1.0.0
* Include advanced semi-analytical orbital models such as the Laplace Resonance (use jreanud90/LINT as a baseline).
    
