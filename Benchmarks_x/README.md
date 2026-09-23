# TidalPy Benchmarks

This directory holds two kinds of benchmark:

- **Validation**: notebooks that show TidalPy reproduces the physics of published models and other packages (Love numbers versus Guo et al. 2004 and Farrell 1972, radial profiles versus Tobie et al. 2005 and Roberts & Nimmo 2008, an equation-of-state comparison against BurnMan, and a reproduction of the high-eccentricity, obliquity, and dual-body results of Renaud et al. 2021).
- **Performance**: a harness under `Performance/` that times the common tasks taught in `Demos_x/` and stores each run as JSON tagged with the date, TidalPy version, git commit, OS, and CPU, so performance can be tracked over time. See `Performance/README.md`.

## Layout

```
Benchmarks_x/
  RadialSolver/          validation of the 1D radial (Love number) solver
  EOS/                   equation-of-state comparison against BurnMan
  Tides/                 tidal dissipation and orbital rates versus Renaud et al. (2021)
  Performance/           the performance-tracking harness, tasks, and trend notebook
```

Reference datasets (for example the Earth structure and Love-number files) are stored next to the notebook that uses them.
