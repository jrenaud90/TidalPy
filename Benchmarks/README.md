# TidalPy Benchmarks

This directory holds two kinds of benchmarks

- **Validation**: notebooks that show TidalPy can reproduce the physics of published models, closed-form results, and other packages:
  - Love numbers versus Guo et al. (2004) and Farrell (1972), and the closed-form homogeneous viscoelastic Love numbers of Love (1911) with both radial-solver methods;
  - Radial profiles versus Tobie et al. (2005) and Roberts & Nimmo (2008);
  - An equation-of-state comparison against BurnMan;
  - Tidal heating and torque of the constant-time-lag model versus the closed forms of Hut (1981) and Levrard et al. (2007), exact in eccentricity and obliquity, which also measures how far each eccentricity truncation can be trusted;
  - Tidal heating of a viscoelastic body versus an independent calculation from the exact Kepler orbit, for the 1D mode sum and the 3D heating;
  - A reproduction of the high-eccentricity, obliquity, and dual-body results of Renaud et al. (2021).
- **Performance**: a scheme under `Performance/` that times the common tasks used in `Demos/` and stores each run as JSON tagged with the date, TidalPy version, git commit, OS, and CPU, so performance can be tracked over time. See `Performance/README.md`.

## Layout

```
Benchmarks/
  RadialSolver/          validation of the 1D radial (Love number) solver
  EOS/                   equation-of-state comparison against BurnMan
  Tides/                 tidal dissipation, torques, and orbital rates versus closed forms, an exact-orbit
                         calculation, and Renaud et al. (2021)
  Performance/           the performance-tracking harness, tasks, and trend notebook
```

Reference datasets (for example the Earth structure and Love-number files) are stored next to the notebook that uses them.
