# The Future TidalPy Structure

TidalPy's internals are being rewritten in C++. The new implementation lives in modules that carry a `_x` suffix
(`structures_x`, `Tides_x`, `RadialSolver_x`, `rheology_x`, and so on) and ships side by side with the classic
modules today. In a future major release the `_x` modules will become the only TidalPy: the classic modules will be
removed and the new ones will drop their suffix. Nothing about the classic API changes until then, but any new
development on TidalPy will happen in the `_x` modules, other than bug fixes. We highly encourage new projects
to start using the `_x` modules or make plans to switch.

This page explains what is different, maps the classic modules to their replacements, and shows how to port common
workflows.

TidalPy announces this transition once per session when the package is imported. The notice can be silenced with:

```python
import warnings
from TidalPy.exceptions import TidalPyDeprecationWarning
warnings.filterwarnings("ignore", category=TidalPyDeprecationWarning)
```

## Why the rework?

* **Performance.** All core physics now runs in C++ (with the Eigen linear algebra library and CyRK integrators),
  wrapped by thin Cython layers. There is no numba JIT warmup, and hot paths avoid Python entirely.
* **Predictability.** The classic system stored planet state on Python objects and propagated changes through
  cascading updates, which was hard to reason about and easy to break. The new classes store configuration and
  return results from `calc_*` methods without mutating states.
* **Consistency.** Every physics module (rheology, cooling, radiogenics, viscosity, partial melting, equations of
  state, tides) follows the pattern: a C++ class hierarchy, a name-based factory (`make_<module>`), direct callable
  functions, vectorized variants, TOML configuration, and binary files which can be saved and loaded from disk for
  fast and accurate reproducibility.

## Module map

| Classic module | Replacement | Documentation |
|----------------|-------------|---------------|
| `TidalPy.structures` (worlds, layers, burnman builds) | `TidalPy.structures_x` | [Worlds](structures_x/worlds/worlds.md), [Layers](structures_x/layers/base_layer.md), [System](structures_x/system/system.md), [TOML schema](structures_x/config/toml_schema.md) |
| `TidalPy.RadialSolver` | `TidalPy.RadialSolver_x` | [Love number solvers](RadialSolver_x/dense_radial_solution.md) |
| `TidalPy.rheology` | `TidalPy.rheology_x` | [Rheology models](rheology_x/rheology_models.md) |
| `TidalPy.tides` | `TidalPy.Tides_x` (mostly via world methods) | [Global tides](Tides_x/global_tides.md), [3D heating](Tides_x/multilayer_3d_heating.md), [Love numbers](Tides_x/love/love_numbers.md) |
| `TidalPy.cooling` | `TidalPy.cooling_x` | [Cooling models](cooling_x/cooling_models.md) |
| `TidalPy.radiogenics` | `TidalPy.radiogenics_x` | [Radiogenics models](radiogenics_x/radiogenics_models.md) |
| viscosity functions in `TidalPy.rheology` | `TidalPy.viscosity_x` | [Viscosity models](viscosity_x/viscosity_models.md) |
| partial melting in `TidalPy.rheology` | `TidalPy.partial_melt_x` | [Partial melt models](partial_melt_x/partial_melt_models.md) |
| burnman equations of state | `TidalPy.Material_x` | [Material EOS](material_x/material_eos.md) |
| `TidalPy.stellar` | `TidalPy.stellar_x` (world-attached luminosity) | [Luminosity](stellar_x/luminosity.md) |
| orbit/spin evolution in `TidalPy.orbit` | `TidalPy.dynamics_x` + the `System` class | [Dynamics](dynamics_x/dynamics.md), [System](structures_x/system/system.md) |
| assorted helpers | `TidalPy.Utilities_x` | [Utilities](utilities_x/index.md) |

## Porting examples

### Building a world and getting Love numbers

Classic (burnman-backed OOP world):

```python
from TidalPy.structures import build_world

world = build_world("earth")
# State-based updates followed: world.orbit = ..., world.update_orbit_spin(), etc.
```

New system (TOML-configured world; solves are explicit method calls):

```python
from TidalPy.structures_x import build_world

world = build_world("earth_simple")
world.solve_eos()

result = world.solve_love_numbers(frequency_rad_s=1.0e-5)
print(result["success"], world.love_number_k)
```

The world builder reads bundled or user TOML files, a file path, or a Python dictionary. See the
[TOML schema](structures_x/config/toml_schema.md) and the `Demos (_x)` notebooks for more details.

### The standalone radial solver

The array-based `radial_solver` call is nearly identical between the two systems:

```python
# Classic
from TidalPy.RadialSolver import radial_solver

# New (same call signature, same helper builders, new dense-output solution class)
from TidalPy.RadialSolver_x import radial_solver
```

The new solver adds dense radial evaluation at arbitrary radii (`get_radial_solution`), implicit CyRK integrators
(`BDF`, `LSODA`, `Radau`), and a surface conditioning diagnostic (`surface_solve_amplification`). See
[dense radial solutions](RadialSolver_x/dense_radial_solution.md).

### Rheology models

```python
# Classic
from TidalPy.rheology.models import Maxwell
rheology = Maxwell()
complex_modulus = rheology(frequency, shear_modulus, viscosity)

# New (note the argument order: modulus, viscosity, frequency)
from TidalPy.rheology_x import Maxwell, make_rheology
rheology = Maxwell()
complex_modulus = rheology.calc_complex_modulus(shear_modulus, viscosity, frequency)

# Or build by name (case-insensitive, alias-aware)
rheology = make_rheology("maxwell")
```

Every model also provides vectorized `calc_*` variants that accept numpy arrays and broadcast.

## Learning the new system

The best introduction is the `Demos (_x)` notebooks in the navigation: `Basics` notebooks
(configuration, world building, save/load), `Physics` notebooks (orbits, tides, rheology, Love numbers,
3D heating, thermal/EOS), and `Systems` notebooks (multi-world systems, coupled thermal-orbital evolution).
The `Benchmarks (_x)` pages validate the new system against published results and track its performance.
