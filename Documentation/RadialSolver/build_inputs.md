# Helper Functions

_Updated: 2026-09-29_

`radial_solver` takes array inputs: a radius grid, the density and complex moduli on it, the planet bulk density, and per-layer assumptions. Two builders assemble these from a layer description, evaluating the complex moduli with `TidalPy.Rheology` models.

| Builder | Use when |
|---|---|
| `build_rs_input_homogeneous_layers` | Each layer has constant density, moduli, and viscosities. |
| `build_rs_input_from_data` | You have radially resolved data, from an EOS solver or a published interior model. |

Both return a `PlanetBuildData` named tuple whose fields are, in order, the positional arguments of `radial_solver`, so `radial_solver(*build_data, **solver_kwargs)` runs the solve.

> [!TIP]
> If you rebuild a planet many times (in an MCMC, say), build the inputs once and change only what is needed, or better, use the world-attached solver, which allocates far less memory.

## Rheology Arguments

`shear_rheology_model_tuple` and `bulk_rheology_model_tuple` each take:

- one `Rheology` model instance, applied to every layer;
- one model name accepted by `make_rheology` (for example `"maxwell"`, `"elastic"`, `"andrade"`);
- a tuple or list with one instance or name per layer.

```python
from TidalPy.Rheology import Andrade, Elastic, Maxwell

shear_rheology_model_tuple = Maxwell()   # Every layer
bulk_rheology_model_tuple = "elastic"

# One model per layer: solid core, liquid ocean, solid crust
shear_rheology_model_tuple = (Maxwell(), Elastic(), Andrade(alpha=0.3, zeta=1.0))
bulk_rheology_model_tuple = (Elastic(), Elastic(), Elastic())
```

Models with parameters (Andrade, Sundberg, Burgers, Voigt) keep the parameters they were built with.

## Planet with Homogeneous Layers

```python
import numpy as np
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Elastic, Maxwell

# Three-layer planet: solid inner core, liquid outer core, solid mantle (MKS units).
planet_radius     = 6000.0e3
forcing_frequency = 2.0 * np.pi / (86400.0 * 7.5)

build_data = build_rs_input_homogeneous_layers(
    planet_radius,
    forcing_frequency,
    density_tuple                 = (8000.0, 5000.0, 3300.0),
    static_bulk_modulus_tuple     = (2.0e11, 1.5e11, 1.0e11),
    static_shear_modulus_tuple    = (1.0e11, 0.0, 5.0e10),
    bulk_viscosity_tuple          = (1.0e18, 1.0e18, 1.0e18),
    shear_viscosity_tuple         = (1.0e20, 1.0e3, 1.0e19),
    layer_type_tuple              = ("solid", "liquid", "solid"),
    layer_is_static_tuple         = (False, True, False),
    layer_is_incompressible_tuple = (False, False, False),
    shear_rheology_model_tuple    = (Maxwell(), Elastic(), Maxwell()),
    bulk_rheology_model_tuple     = Elastic(),          # one model for every layer
    radius_fraction_tuple         = (0.2, 0.55, 1.0),   # layer tops as fractions of the planet radius
    slices_tuple                  = (10, 12, 20),       # or slice_per_layer=10 for all layers
)

solution = radial_solver(*build_data, degree_l=2, solve_for=("tidal",))
print(solution.k, solution.h, solution.l)
```

The liquid layer's zero static shear modulus and elastic shear rheology give it exactly zero complex shear modulus.

Give the layer sizes with exactly one of:

| Argument | Meaning |
|---|---|
| `thickness_fraction_tuple` | Each layer's thickness over the planet radius (sums to 1). |
| `radius_fraction_tuple` | Each layer's upper radius over the planet radius (increasing, last entry 1). |
| `volume_fraction_tuple` | Each layer's share of the planet volume (sums to 1). |

Each layer's grid runs from its base to its top inclusive, with `slices_tuple[i]` (or `slice_per_layer`) evenly spaced slices, at least 5, so interface radii appear twice as the solver requires.

## Planet from Radially Resolved Data Arrays

```python
import numpy as np
from TidalPy.RadialSolver import build_rs_input_from_data, radial_solver
from TidalPy.Rheology import Elastic, Maxwell

# Data from a third-party interior model, in MKS units.
data = np.load("my_interior_model.npz")
radius_array          = data["radius"]
density_array         = data["density"]
static_bulk_array     = data["bulk_modulus"]
static_shear_array    = data["shear_modulus"]
shear_viscosity_array = data["viscosity"]
bulk_viscosity_array  = np.full_like(radius_array, 1.0e18)   # unused by an elastic bulk rheology

forcing_frequency = 2.0 * np.pi / (86400.0 * 1.5)

build_data = build_rs_input_from_data(
    forcing_frequency,
    radius_array,
    density_array,
    static_bulk_array,
    static_shear_array,
    bulk_viscosity_array,
    shear_viscosity_array,
    layer_upper_radius_tuple      = (1.2e6, 1.8e6, radius_array[-1]),  # last entry: planet radius
    layer_type_tuple              = ("solid", "liquid", "solid"),
    layer_is_static_tuple         = (False, True, False),
    layer_is_incompressible_tuple = (False, False, False),
    shear_rheology_model_tuple    = (Maxwell(), Elastic(), Maxwell()),
    bulk_rheology_model_tuple     = Elastic(),
    warnings                      = True,
)

solution = radial_solver(*build_data, degree_l=2)
```

Array arguments accept anything `numpy.asarray` understands, lists included. The builder copies them and repairs the grid where the solver's rules (start at `r = 0`, every interface radius twice, each layer's upper radius a grid point) are not met, logging a warning per repair (`warnings=False` silences them). An inserted layer base copies the layer's first provided slice, an inserted layer top the slice below it. An interface radius given once is the top of the lower layer, and the upper layer's base is inserted above it. The planet bulk density is the mass of the piecewise-constant shells over the planet volume.

## Output

`PlanetBuildData` fields, in `radial_solver` order:

| Field | Type |
|---|---|
| `radius_array` | float64 array [m] |
| `density_array` | float64 array [kg m⁻³] |
| `complex_bulk_modulus_array` | complex128 array [Pa] |
| `complex_shear_modulus_array` | complex128 array [Pa] |
| `frequency` | float [rad s⁻¹] |
| `planet_bulk_density` | float [kg m⁻³] |
| `layer_types` | tuple of str |
| `is_static_bylayer` | tuple of bool |
| `is_incompressible_bylayer` | tuple of bool |
| `upper_radius_bylayer_array` | float64 array [m] |

## Errors

- `TypeError`: a rheology argument is not a `Rheology` model instance or name.
- `ValueError`: per-layer inputs of the wrong length, fractions that do not describe the whole planet, fewer than 5 slices in a layer, not exactly one layer-size description, a non-ascending radius grid, or a last layer upper radius that is not the planet radius.

## Uniform Sphere: `homogeneous_love_numbers`

`homogeneous_love_numbers` builds the arrays for a single uniform solid layer and runs the solve: the quickest Love number for a demo, a benchmark, or a sanity check.

```python
from TidalPy.RadialSolver import homogeneous_love_numbers
from TidalPy.Rheology import Maxwell

mu = Maxwell().calc_complex_modulus(60.0e9, 1.0e15, 4.1e-5)   # complex shear modulus [Pa]
solution = homogeneous_love_numbers(1.8216e6, 3529.0, mu, 4.1e-5)
print(solution.k)
```

| Argument | Default | Meaning |
|---|---|---|
| `planet_radius` | required | Planet radius [m]. |
| `planet_bulk_density` | required | Uniform density [kg m-3]. |
| `complex_shear_modulus` | required | Complex shear modulus at the forcing frequency [Pa], usually from a rheology's `calc_complex_modulus`. |
| `forcing_frequency` | required | Tidal forcing frequency [rad s-1]. |
| `complex_bulk_modulus` | `200e9 + 0j` | Complex bulk modulus [Pa]. |
| `num_slices` | `60` | Slices in the generated grid. |
| `degree_l` | `2` | Harmonic degree. |
| `layer_is_static`, `layer_is_incompressible` | `True`, `False` | Layer assumptions. |
| `**radial_solver_kwargs` | | Passed to `radial_solver`, for example `solve_for`, `love_method`, or the tolerances. |

It returns a [`RadialSolverSolution`](solution_class.md). For a closed-form answer, use `calc_homogeneous_love_numbers` in [`TidalPy.Tides.love`](../Tides/love/love_numbers.md); the two agreed for a static incompressible sphere in testing.
