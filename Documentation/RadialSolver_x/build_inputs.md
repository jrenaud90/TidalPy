# Helper Functions

_Updated: 2026-09-12_

`TidalPy.RadialSolver_x.radial_solver` takes array-based inputs: a radius grid, the density and complex moduli on that grid, the planet bulk density, and a few per-layer descriptors. Two native builders assemble those inputs from a layer description so you do not have to hand-build the arrays:

| Builder | Use when |
|---|---|
| `build_rs_input_homogeneous_layers` | Each layer has constant density, moduli, and viscosities. |
| `build_rs_input_from_data` | You already have radially resolved data (for example from an EOS solver or a published interior model). |

Both return a `PlanetBuildData` named tuple whose fields are, in order, the positional arguments of `radial_solver`, so `radial_solver(*build_data, **solver_kwargs)` runs the solve.

The builders live in C++ (`RadialSolver_x/build_inputs_.hpp`) behind a thin Cython layer and evaluate the complex moduli with the `TidalPy.rheology_x` models.


> [!TIP]
> These helper functions have been built to be very efficient, however they will still cost some performance
> overhead when used. If calculation speed is critical, and you are rebuilding a planet many times (e.g., in
> a MCMC) it may be more performant to manually construct RadialSolver's inputs and only change what is needed
> rather than calling these helpers every time.

## Rheology Arguments

For each of `shear_rheology_model_tuple` and `bulk_rheology_model_tuple` you can pass:

- One `rheology_x` model instance, applied to every layer;
- One model name accepted by `make_rheology` (for example `"maxwell"`, `"elastic"`, `"andrade"`);
- A tuple or list with one instance or name per layer.

```python
from TidalPy.rheology_x import Andrade, Elastic, Maxwell

# Same shear rheology in every layer, elastic bulk response everywhere.
shear_rheology_model_tuple = Maxwell()
bulk_rheology_model_tuple = "elastic"

# Or one model per layer (here: solid core, liquid ocean, solid crust).
shear_rheology_model_tuple = (Maxwell(), Elastic(), Andrade(alpha=0.3, zeta=1.0))
bulk_rheology_model_tuple = (Elastic(), Elastic(), Elastic())
```

Models with parameters (Andrade, Sundberg, Burgers, Voigt) keep whatever parameters they were built with.

## Planet with Homogeneous Layers

If your planet of interest has multiple distinct layers, but you assume each layer is homogenous in composition and viscoelastic properties (and density), then this helper can build rs inputs with a minimal amount of information from the user.

```python
import numpy as np
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.rheology_x import Elastic, Maxwell

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

Layer sizes are given by exactly one of:

| Argument | Meaning |
|---|---|
| `thickness_fraction_tuple` | Each layer's thickness over the planet radius (sums to 1). |
| `radius_fraction_tuple` | Each layer's upper radius over the planet radius (increasing, last entry 1). |
| `volume_fraction_tuple` | Each layer's share of the planet volume (sums to 1). |

Each layer's grid runs from its base to its top (both inclusive) with `slices_tuple[i]` (or `slice_per_layer`) evenly spaced slices, so interface radii appear twice in `radius_array` as the solver requires. Every layer needs at least 5 slices. The liquid layer above has zero static shear modulus and an elastic shear rheology, which gives it exactly zero complex shear modulus.

## Planet from Radially Resolved Data Arrays

Useful if utilizing a 3rd party EOS solver or when utilizing data found in the literature.

```python
import numpy as np
from TidalPy.RadialSolver_x import build_rs_input_from_data, radial_solver
from TidalPy.rheology_x import Elastic, Maxwell

# Data from a third-party interior model, for example loaded from disk. 
# This must be provided as MKS units to the helper.
data = np.load("my_interior_model.npz")
radius_array          = data["radius"]
density_array         = data["density"]
static_bulk_array     = data["bulk_modulus"]   # This is the _static_ bulk modulus
static_shear_array    = data["shear_modulus"]  # This is the _static_ shear modulus
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
    bulk_rheology_model_tuple     = Elastic(),  # Same for all three layers in this case
    warnings                      = True,
)

solution = radial_solver(*build_data, degree_l=2)
```

The solver requires that the grid starts at `r = 0`, that every interface radius appears twice (top of the lower layer and base of the upper layer), and that each layer's upper radius is a grid point. The builder copies your arrays and repairs them where needed, giving each inserted slice the properties of the neighbouring provided slice: an inserted layer base copies the layer's first provided slice, an inserted layer top copies the slice below it. An interface radius that appears only once is taken as the top of the lower layer, with the properties it carries, and the upper layer's base is inserted above it (the same convention as the classic builder). Each repair is logged as a warning through the TidalPy C++ logger (see [logging](../utilities_x/index.md)); pass `warnings=False` to silence them. The planet bulk density is the mass of the piecewise-constant shells divided by the planet volume.

Array arguments accept anything `numpy.asarray` understands (lists included); they are converted to contiguous float64 arrays.

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

- `TypeError`: a rheology argument is not a `rheology_x` model instance or model name (classic
`TidalPy.rheology` models included).
- `ValueError`: per-layer inputs with the wrong length, fractions that do not describe the whole planet,
fewer than 5 slices in a layer, more (or fewer) than one layer-size description, a non-ascending radius grid, or a last layer upper radius that is not the planet radius.

`perform_checks` is accepted for signature compatibility with the classic builders; the native builders always validate their inputs.

## A uniform sphere in one call: `homogeneous_love_numbers`

When the interior does not matter, `homogeneous_love_numbers` builds the arrays for a single uniform solid layer and runs the solve for you. It is the quickest way to a Love number for a demo, a benchmark, or a sanity check against the closed-form result.

```python
from TidalPy.RadialSolver_x import homogeneous_love_numbers
from TidalPy.rheology_x import Maxwell

mu = Maxwell().calc_complex_modulus(60.0e9, 1.0e15, 4.1e-5)   # complex shear modulus [Pa]
solution = homogeneous_love_numbers(1.8216e6, 3529.0, mu, 4.1e-5)
print(solution.k)
```

| Argument | Default | Meaning |
|---|---|---|
| `planet_radius` | required | Planet radius [m]. |
| `planet_bulk_density` | required | Uniform density [kg m-3]. |
| `complex_shear_modulus` | required | Complex shear modulus at the forcing frequency [Pa], usually from a rheology model's `calc_complex_modulus`. |
| `forcing_frequency` | required | Tidal forcing frequency [rad s-1]. |
| `complex_bulk_modulus` | `200e9 + 0j` | Complex bulk modulus [Pa]. |
| `num_slices` | `60` | Slices in the generated grid. |
| `degree_l` | `2` | Harmonic degree. |
| `layer_is_static`, `layer_is_incompressible` | `True`, `False` | Layer assumptions. |
| `**radial_solver_kwargs` | | Anything else goes straight to `radial_solver`, for example `solve_for`, `love_method`, or the integration tolerances. |

It returns the same [`RadialSolverSolution`](solution_class.md) as any other solve. For the closed-form answer without any integration at all, use `calc_homogeneous_love_numbers` in [`TidalPy.Tides_x.love`](../Tides_x/love/love_numbers.md); for a static incompressible sphere the two agree to machine precision (measured at a few parts in 1e15).

## Migrating from the classic builders

| Classic (`TidalPy.RadialSolver.helpers`) | New (`TidalPy.RadialSolver_x`) |
|---|---|
| `from TidalPy.rheology.models import Maxwell` | `from TidalPy.rheology_x import Maxwell` |
| `shear_rheology_model_tuple=(Maxwell(), Maxwell())` | `shear_rheology_model_tuple=Maxwell()` (or the tuple) |
| Warnings through the Python `TidalPy` logger | Warnings through the C++ logger (same `[logging]` config) |
| Errors raise `ArgumentException` | Errors raise `ValueError` / `TypeError` |
