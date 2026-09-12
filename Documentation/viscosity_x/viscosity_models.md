# Viscosity Models (`viscosity_x`)

_Updated: 2026-09-12_

A viscosity model returns a material's dynamic viscosity $\eta$ [Pa s] as a function of temperature [K] and pressure [Pa]. This is the pre-melt, or "solid", viscosity: the value a material would show with no melt present, which the [partial-melt](../partial_melt_x/partial_melt_models.md) step then weakens. Both are frequency-independent, so both are resolved once per equation-of-state solve and reused across every tidal forcing frequency.

All three models share the same physical picture. Solid-state creep is thermally activated, so viscosity falls exponentially with temperature through a Boltzmann factor $\exp(E_a / RT)$, and rises with pressure through an activation volume that makes the same creep harder to accommodate at depth. The models differ only in how that exponential is anchored: to an absolute flow-law prefactor, to a measured reference viscosity, or not at all.

Quantities are MKS throughout. The math mirrors the validated classic implementation in `TidalPy/rheology/viscosity/viscosity_models.py`, and the molar gas constant $R$ comes from the shared TidalPy configuration rather than being hard-coded.

## The three models

| Model | Aliases | Viscosity $\eta$ [Pa s] |
|---|---|---|
| `ArrheniusViscosity` | `arrhenius`, `arr` | $A \, \sigma^{1-n} d^{\,m} \exp\left(\dfrac{E_a + P V_a}{R T}\right)$, multiplied by $T$ when `additional_temp_dependence` is set |
| `ReferenceViscosity` | `reference`, `ref` | $\eta_\mathrm{ref} \exp\left[\dfrac{E_a + P V_a}{R} \left(\dfrac{1}{T} - \dfrac{1}{T_\mathrm{ref}}\right)\right]$ |
| `ConstantViscosity` | `constant`, `const` | $\eta_\mathrm{ref}$, independent of temperature and pressure |

| Parameter | Symbol | Default | Units | Used by |
|---|---|---|---|---|
| `reference_viscosity` | $\eta_\mathrm{ref}$ | 1.0e22 | Pa s | Reference, Constant |
| `reference_temperature` | $T_\mathrm{ref}$ | 1000.0 | K | Reference |
| `molar_activation_energy` | $E_a$ | 3.0e5 | J mol$^{-1}$ | Arrhenius, Reference |
| `molar_activation_volume` | $V_a$ | 0.0 | m$^3$ mol$^{-1}$ | Arrhenius, Reference |
| `arrhenius_coeff` | $A$ | 1.0 | model-dependent | Arrhenius |
| `stress` | $\sigma$ | 1.0 | Pa | Arrhenius |
| `stress_expo` | $n$ | 1.0 | — | Arrhenius |
| `grain_size` | $d$ | 1.0e-3 | m | Arrhenius |
| `grain_size_expo` | $m$ | 0.0 | — | Arrhenius |
| `additional_temp_dependence` | — | `False` | — | Arrhenius |

The stress exponent $n$ distinguishes creep regimes: $n = 1$ is diffusion creep, where the flow law is linear and the stress term drops out, and $n > 1$ is dislocation creep, where the material shears more readily the harder it is pushed. The grain-size exponent plays the same role for grain-boundary processes. Setting `additional_temp_dependence` adds the explicit factor of $T$ that some published diffusion-creep flow laws carry in front of the exponential.

### Behavior at the limits

Every model returns infinity at or below zero temperature. The cold limit of a thermally activated fluid is a solid that does not flow, and an infinite viscosity is exactly what the rheology models need to return a purely elastic response there. A reference model with a non-positive reference temperature returns infinity for the same reason. Very cold but positive temperatures reach the same place by overflowing the exponential, which is deliberate rather than guarded against.

## Choosing a model

`ReferenceViscosity` is the usual choice for a planetary interior. It is anchored to a viscosity someone actually measured or inferred at a stated temperature, so its one strong assumption, the activation energy, is the only thing you have to defend. It is also the easiest to reason about: doubling $E_a$ visibly steepens the response to temperature without moving the anchor point.

`ArrheniusViscosity` is the choice when you are reproducing a published flow law, which is usually quoted as an absolute prefactor with stress and grain-size exponents. It gives full control at the cost of needing every term to be right, since nothing pins the result to a known viscosity.

`ConstantViscosity` is for tests, for benchmark comparisons against analytic results, and for layers whose temperature you have no reason to trust. It is not a physical claim, and using it removes the feedback between heating and viscosity that makes tidal evolution interesting.

## Python API

```python
from TidalPy.viscosity_x import (
    ArrheniusViscosity, ConstantViscosity, ReferenceViscosity, make_viscosity)

viscosity_model = ReferenceViscosity(
    reference_viscosity=1.0e22, reference_temperature=1000.0,
    molar_activation_energy=3.0e5, molar_activation_volume=0.0)

eta = viscosity_model.calc_viscosity(temperature=1500.0, pressure=1.0e10)   # Pa s

# Name factory: case-insensitive, aliases accepted.
arrhenius_model = make_viscosity("arr", {"arrhenius_coeff": 1.1e7,
                                         "grain_size_expo": 2.0,
                                         "additional_temp_dependence": True})
```

Constructors take every parameter their model uses as a keyword with the default from the table above: `ConstantViscosity(reference_viscosity=1.0e22)`, `ReferenceViscosity(reference_viscosity=1.0e22, reference_temperature=1000.0, molar_activation_energy=3.0e5, molar_activation_volume=0.0)`, and `ArrheniusViscosity(arrhenius_coeff=1.0, stress=1.0, stress_expo=1.0, grain_size=1.0e-3, grain_size_expo=0.0, molar_activation_energy=3.0e5, molar_activation_volume=0.0, additional_temp_dependence=False)`.

`make_viscosity(model_name, config=None)` resolves a name or alias case-insensitively and reads the parameters it recognizes from `config`. Absent keys fall back to the model default, keys a model does not use are ignored, and an unrecognized name raises `ValueError`.

| Member | Returns | Description |
|---|---|---|
| `calc_viscosity(temperature, pressure=0.0)` | `float` [Pa s] | Dynamic viscosity at those conditions. |
| `model_name` | `str` | The resolved model name (`arrhenius`, `reference`, `constant`). |
| `get_config_dict()` | `dict` | `model` plus every parameter the model carries. |
| `save_config(path)` | — | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path, force=False)` | — | TidalPy binary format; see [Binary serialization](../utilities_x/binary_x.md). |

Parameters are read-only properties: `reference_viscosity` on the constant model; `reference_viscosity`, `reference_temperature`, `molar_activation_energy`, and `molar_activation_volume` on the reference model; and `arrhenius_coeff`, `molar_activation_energy`, and `additional_temp_dependence` on the Arrhenius model. Every parameter, including the ones without a property, appears in `get_config_dict()`.

## Attaching a viscosity model to a layer

```python
from TidalPy.viscosity_x import make_viscosity
from TidalPy.structures_x.layers.physics import PhysicsLayer

mantle = PhysicsLayer("mantle", 0, 0.0, 1.0e6, 2.1e19,
                      shear_modulus_static=50.0e9, bulk_modulus_static=100.0e9)

mantle.set_shear_viscosity(make_viscosity("reference", {"reference_viscosity": 1.0e20}))
mantle.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e20}))
```

Ownership of the C++ model transfers into the layer, exactly as it does for a rheology. During the world's equation-of-state solve each radial slice's temperature and pressure are pushed through the model, an equation of state that supplies its own viscosity profile overrides the result slice by slice, and the partial-melt model then weakens what remains. Read the outcome back with the layer's `get_shear_viscosity(radius)` and `get_premelt_shear_viscosity(radius)`. The declarative form is a `[layers.<name>.shear_viscosity]` table in a world's TOML; see the [TOML schema](../structures_x/config/toml_schema.md).

## C++ API

The C++ layer is canonical and the Cython classes are thin adapters over it.

`c_ViscosityConfig` holds every parameter for every model in one struct, with the defaults listed above. `c_ViscosityBase : c_PhysicsBase` (in `viscosity_base_.hpp`) declares `calc_viscosity(double temperature, double pressure) const` pure virtual and adds `calc_viscosity_vectorize(temperature, pressure, out_viscosity)`, the radial sweep that fills one entry per slice and is the path the equation-of-state solve actually uses. The concrete models `c_ArrheniusViscosity`, `c_ReferenceViscosity`, and `c_ConstantViscosity` live in `viscosity_.hpp`, each with a getter per parameter, and occupy binary class ids 801 through 803.

The factory mirrors the rheology one: `c_viscosity_model_from_name(name)` maps a name or alias onto the `c_ViscosityModel` enum and throws `std::invalid_argument` for an unknown name, `c_find_viscosity(model, config)` returns a `std::unique_ptr<c_ViscosityBase>` (a name overload does both steps), and `c_viscosity_from_binary(stream, force=false)` peeks the class id, builds the matching model, and calls `read_binary`.

## Adding a new model

1. Add the model's parameters to `c_ViscosityConfig` in `viscosity_.hpp`, with defaults.
2. Add `c_<Name>Viscosity : c_ViscosityBase` implementing `calc_viscosity`, `append_config_entries`, `write_binary`, and `read_binary`.
3. Reserve the next free `BinaryClassID` in the 800 block in `Utilities_x/binary_x/binary_.hpp`.
4. Register a `c_ViscosityModel::<Name>` enum value and wire it into `c_viscosity_model_from_name`, `c_find_viscosity`, and `c_viscosity_from_binary`.
5. Add the Cython `cdef class` with its parameter properties, the adoption branch in `make_viscosity`, tests in `Tests/Test_Viscosity_x/`, and an entry on this page.

## References

- Moore, W. B. (2006). Thermal equilibrium in Europa's ice shell. *Icarus*, 180(1), 141-146. Arrhenius flow law with activation energy and volume.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015. Reference-viscosity (relative activation) law.
