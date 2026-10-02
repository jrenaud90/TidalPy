# Viscosity Models (`Viscosity`)

_Updated: 2026-10-01_

A viscosity model returns a material's dynamic viscosity $\eta$ \[Pa s\] as a function of temperature \[K\] and pressure \[Pa\]. This is the pre-melt viscosity of one phase: the value the phase would show with no melt present, which the material's [melt weakening](../PartialMelt/partial_melt_models.md#melt-weakening) then lowers. Both are frequency-independent, so both are resolved once per equation-of-state solve and reused across every tidal forcing frequency.

Solid-state creep is thermally activated. Viscosity falls exponentially with temperature through the Boltzmann factor $\exp(E_a / RT)$ and rises with pressure through an activation volume. The models differ in how that exponential is anchored: to an absolute flow-law prefactor, to a reference viscosity at a reference temperature, or not at all.

## Models

| Model | Aliases | Viscosity $\eta$ \[Pa s\] |
|---|---|---|
| `ArrheniusViscosity` | `arrhenius`, `arr` | $A \, \sigma^{1-n} d^{\,m} \exp\left(\dfrac{E_a + P V_a}{R T}\right)$, multiplied by $T$ when `additional_temp_dependence = True` |
| `ReferenceViscosity` | `reference`, `ref` | $\eta_\mathrm{ref} \exp\left[\dfrac{E_a + P V_a}{R T} - \dfrac{E_a + P_\mathrm{ref} V_a}{R T_\mathrm{ref}}\right]$ |
| `ConstantViscosity` | `constant`, `const` | $\eta_\mathrm{ref}$, independent of temperature and pressure |
| `InterpolatedViscosity` | `interpolate`, `interp`, `interpolated` | Linear in radius between the points of a `radius_m` and `viscosity_pas` table, held at the end values beyond it |
| `CompositeViscosity` | `composite`, `parallel` | $\left(\sum_i 1 / \eta_i\right)^{-1}$ over its `mechanisms`, each a viscosity model: deformation mechanisms acting in parallel at a common stress, so the weakest dominates |

Each parameter carries two names: the constructor keyword, which also reads as an attribute, and the config key used in a TOML table, a `make_viscosity` config dictionary, and `get_config_dict()`. A dimensional config key ends in its unit while the code name does not.

| Parameter | Config key | Symbol | Default | Units | Used by |
|---|---|---|---|---|---|
| `reference_viscosity` | `reference_viscosity_pas` | $\eta_\mathrm{ref}$ | 1.0e22 | Pa s | Reference, Constant |
| `reference_temperature` | `reference_temperature_k` | $T_\mathrm{ref}$ | 1000.0 | K | Reference |
| `molar_activation_energy` | `molar_activation_energy_j_mol` | $E_a$ | 3.0e5 | J mol$^{-1}$ | Arrhenius, Reference |
| `molar_activation_volume` | `molar_activation_volume_m3_mol` | $V_a$ | 0.0 | m$^3$ mol$^{-1}$ | Arrhenius, Reference |
| `reference_pressure` | `reference_pressure_pa` | $P_\mathrm{ref}$ | 0.0 | Pa | Reference |
| `radius` | `radius_m` | $r$ | `[0.0]` | m | Interpolated |
| `viscosity` | `viscosity_pas` | $\eta$ | `[1.0e22]` | Pa s | Interpolated |
| `arrhenius_coeff` | `arrhenius_coeff` | $A$ | 1.0 | model-dependent | Arrhenius |
| `stress` | `stress_pa` | $\sigma$ | 1.0 | Pa | Arrhenius |
| `stress_expo` | `stress_expo` | $n$ | 1.0 | - | Arrhenius |
| `grain_size` | `grain_size_m` | $d$ | 1.0e-3 | m | Arrhenius |
| `grain_size_expo` | `grain_size_expo` | $m$ | 0.0 | - | Arrhenius |
| `additional_temp_dependence` | `additional_temp_dependence` | - | `False` | - | Arrhenius |

The stress exponent $n$ distinguishes creep regimes. With $n = 1$ the material is in diffusion creep, where the flow law is linear and the stress term drops out, and with $n > 1$ it is in dislocation creep, where the viscosity decreases as the stress increases. The grain-size exponent plays the same role for grain-boundary processes. Setting `additional_temp_dependence = True` adds the explicit factor of $T$ that some published diffusion-creep flow laws carry in front of the exponential.

### Behavior at the Limits

Published flow laws for ice and olivine combine several mechanisms (diffusion creep, dislocation creep, grain-boundary sliding), and a law whose activation energy switches at a temperature is two mechanisms; `CompositeViscosity(diffusion_law, dislocation_law)` or a TOML table `model = "composite"` with `mechanisms = [{ model = ... }, ...]` combines them. A `ConstantViscosity` may be infinite, a rigid (purely elastic) material.

The reference viscosity $\eta_\mathrm{ref}$ is the viscosity at $T_\mathrm{ref}$ and $P_\mathrm{ref}$. With the default $P_\mathrm{ref} = 0$ a positive activation volume always raises the viscosity with pressure, as it does in the Arrhenius law; a law for a deep layer with a large activation volume is better anchored at a pressure inside the layer. The interpolated model reads only the radius, so its layer must keep its volume. A seismic profile's quality factor travels in this table for the `seismic_q` rheology, which reads its viscosity input as $Q$.

The Arrhenius and reference models return infinity at or below zero temperature (`ConstantViscosity` returns its constant at any temperature, `InterpolatedViscosity` its table value). The cold limit of a thermally activated fluid is a solid that does not flow, and an infinite viscosity makes the rheology models return a purely elastic response. A non-positive reference temperature is refused when the model is built. Very cold, but positive, temperatures can also return infinity when the exponential overflows. This is intended behavior.

### Choosing a Model

`ReferenceViscosity` is the usual choice for a planetary interior. It is anchored to a viscosity measured or inferred at a stated temperature, and its one strong assumption is the activation energy. Doubling $E_a$ steepens the response to temperature without moving the anchor point.

`ArrheniusViscosity` is the choice when reproducing a published flow law, which is usually quoted as an absolute prefactor with stress and grain-size exponents. It gives full control at the cost of several additional terms that must be set.

`ConstantViscosity` is for tests, for benchmark comparisons against analytic results, and for layers whose temperature is unknown.

## Python API

```python
from TidalPy.Viscosity import (
    ArrheniusViscosity,
    ConstantViscosity,
    ReferenceViscosity,
    make_viscosity
)

viscosity_model = ReferenceViscosity(
    reference_viscosity=1.0e22,
    reference_temperature=1000.0,
    molar_activation_energy=3.0e5,
    molar_activation_volume=0.0
)

eta = viscosity_model.calc_viscosity(temperature=1500.0, pressure=1.0e10)   # Pa s

# Name factory: case-insensitive, aliases accepted.
arrhenius_model = make_viscosity(
    "arr",
    {
        "arrhenius_coeff": 1.1e7,
        "grain_size_expo": 2.0,
        "additional_temp_dependence": True
    }
)
```

Constructors take every parameter their model uses as a keyword, by its argument name or its config key, with the defaults from the table above. `make_viscosity(model_name, config=None)` resolves a name or alias case-insensitively and builds that model from `config`, keyed by config key. Absent keys (all of them for `config=None`) take the model's defaults. An unknown name, a key the model does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted name or key.

| Member | Returns | Description |
|---|---|---|
| `calc_viscosity(temperature, pressure=0.0, radius=nan)` | `float` or `np.ndarray` \[Pa s\] | Dynamic viscosity; floats give a float, arrays broadcast together. The radius is read only by the interpolated model. |
| `model_name` | `str` | The canonical model name (`arrhenius`, `reference`, `constant`, `interpolate`, `composite`). |
| `parameters` | `dict` | Every parameter by argument name. |
| `get_parameter(name)` | value | One parameter by argument name or config key. |
| `get_parameter_info()` | `list` of `dict` | Each parameter's name, config key, kind, default, bounds, and description. |
| `with_parameters(**changes)` | model | A new model with some parameters changed; this one is unchanged. |
| `get_config_dict()` | `dict` | `model` plus every parameter under its config key, ready for `make_viscosity`. |
| `save_config(path)` | - | Writes that dict to a TOML file. |

Parameters also read as attributes under their argument names (`model.reference_viscosity`, `model.stress_expo`). Models are not changed in place, so one model can be attached to several layers.

```python
import numpy as np

# A temperature sweep at one pressure
eta = viscosity_model.calc_viscosity(
    np.linspace(1200.0, 1800.0, 50),
    1.0e9
)

# A copy with a dislocation-creep stress exponent; arrhenius_model is unchanged
dislocation_model = arrhenius_model.with_parameters(
    stress_expo=3.5
)
```

### Attaching a Model to a `Layer`

```python
from TidalPy.Material import Material, Phase
from TidalPy.Structures.layers import Layer
from TidalPy.Viscosity import make_viscosity

rock = Phase(
    eos={"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 1.0e11},
    shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10},
    shear_viscosity=make_viscosity(
        "reference",
        {"reference_viscosity_pas": 1.0e20}),
    bulk_viscosity=make_viscosity(
        "constant",
        {"reference_viscosity_pas": 1.0e20}),
    shear_rheology="maxwell")

mantle = Layer(
    "mantle",
    0,
    0.0,
    1.0e6,
    material=Material(solid=rock),
    temperature=1200.0)
print(mantle.calc_state(1.0e9)["shear_viscosity"])   # [Pa s] at 1 GPa and the layer's temperature
```

A viscosity model belongs to a phase of the layer's material, which shares it rather than copying it. The world's equation-of-state solve evaluates the material, and so the model, at the local temperature and pressure as it integrates, and the material's melt weakening then lowers what it gives wherever melt is present. Read the outcome back with `get_shear_viscosity(radius)` on the layer or the world. The declarative form is a `[layers.<name>.material.solid.shear_viscosity]` table in a world's TOML. See [Phases and Materials](../Material/materials.md) and the [TOML schema](../Structures/config/toml_schema.md).

## C++ API

The C++ layer is canonical and the Cython classes are thin adapters over it.

`c_ViscosityBase : c_PhysicsBase` (in `viscosity_base_.hpp`) declares `calc_viscosity(const c_ThermoPoint& point) const noexcept` pure virtual (the point holds the pressure, temperature, and radius) and adds `calc_viscosity_vectorize(temperature, pressure, radius, out_viscosity)`, which broadcasts a single value against a vector. The concrete models `c_ArrheniusViscosity`, `c_ReferenceViscosity`, `c_ConstantViscosity`, `c_InterpolatedViscosity`, and `c_CompositeViscosity` live in `viscosity_.hpp`. Each derives from `c_SpecModel<Model, c_ViscosityBase>` (`Utilities/classes/spec_model_.hpp`) and declares its parameters once, in `parameter_specs()`: argument name, config key, member, default, bounds, and a one-line description. Construction from a `c_ParamMap`, validation, the config entries, the binary record (written by key), copies, and `with_parameters` all follow from that table.

`c_viscosity_registry()` lists each model's names (canonical first, then aliases), binary class id, and constructor. The family's entry points are one line each over the generic registry functions of `Utilities/classes/registry_.hpp`: `c_find_viscosity(name, params)` returns a `std::unique_ptr<c_ViscosityBase>` and throws `std::invalid_argument` for an unknown name or parameter, `c_viscosity_from_binary(stream, force=false)` builds the model a record names and reads it, and `c_viscosity_canonical_name` and `c_viscosity_model_names` resolve names.

## Adding a New Model

1. Add `c_<Name>Viscosity : c_SpecModel<c_<Name>Viscosity, c_ViscosityBase>` to `viscosity_.hpp`: its `parameter_specs()` table, `C_CLASS_ID`, two constructors that call `p_initialize`, and `calc_viscosity`. Override `p_validate` for checks across parameters and `p_update_derived` for cached values.
2. Reserve the next free `BinaryClassID` in the 800 block in `Utilities/binary/binary_.hpp`.
3. Add one row to `c_viscosity_registry()`.
4. Add a two-line Cython subclass (a docstring and `MODEL_NAME`) to `viscosity.pyx`, include it in `_VISCOSITY_CLASSES`, export it from `__init__.py`, add its physics tests to `Tests/Test_Viscosity/`, and add it to this page. The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover its parameters, config, binary record, and errors without changes.

## References

- Moore, W. B. (2006). Thermal equilibrium in Europa's ice shell. *Icarus*, 180(1), 141-146. Arrhenius flow law with activation energy and volume.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015. Reference-viscosity (relative activation) law.
