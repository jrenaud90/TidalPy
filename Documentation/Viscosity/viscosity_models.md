# Viscosity Models (`Viscosity`)

_Updated: 2026-10-07_

A viscosity model returns a phase's dynamic viscosity $\eta$ \[Pa s\] at a temperature \[K\] and pressure \[Pa\]. This is the pre-melt viscosity, which the material's [melt weakening](../PartialMelt/partial_melt_models.md#melt-weakening) then lowers. Both are frequency-independent, so they are resolved once per equation-of-state solve and reused at every tidal forcing frequency.

Solid-state creep is thermally activated: viscosity falls exponentially with temperature through the Boltzmann factor $\exp(E_a / RT)$ and rises with pressure through an activation volume. The models differ in how that exponential is anchored: to an absolute flow-law prefactor, to a reference viscosity at a reference temperature, or not at all.

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

## Choosing a Model

- `ReferenceViscosity` is the usual choice for a planetary interior. It is anchored to a viscosity measured or inferred at a stated temperature, and its one strong assumption is the activation energy: doubling $E_a$ steepens the response to temperature without moving the anchor point.
- `ArrheniusViscosity` reproduces a published flow law, usually quoted as an absolute prefactor with stress and grain-size exponents. It gives full control at the cost of more terms to set.
- `ConstantViscosity` is for tests, benchmarks against analytic results, and layers whose temperature is unknown.
- `InterpolatedViscosity` carries a viscosity that varies with radius inside a layer, such as a seismic profile's.
- `CompositeViscosity` combines mechanisms. Published flow laws for ice and olivine combine several (diffusion creep, dislocation creep, grain-boundary sliding), and a law whose activation energy switches at a temperature is two mechanisms.

## Models

| Model | Aliases | Viscosity $\eta$ \[Pa s\] |
|---|---|---|
| `ArrheniusViscosity` | `arrhenius`, `arr` | $A \, \sigma^{1-n} d^{\,m} \exp\left(\dfrac{E_a + P V_a}{R T}\right)$, multiplied by $T$ when `additional_temp_dependence = True` |
| `ReferenceViscosity` | `reference`, `ref` | $\eta_\mathrm{ref} \exp\left[\dfrac{E_a + P V_a}{R T} - \dfrac{E_a + P_\mathrm{ref} V_a}{R T_\mathrm{ref}}\right]$ |
| `ConstantViscosity` | `constant`, `const` | $\eta_\mathrm{ref}$, independent of temperature and pressure |
| `InterpolatedViscosity` | `interpolate`, `interp`, `interpolated` | Linear in radius between the points of a `radius_m` and `viscosity_pas` table, held at the end values beyond it |
| `CompositeViscosity` | `composite`, `parallel` | $\left(\sum_i 1 / \eta_i\right)^{-1}$ over its `mechanisms`, each a viscosity model: deformation mechanisms acting in parallel at a common stress, so the weakest dominates |

Each parameter has a constructor keyword (also an attribute) and a config key (TOML tables, `make_viscosity` configs, `get_config_dict()`); a dimensional config key ends in its unit.

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

The stress exponent $n$ distinguishes creep regimes. With $n = 1$ the material is in diffusion creep, where the flow law is linear and the stress term drops out; with $n > 1$ it is in dislocation creep, where the viscosity decreases as the stress increases. The grain-size exponent plays the same role for grain-boundary processes. `additional_temp_dependence = True` adds the explicit factor of $T$ that some published diffusion-creep laws carry in front of the exponential.

$\eta_\mathrm{ref}$ is the viscosity at $T_\mathrm{ref}$ and $P_\mathrm{ref}$. With the default $P_\mathrm{ref} = 0$ a positive activation volume always raises the viscosity with pressure, as in the Arrhenius law; a deep layer with a large activation volume is better anchored at a pressure inside the layer.

`CompositeViscosity(diffusion_law, dislocation_law)`, or a TOML table `model = "composite"` with `mechanisms = [{ model = ... }, ...]`, combines mechanisms. The interpolated model reads only the radius, so its layer must keep its volume. A seismic profile's quality factor travels in its table for the `seismic_q` rheology, which reads its viscosity input as $Q$.

## Behavior at the Limits

- The Arrhenius and reference models return infinity at or below zero temperature, and at very cold positive temperatures where the exponential overflows. This is intended: the cold limit of a thermally activated fluid is a solid that does not flow, and an infinite viscosity makes the rheology models return a purely elastic response.
- `ConstantViscosity` returns its constant at any temperature, and may be infinite (a rigid, purely elastic material). `InterpolatedViscosity` returns its table value.
- A non-positive reference temperature is refused when the model is built.

## Python API

Constructors take every parameter their model uses as a keyword, by argument name or config key, with the defaults above. `make_viscosity(model_name, config=None)` resolves a name or alias case-insensitively and builds that model from `config`, keyed by config key; absent keys (all of them for `config=None`) take the defaults. An unknown name, a key the model does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted name or key.

| Member | Returns | Description |
|---|---|---|
| `calc_viscosity(temperature, pressure=0.0, radius=nan)` | `float` or `np.ndarray` \[Pa s\] | Dynamic viscosity; floats give a float, arrays broadcast together. The radius is read only by the interpolated model. |
| `model_name` | `str` | The canonical model name (`arrhenius`, `reference`, `constant`, `interpolate`, `composite`). |
| `save_config(path)` | - | Writes `get_config_dict()` to a TOML file. |

`parameters`, `get_parameter`, `get_parameter_info`, `with_parameters`, and `get_config_dict` are shared by every physics model (see [`PhysicsBase`](../Utilities/classes.md#physicsbase)). Parameters read as attributes under their argument names (`model.reference_viscosity`). Models are not changed in place, so one model can serve several layers.

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

A viscosity model fills a phase's `shear_viscosity` or `bulk_viscosity` slot (`Phase(shear_viscosity=make_viscosity("reference", {...}))`, or a `[layers.<name>.material.solid.shear_viscosity]` table in a world's TOML), and the phase shares it rather than copying it (see [Building a Material](../Material/materials.md#building-a-material)). The world's EOS solve evaluates it at the local temperature and pressure, and melt weakening then lowers it wherever melt is present. Read the result with `get_shear_viscosity(radius)` on the layer or the world, or at one point with `layer.calc_state(pressure)["shear_viscosity"]` (at the layer's temperature).

## C++ API

`c_ViscosityBase` (`viscosity_base_.hpp`) declares `calc_viscosity(const c_ThermoPoint& point) const noexcept` (the point holds the pressure, temperature, and radius) and provides `calc_viscosity_vectorize(temperature, pressure, radius, out_viscosity)`, which broadcasts a single value against a vector. The models `c_ArrheniusViscosity`, `c_ReferenceViscosity`, `c_ConstantViscosity`, `c_InterpolatedViscosity`, and `c_CompositeViscosity` are in `viscosity_.hpp`, and the Cython classes wrap them. `c_find_viscosity(name, params)` returns a `std::unique_ptr<c_ViscosityBase>` and throws `std::invalid_argument` for an unknown name or parameter; `c_viscosity_from_binary(stream, force=false)` rebuilds a saved model; `c_viscosity_canonical_name` and `c_viscosity_model_names` resolve names.

## Adding a New Model

1. Add `c_<Name>Viscosity : c_SpecModel<c_<Name>Viscosity, c_ViscosityBase>` to `viscosity_.hpp` with its `parameter_specs()` table (argument name, config key, member, default, bounds, description), `C_CLASS_ID`, two constructors that call `p_initialize`, and `calc_viscosity`. Override `p_validate` for cross-parameter checks and `p_update_derived` for cached values.
2. Reserve the next free `BinaryClassID` in the 800 block of `Utilities/binary/binary_.hpp`, and add a row to `c_viscosity_registry()`.
3. Add a Cython subclass (a docstring and `MODEL_NAME`) to `viscosity.pyx`, include it in `_VISCOSITY_CLASSES`, and export it from `__init__.py`.
4. Add physics tests to `Tests/Test_Viscosity/` and document the model here. The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover parameters, config, binary record, and errors.

## References

- Moore, W. B. (2006). Thermal equilibrium in Europa's ice shell. *Icarus*, 180(1), 141-146. Arrhenius flow law with activation energy and volume.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015. Reference-viscosity (relative activation) law.
