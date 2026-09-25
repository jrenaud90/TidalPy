# Viscosity Models (`viscosity_x`)

_Updated: 2026-09-23_

A viscosity model returns a material's dynamic viscosity $\eta$ \[Pa s\] as a function of temperature \[K\] and pressure \[Pa\]. This is the pre-melt, or "solid", viscosity: the value a material would show with no melt present, which the [partial-melt](../partial_melt_x/partial_melt_models.md) step then weakens. Both are frequency-independent, so both are resolved once per equation-of-state solve and reused across every tidal forcing frequency.

Solid-state creep is thermally activated: viscosity falls exponentially with temperature through the Boltzmann factor $\exp(E_a / RT)$ and rises with pressure through an activation volume. The models differ in how that exponential is anchored: to an absolute flow-law prefactor, to a reference viscosity at a reference temperature, or not at all.

## Models

| Model | Aliases | Viscosity $\eta$ [Pa s] |
|---|---|---|
| `ArrheniusViscosity` | `arrhenius`, `arr` | $A \, \sigma^{1-n} d^{\,m} \exp\left(\dfrac{E_a + P V_a}{R T}\right)$, multiplied by $T$ when `additional_temp_dependence` is set |
| `ReferenceViscosity` | `reference`, `ref` | $\eta_\mathrm{ref} \exp\left[\dfrac{E_a}{R} \left(\dfrac{1}{T} - \dfrac{1}{T_\mathrm{ref}}\right) + \dfrac{P V_a}{R T}\right]$ |
| `ConstantViscosity` | `constant`, `const` | $\eta_\mathrm{ref}$, independent of temperature and pressure |

Each parameter carries two names: the constructor keyword, which is also the read-only property, and the config key used in a TOML table, a `make_viscosity` config dictionary, and `get_config_dict()`. A dimensional config key ends in its unit; the code name does not.

| Parameter | Config key | Symbol | Default | Units | Used by |
|---|---|---|---|---|---|
| `reference_viscosity` | `reference_viscosity_pas` | $\eta_\mathrm{ref}$ | 1.0e22 | Pa s | Reference, Constant |
| `reference_temperature` | `reference_temperature_k` | $T_\mathrm{ref}$ | 1000.0 | K | Reference |
| `molar_activation_energy` | `molar_activation_energy_j_mol` | $E_a$ | 3.0e5 | J mol$^{-1}$ | Arrhenius, Reference |
| `molar_activation_volume` | `molar_activation_volume_m3_mol` | $V_a$ | 0.0 | m$^3$ mol$^{-1}$ | Arrhenius, Reference |
| `arrhenius_coeff` | `arrhenius_coeff` | $A$ | 1.0 | model-dependent | Arrhenius |
| `stress` | `stress_pa` | $\sigma$ | 1.0 | Pa | Arrhenius |
| `stress_expo` | `stress_expo` | $n$ | 1.0 | - | Arrhenius |
| `grain_size` | `grain_size_m` | $d$ | 1.0e-3 | m | Arrhenius |
| `grain_size_expo` | `grain_size_expo` | $m$ | 0.0 | - | Arrhenius |
| `additional_temp_dependence` | `additional_temp_dependence` | - | `False` | - | Arrhenius |

The stress exponent $n$ distinguishes creep regimes: $n = 1$ is diffusion creep, where the flow law is linear and the stress term drops out, and $n > 1$ is dislocation creep, where the material shears more readily the harder it is pushed. The grain-size exponent plays the same role for grain-boundary processes. Setting `additional_temp_dependence` adds the explicit factor of $T$ that some published diffusion-creep flow laws carry in front of the exponential.

### Behavior at the Limits

The reference viscosity $\eta_\mathrm{ref}$ is the viscosity at $T_\mathrm{ref}$ and zero pressure, so a positive activation volume always raises the viscosity with pressure, as it does in the Arrhenius law.

The Arrhenius and reference models return infinity at or below zero temperature (`ConstantViscosity` returns its constant at any temperature): the cold limit of a thermally activated fluid is a solid that does not flow, and an infinite viscosity makes the rheology models return a purely elastic response. A reference model with a non-positive reference temperature also returns infinity. Very cold but positive temperatures reach infinity by overflowing the exponential, which is deliberate rather than guarded against.

### Choosing a Model

`ReferenceViscosity` is the usual choice for a planetary interior. It is anchored to a viscosity measured or inferred at a stated temperature, and its one strong assumption is the activation energy. Doubling $E_a$ steepens the response to temperature without moving the anchor point.

`ArrheniusViscosity` is the choice when reproducing a published flow law, which is usually quoted as an absolute prefactor with stress and grain-size exponents. It gives full control at the cost of several additional terms that must be set.

`ConstantViscosity` is for tests, for benchmark comparisons against analytic results, and for layers whose temperature is unknown.

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

Constructors take every parameter their model uses as a keyword with the default from the table above:

`ConstantViscosity(reference_viscosity=1.0e22)`

`ReferenceViscosity(reference_viscosity=1.0e22, reference_temperature=1000.0, molar_activation_energy=3.0e5, molar_activation_volume=0.0)`

`ArrheniusViscosity(arrhenius_coeff=1.0, stress=1.0, stress_expo=1.0, grain_size=1.0e-3, grain_size_expo=0.0, molar_activation_energy=3.0e5, molar_activation_volume=0.0, additional_temp_dependence=False)`

`make_viscosity(model_name, config=None)` resolves a name or alias case-insensitively and reads the parameters it recognizes from `config`. Absent keys fall back to the model default and keys another viscosity model uses are ignored. An unrecognized name raises `ValueError`, and so does a key that no viscosity model reads, with the closest accepted key named in the message.

| Member | Returns | Description |
|---|---|---|
| `calc_viscosity(temperature, pressure=0.0)` | `float` [Pa s] | Dynamic viscosity at those conditions. |
| `model_name` | `str` | The resolved model name (`arrhenius`, `reference`, `constant`). |
| `get_config_dict()` | `dict` | `model` plus every parameter the model carries. |
| `save_config(path)` | - | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path, force=False)` | - | TidalPy binary format; see [Binary serialization](../utilities_x/binary_x.md). |

Parameters are read-only properties under their code names: `reference_viscosity` on the constant model; `reference_viscosity`, `reference_temperature`, `molar_activation_energy`, and `molar_activation_volume` on the reference model; and every constructor keyword on the Arrhenius model: `arrhenius_coeff`, `stress`, `stress_expo`, `grain_size`, `grain_size_expo`, `molar_activation_energy`, `molar_activation_volume`, and `additional_temp_dependence`. `get_config_dict()` emits the config keys, so a dictionary read back from a model or a TOML file feeds straight into `make_viscosity`.

### Attaching a Model to a `Layer`

```python
from TidalPy.Material_x.eos import ConstantDensityEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.structures_x.layers.physics import PhysicsLayer

mantle = PhysicsLayer("mantle", 0, 0.0, 1.0e6, 2.1e19)
mantle.set_eos(ConstantDensityEOS(shear_modulus_static=50.0e9, bulk_modulus_static=100.0e9))

mantle.set_shear_viscosity(make_viscosity("reference", {"reference_viscosity_pas": 1.0e20}))
mantle.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e20}))
```

A viscosity model belongs to the layer's material, which is defined via its EOS model. The layer's `set_shear_viscosity` and `set_bulk_viscosity` are helpers that hand the model to the attached EOS (so attach the EOS first), and the same two methods are on the EOS model itself. Ownership of the C++ model transfers, as it does for a rheology. The world's equation-of-state solve evaluates the model at the local temperature and pressure as it integrates, a table the EOS model carries overrides the result, and the partial-melt model then weakens what remains. Read the outcome back with `get_shear_viscosity(radius)` on the layer or the world. The declarative form is a `[layers.<name>.material.shear_viscosity]` table in a world's TOML; see the [TOML schema](../structures_x/config/toml_schema.md).

## C++ API

The C++ layer is canonical and the Cython classes are thin adapters over it.

`c_ViscosityConfig` holds every parameter for every model in one struct, with the defaults listed above. `c_ViscosityBase : c_PhysicsBase` (in `viscosity_base_.hpp`) declares `calc_viscosity(double temperature, double pressure) const` pure virtual and adds `calc_viscosity_vectorize(temperature, pressure, out_viscosity)`, the radial sweep the equation-of-state solve uses. The concrete models `c_ArrheniusViscosity`, `c_ReferenceViscosity`, and `c_ConstantViscosity` live in `viscosity_.hpp`, each with a getter per parameter, and occupy binary class ids 801 through 803.

The factory mirrors the rheology one: `c_viscosity_model_from_name(name)` maps a name or alias onto the `c_ViscosityModel` enum and throws `std::invalid_argument` for an unknown name, `c_find_viscosity(model, config)` returns a `std::unique_ptr<c_ViscosityBase>` (a name overload does both steps), and `c_viscosity_from_binary(stream, force=false)` peeks the class id, builds the matching model, and calls `read_binary`.

## Adding a New Model

1. Add the model's parameters to `c_ViscosityConfig` in `viscosity_.hpp`, with defaults.
2. Add `c_<Name>Viscosity : c_ViscosityBase` implementing `calc_viscosity`, `append_config_entries`, `write_binary`, and `read_binary`.
3. Reserve the next free `BinaryClassID` in the 800 block in `Utilities_x/binary_x/binary_.hpp`.
4. Register a `c_ViscosityModel::<Name>` enum value and wire it into `c_viscosity_model_from_name`, `c_find_viscosity`, and `c_viscosity_from_binary`.
5. Add the Cython `cdef class` with its parameter properties, the adoption branch in `make_viscosity`, tests in `Tests/Test_Viscosity_x/`, and an entry on this page.

## References

- Moore, W. B. (2006). Thermal equilibrium in Europa's ice shell. *Icarus*, 180(1), 141-146. Arrhenius flow law with activation energy and volume.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015. Reference-viscosity (relative activation) law.
