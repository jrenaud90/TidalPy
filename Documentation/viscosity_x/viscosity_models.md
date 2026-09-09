# Viscosity Models (`viscosity_x`)

A **viscosity** model returns a material's dynamic viscosity [Pa·s] as a function
of temperature [K] and pressure [Pa]. This is the pre-melt ("solid") viscosity
that the [partial-melt](../partial_melt_x/partial_melt_models.md) step weakens;
both are frequency-independent and cached once per EOS solve in the whole-planet
EOS pipeline.

All quantities are **MKS**. The math mirrors the validated legacy implementation
in `TidalPy/rheology/viscosity/viscosity_models.py`. The molar gas constant `R`
comes from the shared TidalPy config (`tidalpy_config_ptr->d_R`).

---

## Models

| Model | Alias | Formula |
|-------|-------|---------|
| `ArrheniusViscosity` | `arr` | `η = A·σ^(1−n)·d^m·exp((E_a + P·V_a)/(R·T))` (× `T` if `additional_temp_dependence`). |
| `ReferenceViscosity` | `ref` | `η = η_ref·exp(((E_a + P·V_a)/R)·(1/T − 1/T_ref))`. |
| `ConstantViscosity` | `const` | `η = reference_viscosity` (independent of T, P). |

where `E_a` is the molar activation energy [J/mol], `V_a` the molar activation
volume [m³/mol], `σ` the applied stress [Pa] with exponent `n`, and `d` the grain
size [m] with exponent `m`.

**References:** Moore (2006) — Arrhenius flow law; Henning (2009) —
reference-viscosity law.

---

## Python API

```python
from TidalPy.viscosity_x import (
    ArrheniusViscosity, ReferenceViscosity, ConstantViscosity, make_viscosity)

visc = ReferenceViscosity(
    reference_viscosity=1.0e22, reference_temperature=1000.0,
    molar_activation_energy=3.0e5, molar_activation_volume=0.0)

eta = visc.calc_viscosity(temperature_k=1500.0, pressure_pa=1.0e10)   # [Pa·s]

# Name factory (aliases, case-insensitive):
visc = make_viscosity("arr", {"arrhenius_coeff": 1.1e7, "grain_size_expo": 2.0,
                              "additional_temp_dependence": True})
```

**Constructors / parameters**

- `ConstantViscosity(reference_viscosity=1e22)`
- `ReferenceViscosity(reference_viscosity=1e22, reference_temperature=1000,
  molar_activation_energy=3e5, molar_activation_volume=0)`
- `ArrheniusViscosity(arrhenius_coeff=1.0, stress=1.0, stress_expo=1.0,
  grain_size=1e-3, grain_size_expo=0.0, molar_activation_energy=3e5,
  molar_activation_volume=0.0, additional_temp_dependence=False)`

**Methods / properties**

| Member | Returns | Description |
|--------|---------|-------------|
| `calc_viscosity(T, P=0.0)` | float [Pa·s] | Dynamic viscosity. |
| `get_config_dict()` | dict | Model name + parameters. |
| `save_config(path)` / `save_binary(path)` / `load_binary(path)` | — | Inherited from `TidalPyBaseClass`. |

`make_viscosity(name, config=None)` returns the concrete subclass; absent config
keys fall back to the C++ defaults. Unknown names raise `ValueError`.

---

## C++ API

The C++ layer is canonical; the Cython classes are thin adapters.

**Config struct** `tidalpy::c_ViscosityConfig` (defaults): `reference_viscosity`
(1e22), `reference_temperature` (1000), `molar_activation_energy` (3e5),
`molar_activation_volume` (0), `arrhenius_coeff` (1.0), `stress` (1.0),
`stress_expo` (1.0), `grain_size` (1e-3), `grain_size_expo` (0.0),
`additional_temp_dependence` (false).

**Base class** `c_ViscosityBase : c_PhysicsBase` (in `viscosity_base_.hpp`):
pure virtual `calc_viscosity(double temperature_k, double pressure_pa) const`,
plus `calc_viscosity_vectorize(temperature, pressure, out_viscosity)` (the primary
radial sweep — one entry per slice).

**Concrete models** (`viscosity_.hpp`): `c_ArrheniusViscosity`,
`c_ReferenceViscosity`, `c_ConstantViscosity` (each with parameter getters). Binary
class ids **801–803**.

**Enum factory:**
`c_viscosity_model_from_name(const std::string&) -> c_ViscosityModel`
(`Arrhenius`/`Reference`/`Constant`, throws `std::invalid_argument` on unknown
names); `c_find_viscosity(c_ViscosityModel, const c_ViscosityConfig&) ->
std::unique_ptr<c_ViscosityBase>` (and a name overload);
`c_viscosity_from_binary(std::istream&, bool force=false)` (peek class id → build →
`read_binary`).

---

## Adding a new model

1. Add a `c_<Name>Viscosity : c_ViscosityBase` in `viscosity_.hpp` implementing
   `calc_viscosity`, `write_binary`, `read_binary`.
2. Add its parameters to `c_ViscosityConfig`.
3. Register a `c_ViscosityModel::<Name>` enum value and wire it into
   `c_viscosity_model_from_name`, `c_find_viscosity`, and `c_viscosity_from_binary`.
4. Reserve a unique `BinaryClassID` (next free in the 800-block) in
   `Utilities_x/binary_x/binary_.hpp`.
5. Add the Cython `cdef class`, factory branch, config-dict, and tests + this doc.
