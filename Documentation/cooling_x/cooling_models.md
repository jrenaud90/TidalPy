# Cooling Models (`cooling_x`)

_Updated: 2026-09-12_

A cooling model maps a layer's physical state onto a **cooling result**: the surface heat flux $q$ [W m$^{-2}$], the thermal boundary-layer thickness [m], and the Rayleigh and Nusselt numbers. The heat flux is what drives a layer's thermal evolution, and the boundary-layer thickness is what makes convective transport so much more effective than conduction: the same temperature drop is squeezed across a thin layer at the top instead of the whole interior.

Each `SolidLiquidLayer` can hold one cooling model. The math mirrors the validated classic implementation in `TidalPy/cooling/cooling_models.py`.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_CoolingBase   (abstract)
              ├── c_OffCooling         aliases "off", "none"
              ├── c_ConvectiveCooling  aliases "convection", "convective"
              └── c_ConductiveCooling  aliases "conduction", "conductive"
```

`c_CoolingBase` declares `calc_cooling(inputs)` pure virtual and supplies the vectorized loops, the configuration export, and the binary encoding.

## Inputs and result

A cooling evaluation depends on eight physical quantities. In C++ they are bundled into a `c_CoolingInputs` struct, following the style guide's rule against long argument lists; the Python `calc_cooling` takes them as explicit arguments in the same order.

| Input | Units | Meaning |
|---|---|---|
| `delta_temp` | K | Temperature drop across the layer. |
| `thickness` | m | Layer or sub-layer thickness. |
| `gravity` | m s$^{-2}$ | Gravitational acceleration. |
| `density` | kg m$^{-3}$ | Bulk density. |
| `viscosity` | Pa s | Dynamic viscosity. |
| `thermal_conductivity` | W m$^{-1}$ K$^{-1}$ | Thermal conductivity. |
| `thermal_diffusivity` | m$^2$ s$^{-1}$ | Thermal diffusivity. |
| `thermal_expansion` | K$^{-1}$ | Thermal expansivity. |

The result is a `CoolingResult` carrying `cooling_flux` [W m$^{-2}$], `boundary_layer_thickness` [m], `rayleigh`, and `nusselt`. Each field is a Python `float` for a scalar evaluation and a `float64` array for a vectorized one. The result supports `to_dict()` and unpacks in field order.

```python
from TidalPy.cooling_x import ConvectiveCooling

model = ConvectiveCooling()
result = model.calc_cooling(1000.0, 1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)
flux, boundary_layer, rayleigh, nusselt = result
```

## The three models

| Model | Cooling flux $q$ [W m$^{-2}$] | Boundary layer | Ra | Nu |
|---|---|---|---|---|
| `OffCooling` | $0$ | $0.5 \times$ thickness | 0 | 1 |
| `ConductiveCooling` | $k \, \Delta T / d$ | thickness $d$ | 0 | 1 |
| `ConvectiveCooling` | $k \, \Delta T / \delta$ | $\delta = d / \mathrm{Nu}$ | see below | see below |

`OffCooling` holds the heat in. It is the way to run a layer adiabatically, or to isolate the effect of a heating term by removing the sink. `ConductiveCooling` transports heat by diffusion alone, which is the right description for a stagnant lid, a thin ice shell, or any layer whose Rayleigh number sits below critical. `ConvectiveCooling` is the parameterized boundary-layer model and the usual choice for a mantle.

For convection,

$$\mathrm{Ra} = \frac{\alpha_\mathrm{th} \rho g \, \Delta T \, d^3}{\eta \kappa}, \qquad \mathrm{Nu} = \max\left[ a \left( \frac{\mathrm{Ra}}{\mathrm{Ra}_\mathrm{crit}} \right)^{b},\; 2 \right], \qquad \delta = \frac{d}{\mathrm{Nu}}, \qquad q = \frac{k \, \Delta T}{\delta}$$

with $\alpha_\mathrm{th}$ the thermal expansivity, $\kappa$ the thermal diffusivity, and the fitted constants $a$ = `convection_alpha` (default 1.0), $b$ = `convection_beta` (default 1/3), and $\mathrm{Ra}_\mathrm{crit}$ = `critical_rayleigh` (default 1100.0).

The Rayleigh number is the ratio of the buoyancy driving a hot parcel upward to the diffusion bleeding its heat away; above the critical value convection sets in. The Nusselt number is how many times more heat that convection carries than conduction alone would, and the boundary layer is thinned by exactly that factor. The exponent of 1/3 is the classical boundary-layer result, which has the useful consequence that the convective flux is independent of the layer thickness.

The Nusselt number is floored at 2, which is the stagnant-lid limit: even a barely convecting layer moves about twice the conductive flux. Degenerate inputs take a defined path rather than producing a division by zero. A non-positive temperature drop, or a thickness below the shared `minimum_layer_thickness` configuration floor, sets the Rayleigh number to zero and the Nusselt number to 2, matching the classic implementation's edge behavior.

## Building and evaluating a model

```python
from TidalPy.cooling_x import ConvectiveCooling, make_cooling

convective_cooling = ConvectiveCooling(convection_alpha=1.0, convection_beta=1.0/3.0,
                                       critical_rayleigh=1100.0)

# delta_temp, thickness, gravity, density, viscosity, conductivity, diffusivity, expansion
result = convective_cooling.calc_cooling(1000.0, 1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)

# Name factory: case-insensitive, aliases accepted.
conductive_cooling = make_cooling("conductive")
```

`make_cooling(model_name, config=None)` resolves a name or alias case-insensitively and reads `convection_alpha`, `convection_beta`, and `critical_rayleigh` from `config`. Unknown names raise `ValueError`. The convection parameters are read-only properties on `ConvectiveCooling`; the other two models carry none. The resolved `model_name` is `off`, `conduction`, or `convection`, and that is the name written to a configuration dict.

### Factory internals

At the C++ level the factory is enum-based, matching the other physics modules. `c_CoolingModel` names one value per model, `c_cooling_model_from_name(name)` maps a case-insensitive name or alias onto it and throws `std::invalid_argument` for an unknown name, and `c_find_cooling(model, config)` returns a `std::unique_ptr<c_CoolingBase>`. The Python `make_cooling` wraps both and adopts the returned pointer into the matching wrapper.

## Vectorized evaluation

Two of the eight inputs change as a layer evolves thermally: the temperature drop and the viscosity. The vectorized methods, defined once on the base class so every model inherits them, sweep exactly those and hold the rest fixed.

| Method | Sweeps |
|---|---|
| `calc_cooling_vectorize_temperature(delta_temp[], ...)` | A temperature-drop array at fixed viscosity. |
| `calc_cooling_vectorize_viscosity(..., viscosity[], ...)` | A viscosity array at fixed temperature drop. |
| `calc_cooling_vectorize_all(delta_temp[], ..., viscosity[], ...)` | Two equal-length arrays, element by element. |

```python
import numpy as np
from TidalPy.cooling_x import ConvectiveCooling

model = ConvectiveCooling()
sweep = model.calc_cooling_vectorize_temperature(np.linspace(100.0, 2000.0, 50),
                                                 1.0e6, 9.8, 3300.0, 1.0e21,
                                                 4.0, 1.0e-6, 3.0e-5)
flux_profile = sweep.cooling_flux   # float64 array
```

At the C++ level each fills a caller-supplied `std::vector<c_CoolingResult>`, copying the base inputs and overriding the swept field per element. Mismatched lengths raise `ValueError` from Python. The Cython wrappers return one `CoolingResult` whose four fields are arrays.

## Convenience functions

For a one-shot evaluation with no model object left behind:

```python
import numpy as np
from TidalPy.cooling_x import conductive, convective, cooling_off

result = convective(1000.0, 1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)

# delta_temp and viscosity may be floats or arrays and are broadcast together;
# the remaining inputs are scalar constants.
sweep = convective(np.linspace(100.0, 2000.0, 50), 1.0e6, 9.8, 3300.0,
                   1.0e21, 4.0, 1.0e-6, 3.0e-5)
```

The signatures are `cooling_off(delta_temp, thickness)`, `conductive(delta_temp, thickness, thermal_conductivity)`, and `convective(delta_temp, thickness, gravity, density, viscosity, thermal_conductivity, thermal_diffusivity, thermal_expansion, convection_alpha=1.0, convection_beta=1/3, critical_rayleigh=1100.0)`. Each builds a stack-allocated C++ model, picks the most specific vectorized routine for the input pattern, and returns a `CoolingResult`. The off and conduction functions take only the inputs they actually use, which is why their argument lists are shorter than `calc_cooling`'s fixed eight.

## Serialization

| Call | Result |
|---|---|
| `get_config_dict()` | `model` plus any model parameters. |
| `save_config(path)` | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path, force=False)` | TidalPy binary format, class ids 401 through 403. Off and conduction write no parameters; convection writes its three scalars. |

A cooling model attached to a layer is written as part of that layer's binary record and rebuilt recursively on load. See [Binary serialization](../utilities_x/binary_x.md).

## Adding a new cooling model

To add a cooling model named `Foo`:

**C++ (`TidalPy/cooling_x/cooling_.hpp`)**

1. If `Foo` needs new parameters, add them to `c_CoolingConfig` with defaults.
2. Add a `cool_foo(const c_CoolingInputs&[, const c_CoolingConfig&])` free function implementing the heat-transport law. Use the `cool_guard` floor on any denominator.
3. Add the class `c_Foo : public c_CoolingBase` with constructors `c_Foo()` and `explicit c_Foo(const c_CoolingConfig&)` that pass a model-name string to the base and copy any parameters into `p_*` members, a `get_*` accessor per parameter, an override of `calc_cooling(const c_CoolingInputs&)` returning a `c_CoolingResult`, and overrides of `write_binary` / `read_binary` built on the base helpers.

**C++ (`TidalPy/Utilities_x/binary_x/binary_.hpp`)**

4. Add a unique `BinaryClassID::Foo` value. Cooling models occupy the 40X range.

**C++ factory (`cooling_.hpp`)**

5. Add `Foo` to the `c_CoolingModel` enum, map its name and aliases in `c_cooling_model_from_name`, and add a `case` to `c_find_cooling`. Add a `case BinaryClassID::Foo` to `c_cooling_from_binary` as well, so a `Foo` attached to a layer can be rebuilt when the layer is loaded.

**Cython (`cooling.pxd` / `cooling.pyx`)**

6. Declare `c_Foo` (constructors and getters) in `cooling.pxd` and add the enum value to the `c_CoolingModel` cimport.
7. Add the `cdef class Foo(CoolingBase)` wrapper in `cooling.pyx` with its parameter properties, the adoption branch in `make_cooling`, and a lower-case `foo(...)` convenience function. The config dict comes from the C++ `append_config_entries` override.

**Package, tests, and documentation**

8. Export `Foo` and `foo` from `TidalPy/cooling_x/__init__.py`, and the C++ names from `__init__.pxd`.
9. Add `Foo` to the lists in `Tests/Test_Cooling_x/test_cooling_01.py`, which cover the model name, the result against an independent reference, the factory and aliases, vectorization, the config dict, the binary round trip, and `isinstance`.
10. Document the model here with its formula, parameters, and references.

No build-system change is needed; `cooling_x.cooling` is already registered in `cython_extensions.json`.

## References

- Turcotte, D. L., and Schubert, G. (2002). *Geodynamics*, second edition. Rayleigh and Nusselt convection scaling, and conduction.
- Solomatov, V. S. (1995). Scaling of temperature- and stress-dependent viscosity convection. *Physics of Fluids*, 7(2), 266-274. Stagnant-lid scaling.
- Schubert, G., Turcotte, D. L., and Olson, P. (2001). *Mantle Convection in the Earth and Planets*. Boundary-layer theory.
