# Cooling Models (`Cooling`)

_Updated: 2026-10-02_

A cooling model says how heat moves through a layer. Its flux law maps the layer's physical state onto a **cooling result**: the surface heat flux $q$ [W m$^{-2}$], the thermal boundary-layer thickness [m], and the Rayleigh and Nusselt numbers. In a world's thermal solve the same model builds its layer's temperature profile from that law. The boundary-layer thickness is what makes convective transport so much more effective than conduction: the same temperature drop is squeezed across a thin layer at the top instead of the whole interior.

Each `Layer` can hold one cooling model.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_CoolingBase   (abstract)
              ├── c_OffCooling         aliases "off", "none"
              ├── c_ConvectiveCooling  aliases "convection", "convective"
              └── c_ConductiveCooling  aliases "conduction", "conductive"
```

`c_CoolingBase` declares `get_temperature_kind()`, `calc_cooling(inputs)`, `build_profile(context, probe, out)`, and `calc_cooling_vectorize(...)` pure virtual. Every concrete model derives from `c_SpecModel`, which gives it its parameters, config dict, and binary record from one table. The canonical names are `off`, `conduction`, and `convection`, and those are the names written to a config dict.

## Inputs and Result

A flux-law evaluation depends on eight physical quantities and one flag. In C++ they are bundled into a `c_CoolingInputs` struct; the Python `calc_cooling` takes them as explicit arguments in the same order.

| Input | Units | Meaning |
|---|---|---|
| `delta_temp` | K | Temperature drop across the layer, the sum of the drops across its boundary layers. |
| `thickness` | m | Layer or sub-layer thickness. |
| `gravity` | m s$^{-2}$ | Gravitational acceleration. |
| `density` | kg m$^{-3}$ | Bulk density. |
| `viscosity` | Pa s | Dynamic viscosity. |
| `thermal_conductivity` | W m$^{-1}$ K$^{-1}$ | Thermal conductivity. |
| `thermal_diffusivity` | m$^2$ s$^{-1}$ | Thermal diffusivity. |
| `thermal_expansion` | K$^{-1}$ | Thermal expansivity. |
| `liquid` | - | The convecting interior is liquid (a magma ocean), so the convection model takes its liquid scaling. Default `False`. |

The result is a `CoolingResult` carrying `cooling_flux` [W m$^{-2}$], `boundary_layer_thickness` [m], `rayleigh`, and `nusselt`. Each field is a Python `float` for a scalar evaluation and a `float64` array for a vectorized one. The result supports `to_dict()` and unpacks in field order.

```python
from TidalPy.Cooling import ConvectiveCooling

model = ConvectiveCooling()
# delta_temp, thickness, gravity, density, viscosity, conductivity, diffusivity, expansion
result = model.calc_cooling(1000.0, 1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)
flux, boundary_layer, rayleigh, nusselt = result
```

## Current Models

| Model | Cooling flux $q$ [W m$^{-2}$] | Boundary layer | Ra | Nu | Profile in a thermal solve |
|---|---|---|---|---|---|
| `OffCooling` | $0$ | $0.5 \times$ thickness | 0 | 1 | Isothermal |
| `ConductiveCooling` | $k \, \Delta T / d$ | thickness $d$ | 0 | 1 | Two conducting halves |
| `ConvectiveCooling` | $k \, \Delta T / \delta$ | $\delta = d / \mathrm{Nu}$ | see below | see below | Boundary layers around an adiabatic interior |

`OffCooling` turns off cooling (no heat transfer inside the layer). `ConductiveCooling` transports heat by diffusion, applicable for a stagnant lid, a thin ice shell, or any layer whose Rayleigh number sits below critical. `ConvectiveCooling` is the parameterized boundary-layer model and the usual choice for a mantle.

### Convection

$$\mathrm{Ra} = \frac{\alpha_\mathrm{th} \rho g \, \Delta T \, d^3}{\eta \kappa}, \qquad \mathrm{Nu} = \max\left[ a \left( \frac{\mathrm{Ra}}{\mathrm{Ra}_\mathrm{crit}} \right)^{b},\; \mathrm{Nu}_\mathrm{min} \right], \qquad \delta = \frac{d}{\mathrm{Nu}}, \qquad q = \frac{k \, \Delta T}{\delta} = \mathrm{Nu} \, \frac{k \, \Delta T}{d}$$

with $\alpha_\mathrm{th}$ the thermal expansivity, $\kappa$ the thermal diffusivity, $\mathrm{Nu}_\mathrm{min}$ the `minimum_nusselt` setting of the `[numerical]` configuration (default 1), and the fitted constants $a$ = `convection_alpha` (default 1.0), $b$ = `convection_beta` (default 1/3), and $\mathrm{Ra}_\mathrm{crit}$ = `critical_rayleigh` (default 1100.0).

The Rayleigh number is the ratio of the buoyancy driving a hot parcel upward to the diffusion bleeding its heat away; above the critical value convection sets in. The Nusselt number is how many times more heat that convection carries than conduction across the whole layer would, so $\mathrm{Nu} = 1$ is conduction. The boundary layer $\delta$ is the conducting thickness that carries the flux across the whole drop $\Delta T$. A layer cooled at its top and heated at its base splits $\Delta T$ between a boundary layer at each end, which carry the same flux, so each is $\delta / 2$ thick. The exponent of 1/3 is the classical boundary-layer result, which has the useful consequence that the convective flux is independent of the layer thickness.

### Magma Ocean

A liquid interior convects in the soft-turbulence regime of a low-viscosity fluid, not as a solid mantle does. With `liquid` set the convection model takes

$$\mathrm{Nu} = a_\mathrm{liquid} \, \mathrm{Ra}^{b_\mathrm{liquid}}$$

($\mathrm{Ra}_\mathrm{crit} = 1$), with $a_\mathrm{liquid}$ = `liquid_convection_alpha` (default 0.089) and $b_\mathrm{liquid}$ = `liquid_convection_beta` (default 1/3) (Solomatov 2000; Lebrun et al. 2013). With a melt's viscosity of order 0.1 Pa s the Rayleigh number reaches $10^{25}$ or more, and the boundary layers are millimeters thick. In a thermal solve the interior is liquid where the Love solve takes it as liquid (see [Convective Reference State](#convective-reference-state)), and the solve reports it as `layer_magma_ocean`.

### Thermal Profiles

During a thermal solve each model builds its layer's profile on every pass, against the structure that pass solved (`build_profile`). The model reads the solved gravity and pressure and the layer's material, with the layer's physics switches, through the network.

- `off` (and a layer with no cooling model): one temperature throughout and no thermal resistance, so the layer conducts perfectly.
- `conduction`: two conducting halves meeting at the mid-radius, where the layer's temperature applies. Each half is a spherical shell of resistance $R = (1 / 4 \pi \bar k)(1/r_a - 1/r_b)$ \[K W$^{-1}$\], where $\bar k$ is the material's conductivity averaged over the temperature between the shell's two ends (the Kirchhoff transform of steady conduction, exact for a conductivity that follows the temperature, such as ice's). The ends are the layer's own temperature and the interface temperature of the last pass. The average is split where the shell's material starts or finishes melting, where its conductivity kinks.
- `convection`: a conducting boundary layer at the base and at the top around an adiabatic interior. The flux law's $\delta$ carries the flux across the whole drop, so each of two boundary layers takes $\delta / 2$, and a layer whose base carries no heat (the innermost layer, or one above a layer outside the network) has only the upper one, of thickness $\delta$. Each is at most 40 percent of the layer. The layer's temperature applies at the top of the interior, the upper-mantle temperature of parameterized convection (Stevenson et al. 1983), and the interior warms downward from it along the adiabat $dT/dr = -(\alpha + \alpha_L) g T / c_p$, with the material's expansivity $\alpha$, latent expansion $\alpha_L$ (see [Latent Heat](../Material/materials.md#latent-heat)), and heat capacity $c_p$. The network finds the base temperature by marching that adiabat down from the top with the material at the temperature reached, so a melting range, which buffers the adiabat, is followed. The march takes Runge-Kutta steps with step doubling, refining wherever the adiabat kinks: where the material starts or finishes melting, and where a melting curve changes branch inside a melting range (the Monteux et al. 2016 peridotite curves do at 20 GPa). The boundary layers' resistances use the same temperature-averaged conductivity as conduction. A converged solve's integrated profile lands on the layer's temperature at the top of its interior and on the surface temperature at the top of the world.

The temperature drop of a convecting layer is the sum of the drops across its two boundary layers: from the top of the layer below (the end of its adiabat, when that layer convects) to the base of this layer's adiabat, and from the layer's temperature to the base of the layer above or to the surface temperature. The adiabat's base comes from the previous pass. The network joins every layer's profile into a chain of resistances and integrates the temperature and heat flow through the planet (see [Temperature and Heat Flow](../Structures/worlds/worlds.md#temperature-and-heat-flow)).

### Convective Reference State

The convection model evaluates the viscosity of its Rayleigh number at its own reference point: the top of its adiabatic interior, under its upper boundary layer, at the layer's own temperature and the solved pressure there. The boundary layer's thickness follows from that viscosity, so the model iterates the two to agreement on each pass's structure, starting from the previous pass's boundary layer (the layer's top before there is one). Near the melt onset, where the viscosity follows the pressure through the melt fraction, this saves the solve tens of passes. That is where the layer's temperature applies, so the viscosity is the one the interior actually has where it is coolest. The gravity, the density, and the thermal constants are the mid-layer values at the layer's temperature, and the diffusivity is $k / (\rho c_p)$ there.

The interior is liquid at the reference point, and takes the magma-ocean scaling, where the Love solve takes it as liquid: everywhere in a liquid layer (a liquid-only material, or `state = "liquid"`), and, in a layer that can change state, where the material is fully molten or its post-melt shear modulus is at or below `[numerical] minimum_solid_rigidity` times the world's $\rho g R$. A melting mantle past the critical melt fraction of its weakening law is therefore a magma ocean even before it is fully molten. The solve reports the reference point as `layer_reference_pressure`, `layer_reference_viscosity`, and `layer_reference_melt_fraction` (NaN for a layer that does not convect).

### Behavior at the Limits

At the default floor of 1 a sub-critical or rigid layer ($\mathrm{Ra} < \mathrm{Ra}_\mathrm{crit}$, or an infinite viscosity) returns the conductive flux $k \, \Delta T / d$, the same as `ConductiveCooling`. Degenerate inputs take a defined path rather than producing a division by zero. A non-positive temperature drop, or a thickness below the shared `minimum_layer_thickness` configuration floor, sets the Rayleigh number to zero and the Nusselt number to $\mathrm{Nu}_\mathrm{min}$, so the boundary layer is $d / \mathrm{Nu}_\mathrm{min}$ (the whole layer below the minimum thickness). The flux is then zero for a zero temperature drop, but the boundary layer still sets the resistance between the layer and its neighbors in a thermal solve.

A NaN input other than the drop or the thickness (usually a NaN viscosity, from a phase with no viscosity law) passes through to a NaN result. In a thermal solve a convecting layer whose flux law gives no usable boundary layer takes the 40 percent limit for each, logs a warning, and reports `layer_boundary_fallback`. A liquid layer with a convection model is a magma ocean throughout, and its liquid phase needs a viscosity law. A well-mixed liquid (an ocean or a liquid core) can instead be left with no cooling model, which holds it at one temperature.

### Choosing a Model

- `convection`: a mantle or a thick ice shell that is warm enough to flow. Give its material a viscosity law that depends on temperature.
- `conduction`: a stagnant lid, a crust, a thin ice shell, or a layer known to sit below the critical Rayleigh number.
- `off`, or no cooling model: a layer whose internal gradient does not matter, a well-mixed liquid (an ocean or a liquid core), or a layer whose temperature is set by hand.

## Python API

```python
from TidalPy.Cooling import ConvectiveCooling, make_cooling

convective_cooling = ConvectiveCooling(
    convection_alpha=1.0,
    convection_beta=1.0 / 3.0,
    critical_rayleigh=1100.0)

# delta_temp, thickness, gravity, density, viscosity, conductivity, diffusivity, expansion
result = convective_cooling.calc_cooling(1000.0, 1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)

# The same layer as a magma ocean, with a melt's viscosity
magma_ocean = convective_cooling.calc_cooling(
    1000.0, 1.0e6, 9.8, 3300.0, 0.1, 4.0, 1.0e-6, 3.0e-5,
    liquid=True)
print(result.boundary_layer_thickness, magma_ocean.boundary_layer_thickness)   # [m]

# Name factory: case-insensitive, aliases accepted.
conductive_cooling = make_cooling("conductive")
```

Constructors take the model's parameters by argument name or config key, as keywords or positionally in the order `get_parameter_info()` lists them. `make_cooling(model_name, config=None)` resolves a name or alias case-insensitively and builds the model from `config`. Absent keys take the model's defaults. An unknown name, a key the model does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted one.

| Parameter | Config key | Default | Model |
|---|---|---|---|
| `convection_alpha` | `convection_alpha` | 1.0 | Convection, solid interior |
| `convection_beta` | `convection_beta` | 1/3 | Convection, solid interior |
| `critical_rayleigh` | `critical_rayleigh` | 1100.0 | Convection, solid interior |
| `liquid_convection_alpha` | `liquid_convection_alpha` | 0.089 | Convection, liquid interior |
| `liquid_convection_beta` | `liquid_convection_beta` | 1/3 | Convection, liquid interior |

The parameters read as attributes, and every model has `model_name`, `parameters`, `get_parameter(name)`, `get_parameter_info()`, `with_parameters(**changes)`, `get_config_dict()`, and `save_config(path)`. `off` and `conduction` have no parameters. `cooling_model_names()` lists the canonical names and `cooling_config_keys(name)` the keys one model reads.

### Factory Internals

At the C++ level `c_cooling_registry()` lists each model's names (canonical first, then aliases), its binary class id, and its constructor. `c_find_cooling(name, params)` builds a model from a `c_ParamMap` and throws `std::invalid_argument` for an unknown name or parameter, and `c_cooling_canonical_name` and `c_cooling_model_names` resolve names. The Python `make_cooling` passes the config dict through and wraps the result in the matching class.

### Vectorized Evaluation

Two of the eight inputs change as a layer evolves thermally: the temperature drop and the viscosity. The vectorized methods, defined once on the base class so every model inherits them, sweep exactly those and hold the rest fixed. Each also takes `liquid`.

| Method | Sweeps |
|---|---|
| `calc_cooling_vectorize_temperature(delta_temp[], ...)` | A temperature-drop array at fixed viscosity. |
| `calc_cooling_vectorize_viscosity(..., viscosity[], ...)` | A viscosity array at fixed temperature drop. |
| `calc_cooling_vectorize_all(delta_temp[], ..., viscosity[], ...)` | Two equal-length arrays, element by element. |

```python
import numpy as np
from TidalPy.Cooling import ConvectiveCooling

model = ConvectiveCooling()
# delta_temp, thickness, gravity, density, viscosity, conductivity, diffusivity, expansion
sweep = model.calc_cooling_vectorize_temperature(
    np.linspace(100.0, 2000.0, 50),
    1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5
)
flux_profile = sweep.cooling_flux   # float64 array
```

At the C++ level `calc_cooling_vectorize(delta_temp, viscosity, base_inputs, num_points, ...)` broadcasts the two vectors (each one value per point or a single value used at every point) and writes into four caller-supplied buffers. Mismatched lengths raise `ValueError` from Python. The Cython wrappers return one `CoolingResult` whose four fields are arrays.

### Convenience Functions

For a one-shot evaluation that does not leave an object behind.

```python
import numpy as np
from TidalPy.Cooling import conductive, convective, cooling_off

result = convective(1000.0, 1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)

# delta_temp and viscosity may be floats or arrays and are broadcast together;
# the remaining inputs are scalar constants.
sweep = convective(
    np.linspace(100.0, 2000.0, 50),
    1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)

# A liquid interior with a stiffer scaling
magma_ocean = convective(
    1000.0, 1.0e6, 9.8, 3300.0, 0.1, 4.0, 1.0e-6, 3.0e-5,
    liquid=True,
    liquid_convection_alpha=0.1)
```

The signatures are `cooling_off(delta_temp, thickness)`, `conductive(delta_temp, thickness, thermal_conductivity)`, and `convective(delta_temp, thickness, gravity, density, viscosity, thermal_conductivity, thermal_diffusivity, thermal_expansion, convection_alpha=None, convection_beta=None, critical_rayleigh=None, liquid=False, liquid_convection_alpha=None, liquid_convection_beta=None)`, where `None` takes the model's default. Each builds a C++ model for the one call, runs the shared broadcasting loop over the inputs (a scalar input is held fixed), and returns a `CoolingResult`. The off and conduction functions take only the inputs they actually use, which is why their argument lists are shorter than `calc_cooling`'s fixed eight.

### Attaching a Model to a `Layer`

```python
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import TerrestrialWorld

core = Layer(
    "core",
    0,
    0.0,
    9.0e5,
    material="simple_iron_core",
    temperature=1900.0)                  # No cooling model: isothermal
mantle = Layer(
    "mantle",
    1,
    9.0e5,
    1.82e6,
    material="peridotite",
    temperature=1650.0,                  # [K] at the top of its adiabatic interior
    use_melting=True,
    use_pressure_melting=True,
    cooling="convection")                # A model, a model name, or a config table

world = TerrestrialWorld("Io-like", 1.82e6, 8.93e22)
world.add_layer(core)
world.add_layer(mantle)
result = world.solve_eos(
    solve_temperature=True,
    surface_temperature=110.0)           # [K] what the outermost layer radiates to

print(result["layer_nusselt_number"])        # The mantle convects
print(result["layer_boundary_thickness"])    # [m] each of the mantle's boundary layers
print(result["layer_reference_viscosity"])   # [Pa s] at the top of the mantle's interior
print(result["layer_magma_ocean"])           # Whether that interior is liquid

world.mantle.cooling = {"model": "convection", "critical_rayleigh": 1000.0}   # Replace it
world.mantle.cooling = None                                                   # Remove it
```

The layer shares the model rather than copying it, so one model can serve several layers. Setting `cooling` on a layer of a solved world makes the world forget its solve. The declarative form is a `[layers.<name>.cooling]` table in a world's TOML, keyed by `model` plus any parameters. See the [TOML schema](../Structures/config/toml_schema.md) and [Layer](../Structures/layers/layer.md).

## C++ API

The models are in `TidalPy/Cooling/cooling_.hpp` and their base and the profile interface in `cooling_base_.hpp` (namespace `tidalpy`, header only).

| Type | Contents |
|---|---|
| `c_CoolingInputs` | The flux law's inputs (see [Inputs and Result](#inputs-and-result)). |
| `c_CoolingResult` | `cooling_flux`, `blt`, `rayleigh_number`, `nusselt_number`. |
| `c_LayerThermalContext` | A layer during one pass: its radii, its own temperature, the temperatures across its lower and upper boundary layers, the base temperature, boundary-layer thickness, and interface temperatures (`inner_node_temperature`, `outer_node_temperature`) of the last pass, and whether its base is insulated. |
| `c_LayerThermalProbe` | The network's read access: `calc_structure(radius, gravity, pressure)` of the solved structure, `calc_transport_state(pressure, temperature, radius, out)` of the layer's material (a `c_TransportState`, liquid where the Love solve takes it as liquid), `calc_shell_resistance(radius_inner, radius_outer, temperature_inner, temperature_outer)` (the temperature-averaged resistance of a conducting shell), and `calc_adiabat_base_temperature(radius_lower, radius_upper, top_temperature)` (the adiabat marched down from its top). |
| `c_LayerThermalProfile` | What a model returns: the boundary-layer thickness, the bottom and top resistances, the base temperature, Ra and Nu, the fallback and magma-ocean flags, the material values it used, and the reference state. |

`c_CoolingBase` declares `get_temperature_kind()` (`Isothermal`, `Conductive`, or `Adiabatic`, known before any structure is solved), `calc_cooling(inputs)`, `build_profile(context, probe, out)`, and `calc_cooling_vectorize(...)`, which each model implements through the shared `p_vectorize_kernel` so its flux law inlines into the loop. The flux laws are also free functions, `cool_off(inputs)`, `cool_conduction(inputs)`, and `cool_convection(inputs, alpha, beta, critical_rayleigh)`. `c_shell_conductivity(radius_inner, radius_outer, resistance, fallback)` gives the uniform conductivity of a shell with a given resistance \[K W$^{-1}$\], and `d_MAX_BOUNDARY_FRACTION` (0.4) caps a boundary layer. `c_find_cooling(name, params)` and `c_cooling_from_binary(stream, force)` build a model by name and rebuild one from a binary record. The world's thermal network (`Structures/worlds/thermal_layout_.hpp`) implements the probe and assembles the profiles; it does not switch on the model type.

## Serialization

| Call | Result |
|---|---|
| `get_config_dict()` | `model` plus any model parameters. |
| `save_config(path)` | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path)` | The model's TidalPy binary record, its parameters written by key. |

A cooling model on a layer is saved and restored with the layer, in its binary record and in its config dict (as its `cooling` table).

## Adding a New Cooling Model

To add a cooling model named `Foo`:

**C++ (`TidalPy/Cooling/cooling_.hpp`)**

1. Add a free function `cool_foo(const c_CoolingInputs&, ...)` for the flux law. Guard any denominator with `c_guard_denominator` (`constants_.hpp`), the shared floor.
2. Add `c_FooCooling : public c_SpecModel<c_FooCooling, c_CoolingBase>`: its `parameter_specs()` table, `C_CLASS_ID`, two constructors that call `p_initialize`, `get_temperature_kind()`, `calc_cooling` (calling `cool_foo`), `build_profile`, and `calc_cooling_vectorize` through `p_vectorize_kernel`. A model whose profile is not one of the three kinds needs a new `c_TemperatureKind` (`Utilities/classes/thermo_point_.hpp`) and its handling in the thermal network.
3. Add one row to `c_cooling_registry()` with the model's names and aliases.

**C++ (`TidalPy/Utilities/binary/binary_.hpp`)**

4. Reserve a unique `BinaryClassID::FooCooling`. Cooling models occupy the 40X range.

**Cython (`cooling.pyx`)**

5. Add `cdef class FooCooling(CoolingBase)` with a docstring and `MODEL_NAME = "foo"`, include it in the `ModelFamily` list, and add a lower-case `foo(...)` convenience function.

**Package, tests, and documentation**

6. Export `FooCooling` and `foo` from `TidalPy/Cooling/__init__.py`.
7. Add its flux law and profile tests to `Tests/Test_Cooling/test_cooling_01.py` and, for a profile, to the world thermal tests in `Tests/Test_Structures/Test_Worlds/`. The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover its parameters, config, binary record, and errors without changes.
8. Document the model here with its formula, parameters, and references.

No build-system change is needed; `Cooling.cooling` is already registered in `cython_extensions.json`.

## References

- Turcotte, D. L., and Schubert, G. (2002). *Geodynamics*, second edition. Rayleigh and Nusselt convection scaling, and conduction.
- Solomatov, V. S. (1995). Scaling of temperature- and stress-dependent viscosity convection. *Physics of Fluids*, 7(2), 266-274. Stagnant-lid scaling.
- Schubert, G., Turcotte, D. L., and Olson, P. (2001). *Mantle Convection in the Earth and Planets*. Boundary-layer theory.
- Stevenson, D. J., Spohn, T., and Schubert, G. (1983). Magnetism and thermal evolution of the terrestrial planets. *Icarus*, 54(3), 466-489. The upper-mantle temperature at the top of the adiabatic interior.
- Solomatov, V. S. (2000). Fluid dynamics of a terrestrial magma ocean. In *Origin of the Earth and Moon*, University of Arizona Press, 323-338. The soft-turbulence scaling of a liquid interior.
- Lebrun, T., Massol, H., Chassefière, E., Davaille, A., Marcq, E., Sarda, P., Leblanc, F., and Brandeis, G. (2013). Thermal evolution of an early magma ocean in interaction with the atmosphere. *Journal of Geophysical Research: Planets*, 118(6), 1155-1176. The same scaling applied to a cooling magma ocean.
