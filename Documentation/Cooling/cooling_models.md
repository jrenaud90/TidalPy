# Cooling Models (`Cooling`)

_Updated: 2026-10-07_

A cooling model says how heat moves through a layer. Its flux law maps the layer's physical state onto a **cooling result**: the surface heat flux $q$ [W m$^{-2}$], the thermal boundary-layer thickness [m], and the Rayleigh and Nusselt numbers. In a world's thermal solve the same model builds its layer's temperature profile. Each `Layer` can hold one cooling model.

| Model | Aliases | Python class | Use for |
|---|---|---|---|
| `off` | `none` | `OffCooling` | A layer whose internal gradient does not matter, a well-mixed liquid (an ocean or a liquid core), or a layer whose temperature is set by hand. No cooling model behaves the same. |
| `conduction` | `conductive` | `ConductiveCooling` | A stagnant lid, a crust, a thin ice shell, or a layer known to sit below the critical Rayleigh number. |
| `convection` | `convective` | `ConvectiveCooling` | A mantle or a thick ice shell warm enough to flow. Give its material a temperature-dependent viscosity law. |

The canonical names (first column) are the ones written to a config dict.

## Quick Example

```python
from TidalPy.Cooling import ConvectiveCooling, make_cooling

convective_cooling = ConvectiveCooling(
    convection_alpha=1.0,
    convection_beta=1.0 / 3.0,
    critical_rayleigh=1100.0)

# delta_temp, thickness, gravity, density, viscosity, conductivity, diffusivity, expansion
result = convective_cooling.calc_cooling(1000.0, 1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)
flux, boundary_layer, rayleigh, nusselt = result

# The same layer as a magma ocean, with a melt's viscosity
magma_ocean = convective_cooling.calc_cooling(
    1000.0, 1.0e6, 9.8, 3300.0, 0.1, 4.0, 1.0e-6, 3.0e-5,
    liquid=True)
print(result.boundary_layer_thickness, magma_ocean.boundary_layer_thickness)   # [m]

# Name factory: case-insensitive, aliases accepted.
conductive_cooling = make_cooling("conductive")
```

The eight positional inputs are described in [Inputs and Result](#inputs-and-result).

## Attaching a Model to a `Layer`

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
print(result["layer_reference_viscosity"])   # [Pa s] over the top of the mantle
print(result["layer_magma_ocean"])           # Whether that interior is liquid

world.mantle.cooling = {"model": "convection", "critical_rayleigh": 1000.0}   # Replace it
world.mantle.cooling = None                                                   # Remove it
```

The layer shares the model rather than copying it, so one model can serve several layers. Setting `cooling` on a layer of a solved world makes the world forget its solve. In a world's TOML the model is a `[layers.<name>.cooling]` table keyed by `model` plus any parameters. See the [TOML schema](../Structures/config/toml_schema.md) and [Layer](../Structures/layers/layer.md).

## Current Models

| Model | Cooling flux $q$ [W m$^{-2}$] | Boundary layer | Ra | Nu | Profile in a thermal solve |
|---|---|---|---|---|---|
| `OffCooling` | $0$ | none (NaN) | 0 | 1 | Isothermal |
| `ConductiveCooling` | $k \, \Delta T / d$ | thickness $d$ | 0 | 1 | Two conducting halves |
| `ConvectiveCooling` | $k \, \Delta T / \delta$ | $\delta = d / \mathrm{Nu}$ | see below | see below | Boundary layers around an adiabatic interior |

`ConductiveCooling` transports heat by diffusion. `ConvectiveCooling` is the parameterized boundary-layer model; its boundary layer is what makes convection so much more effective than conduction, squeezing the same temperature drop across a thin layer instead of the whole interior.

### Convection

$$\mathrm{Ra} = \frac{\alpha_\mathrm{th} \rho g \, \Delta T \, d^3}{\eta \kappa}, \qquad \mathrm{Nu} = \max\left[ a \left( \frac{\mathrm{Ra}}{\mathrm{Ra}_\mathrm{crit}} \right)^{b},\; \mathrm{Nu}_\mathrm{min} \right], \qquad \delta = \frac{d}{\mathrm{Nu}}, \qquad q = \frac{k \, \Delta T}{\delta} = \mathrm{Nu} \, \frac{k \, \Delta T}{d}$$

with $\alpha_\mathrm{th}$ the thermal expansivity, $\kappa$ the thermal diffusivity, $\mathrm{Nu}_\mathrm{min}$ the `[numerical]` setting `minimum_nusselt` (default 1), and the fitted constants $a$ = `convection_alpha` (default 1.0), $b$ = `convection_beta` (default 1/3), and $\mathrm{Ra}_\mathrm{crit}$ = `critical_rayleigh` (default 1100.0).

The Rayleigh number is the ratio of the buoyancy driving a hot parcel upward to the diffusion bleeding its heat away; above the critical value convection sets in. The Nusselt number is how many times more heat convection carries than conduction across the whole layer would, so $\mathrm{Nu} = 1$ is conduction. The boundary layer $\delta$ is the conducting thickness that carries the flux across the whole drop $\Delta T$. A layer cooled at its top and heated at its base splits $\Delta T$ between two boundary layers carrying the same flux, each $\delta / 2$ thick. The exponent 1/3 is the classical boundary-layer result, which makes the convective flux independent of the layer thickness.

### Magma Ocean

A liquid interior convects in the soft-turbulence regime of a low-viscosity fluid, not as a solid mantle does. With `liquid` set the convection model takes

$$\mathrm{Nu} = a_\mathrm{liquid} \, \mathrm{Ra}^{b_\mathrm{liquid}}$$

($\mathrm{Ra}_\mathrm{crit} = 1$), with $a_\mathrm{liquid}$ = `liquid_convection_alpha` (default 0.089) and $b_\mathrm{liquid}$ = `liquid_convection_beta` (default 1/3) (Solomatov 2000; Lebrun et al. 2013). With a melt's viscosity of order 0.1 Pa s the Rayleigh number reaches $10^{25}$ or more, and the boundary layers are millimeters thick. The solve reports a liquid interior as `layer_magma_ocean` (see [Convective Reference State](#convective-reference-state)).

### Thermal Profiles

On every pass of a thermal solve each model builds its layer's profile against the structure that pass solved, reading the solved gravity and pressure and the layer's material.

- `off` (and no cooling model): one temperature throughout and no thermal resistance.
- `conduction`: two conducting halves meeting at the mid-radius, where the layer's temperature applies. Each half is a spherical shell of resistance $R = (1 / 4 \pi \bar k)(1/r_a - 1/r_b)$ \[K W$^{-1}$\], with $\bar k$ the material's conductivity averaged over the temperature between the shell's ends (the layer's temperature and the last pass's interface temperature). This Kirchhoff transform of steady conduction is exact for a conductivity that follows the temperature, such as ice's. The average is split where the material starts or finishes melting, where its conductivity kinks.
- `convection`: conducting boundary layers at the base and the top, each $\delta / 2$ thick, around an adiabatic interior. A layer whose base carries no heat (the innermost layer, or one above a layer outside the network) has only the upper one, $\delta$ thick. Each is at most 40 percent of the layer and uses the same temperature-averaged conductivity as conduction. The layer's temperature applies at the top of the interior, the upper-mantle temperature of parameterized convection (Stevenson et al. 1983), and the interior warms downward along the adiabat $dT/dr = -(\alpha + \alpha_L) g T / c_p$, with the material's expansivity $\alpha$, latent expansion $\alpha_L$ (see [Latent Heat](../Material/materials.md#latent-heat)), and heat capacity $c_p$ at the temperature reached. The march follows a melting range and refines its steps where the adiabat kinks: where the material starts or finishes melting, and where a melting curve changes branch (the Monteux et al. 2016 peridotite curves do at 20 GPa).

A convecting layer's temperature drop is the sum of the drops across its two boundary layers. Only an unstable drop (hotter below) drives convection: a boundary layer the wrong way round, such as a core colder than the mantle base above it, adds nothing. A converged solve lands on the layer's temperature at the top of its interior and on the surface temperature at the top of the world. How the network joins the profiles and integrates the temperature and heat flow is in [Temperature and Heat Flow](../Structures/worlds/worlds.md#temperature-and-heat-flow).

### Convective Reference State

The convection model takes the viscosity of its Rayleigh number over the top of its layer: the logarithmic mean of the material's viscosity at the layer's temperature and the solved pressure, at eight depths spread evenly through the top `viscosity_depth_fraction` of the layer (default 0.05, about 145 km for Earth's mantle). The middle of that range is the reference point the solve reports and the point that decides whether the interior is liquid.

The range is fixed by the layer's geometry on purpose. A viscosity taken where the boundary layer ends would feed back on that layer's thickness, and near melt onset that loop has two self-consistent solutions between which the heat flow can jump. With the range fixed, the heat flow follows the temperature continuously, which a time integration needs.

The interior is liquid at the reference point, and takes the magma-ocean scaling, where the Love solve takes it as liquid: everywhere in a liquid layer (a liquid-only material, or `state = "liquid"`), and, in a layer that can change state, where the material is fully molten or its post-melt shear modulus is at or below `[numerical] minimum_solid_rigidity` times the world's $\rho g R$. With the default threshold of 0, a melting mantle past the end of its weakening law's breakdown, where its shear modulus vanishes, is therefore a magma ocean before it is fully molten. The solve reports `layer_reference_pressure`, `layer_reference_viscosity`, and `layer_reference_melt_fraction` (NaN for a layer that does not convect).

## Inputs and Result

A flux-law evaluation takes eight physical quantities and one flag, in this order in the Python `calc_cooling` (in C++ they are bundled into a `c_CoolingInputs` struct).

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

The result is a `CoolingResult` with `cooling_flux` [W m$^{-2}$], `boundary_layer_thickness` [m], `rayleigh`, and `nusselt`: Python floats for a scalar evaluation, `float64` arrays for a vectorized one. It supports `to_dict()` and unpacks in field order.

## Python API

Constructors take the model's parameters by argument name or config key, as keywords or positionally in `get_parameter_info()` order. `make_cooling(model_name, config=None)` resolves a name or alias case-insensitively and builds the model from `config`; absent keys take the defaults. An unknown name, a key the model does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted one.

| Parameter | Config key | Default | Model |
|---|---|---|---|
| `convection_alpha` | `convection_alpha` | 1.0 | Convection, solid interior |
| `convection_beta` | `convection_beta` | 1/3 | Convection, solid interior |
| `critical_rayleigh` | `critical_rayleigh` | 1100.0 | Convection, solid interior |
| `liquid_convection_alpha` | `liquid_convection_alpha` | 0.089 | Convection, liquid interior |
| `liquid_convection_beta` | `liquid_convection_beta` | 1/3 | Convection, liquid interior |
| `viscosity_depth_fraction` | `viscosity_depth_fraction` | 0.05 | Convection: the top fraction of the layer the interior viscosity is averaged over |

`off` and `conduction` have no parameters. Parameters read as attributes, and every model has `model_name`, `parameters`, `get_parameter(name)`, `get_parameter_info()`, `with_parameters(**changes)`, `get_config_dict()` (`model` plus any parameters), `save_config(path)`, and `save_binary(path)` / `load_binary(path)`. A cooling model is saved and restored with its layer. `cooling_model_names()` lists the canonical names, `canonical_cooling_name(name)` resolves any name or alias, and `cooling_config_keys(name)` gives the keys one model reads.

### Vectorized Evaluation

Two inputs change as a layer evolves thermally: the temperature drop and the viscosity. Every model has three methods that sweep those and hold the rest fixed; each also takes `liquid`.

| Method | Sweeps |
|---|---|
| `calc_cooling_vectorize_temperature(delta_temp[], ...)` | A temperature-drop array at fixed viscosity. |
| `calc_cooling_vectorize_viscosity(..., viscosity[], ...)` | A viscosity array at fixed temperature drop. |
| `calc_cooling_vectorize_all(delta_temp[], ..., viscosity[], ...)` | Two equal-length arrays, element by element. |

They return one `CoolingResult` of arrays. Mismatched lengths raise `ValueError`.

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

### Convenience Functions

One-shot functions that leave no model object behind:

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

The signatures are `cooling_off(delta_temp, thickness)`, `conductive(delta_temp, thickness, thermal_conductivity)`, and `convective(delta_temp, thickness, gravity, density, viscosity, thermal_conductivity, thermal_diffusivity, thermal_expansion, convection_alpha=None, convection_beta=None, critical_rayleigh=None, liquid=False, liquid_convection_alpha=None, liquid_convection_beta=None)`, where `None` takes the model's default. Each returns a `CoolingResult`.

## Limits and Failure Modes

- At the default floor of 1, a sub-critical or rigid layer ($\mathrm{Ra} < \mathrm{Ra}_\mathrm{crit}$, or an infinite viscosity) returns the conductive flux $k \, \Delta T / d$, the same as `ConductiveCooling`.
- A non-positive temperature drop, or a thickness below the shared `minimum_layer_thickness` floor, sets Ra to zero and Nu to $\mathrm{Nu}_\mathrm{min}$, so the boundary layer is $d / \mathrm{Nu}_\mathrm{min}$ (the whole layer below the minimum thickness). The flux is then zero for a zero drop, but the boundary layer still sets the resistance to the neighbors in a thermal solve.
- Any other NaN input (usually a NaN viscosity, from a phase with no viscosity law) gives a NaN result.
- In a thermal solve, a convecting layer whose flux law gives no usable boundary layer takes the 40 percent limit for each, logs a warning, and reports `layer_boundary_fallback`.
- A liquid layer with a convection model is a magma ocean throughout, and its liquid phase needs a viscosity law. A well-mixed liquid (an ocean or a liquid core) can instead have no cooling model, which holds it at one temperature.

## C++ API

The models are in `TidalPy/Cooling/cooling_.hpp`, and their base and the profile interface in `cooling_base_.hpp` (namespace `tidalpy`, header only).

| Type | Contents |
|---|---|
| `c_CoolingInputs` | The flux law's inputs (see [Inputs and Result](#inputs-and-result)). |
| `c_CoolingResult` | `cooling_flux`, `blt`, `rayleigh_number`, `nusselt_number`. |
| `c_LayerThermalContext` | A layer's radii, temperatures, and last-pass state during one thermal pass. |
| `c_LayerThermalProbe` | The thermal network's read access to the solved structure and the layer's material, the shell resistance, and the adiabat's base temperature. |
| `c_LayerThermalProfile` | What `build_profile` returns: boundary-layer thickness, resistances, base temperature, Ra and Nu, flags, and the reference state. |

The abstract `c_CoolingBase` declares `get_temperature_kind()` (`Isothermal`, `Conductive`, or `Adiabatic`), `calc_cooling(inputs)`, `build_profile(context, probe, out)`, and `calc_cooling_vectorize(delta_temp, viscosity, base_inputs, num_points, ...)`, which broadcasts the two vectors (one value per point, or one for all) into four caller-supplied buffers. The concrete classes are `c_OffCooling`, `c_ConductiveCooling`, and `c_ConvectiveCooling`; the flux laws are also free functions `cool_off(inputs)`, `cool_conduction(inputs)`, and `cool_convection(inputs, alpha, beta, critical_rayleigh)`. `c_shell_conductivity(radius_inner, radius_outer, resistance, fallback)` gives the uniform conductivity of a shell with a given resistance \[K W$^{-1}$\]. `c_find_cooling(name, params)` builds a model by name (throwing `std::invalid_argument` for an unknown name or parameter), `c_cooling_from_binary(stream, force)` rebuilds one, and `c_cooling_canonical_name` and `c_cooling_model_names` resolve names.

### Adding a New Cooling Model

A model `Foo` is a flux-law function `cool_foo(const c_CoolingInputs&, ...)` (guard denominators with `c_guard_denominator`) and a `c_SpecModel<c_FooCooling, c_CoolingBase>` in `cooling_.hpp` implementing the base's four methods. A profile that is not one of the three kinds also needs a new `c_TemperatureKind` and thermal-network support. Then add a registry row, a `BinaryClassID` in the 40X range, a Cython class and `foo(...)` function exported from `TidalPy/Cooling/__init__.py`, tests in `Tests/Test_Cooling/`, and a section here.

## References

- Turcotte, D. L., and Schubert, G. (2002). *Geodynamics*, second edition. Rayleigh and Nusselt convection scaling, and conduction.
- Solomatov, V. S. (1995). Scaling of temperature- and stress-dependent viscosity convection. *Physics of Fluids*, 7(2), 266-274. Stagnant-lid scaling.
- Schubert, G., Turcotte, D. L., and Olson, P. (2001). *Mantle Convection in the Earth and Planets*. Boundary-layer theory.
- Stevenson, D. J., Spohn, T., and Schubert, G. (1983). Magnetism and thermal evolution of the terrestrial planets. *Icarus*, 54(3), 466-489. The upper-mantle temperature at the top of the adiabatic interior.
- Solomatov, V. S. (2000). Fluid dynamics of a terrestrial magma ocean. In *Origin of the Earth and Moon*, University of Arizona Press, 323-338. The soft-turbulence scaling of a liquid interior.
- Lebrun, T., Massol, H., Chassefière, E., Davaille, A., Marcq, E., Sarda, P., Leblanc, F., and Brandeis, G. (2013). Thermal evolution of an early magma ocean in interaction with the atmosphere. *Journal of Geophysical Research: Planets*, 118(6), 1155-1176. The same scaling applied to a cooling magma ocean.
