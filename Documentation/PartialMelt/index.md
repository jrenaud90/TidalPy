# Partial Melting (`PartialMelt`)

_Updated: 2026-10-01_

`TidalPy.PartialMelt` contains the melting laws a [material](../Material/materials.md) uses once it begins to melt. Melting curves map a pressure \[Pa\] onto a solidus or liquidus temperature \[K\]. Melt-weakening laws map the melt fraction, the temperature, and both phases' shear moduli \[Pa\] and viscosities \[Pa s\] onto those of the partially molten aggregate. Bulk-mixing laws do the same for the bulk modulus and the bulk viscosity. A material combines two melting curves, one weakening law, and optionally the mixing laws with its solid and liquid phases, and its melt fraction runs linearly between the curves.

Partial melt occurs where some minerals of an assemblage melt while others remain solid. For many materials there is a critical melt fraction where the material switches from behaving like a solid with pockets of melt to a liquid carrying chunks of solid (often between 0.4 and 0.6 for rocks and ices), and the shear modulus and viscosity fall by orders of magnitude across it.

Partial melting is an important feedback that makes solid-body tidal heating self-limiting. Heating raises the temperature, the temperature raises the melt fraction, the melt fraction drops the viscosity and shear modulus by orders of magnitude, and a weaker body deforms more but dissipates less once its Maxwell time falls far below the forcing period. Whether a body runs away to a magma ocean or settles into a warm steady state is largely decided by the shape of the weakening curve near the critical melt fraction, which is why the models differ most sharply there.

| Page | Covers |
|---|---|
| [Melting Laws](partial_melt_models.md) | The melting curves and their pressure slopes, the melt-weakening laws, the bulk-mixing laws, their parameters, the Python and C++ surfaces, and how to add a law. |

```{toctree}
:maxdepth: 1

Melting Laws <partial_melt_models.md>
```

## Where Partial Melt is Used

A material holds its melting laws in a `melting` table beside its solid and liquid phases (see [Phases and Materials](../Material/materials.md)), and a layer applies them when it sets `use_melting` (and `use_pressure_melting` for curves that follow the pressure). During the whole-planet equation-of-state solve the material evaluates its phases at each point, finds the melt fraction from the curves at the local pressure, and passes both phases' values through the weakening and mixing laws. The post-melt values are what [`Rheology`](../Rheology/index.md) turns into a complex modulus. The solve splits a melting layer into solid and liquid zones where its post-melt shear modulus crosses the radial solver's liquid threshold, and a convecting layer whose interior is liquid takes the magma-ocean scaling of its [cooling model](../Cooling/cooling_models.md).

Like viscosity, melt weakening is frequency-independent and therefore resolved once per equation-of-state solve rather than once per tidal mode.

## Examples

`Demos/Physics/10_thermal_eos.ipynb` sweeps a melting mantle's temperature through its melting range, and `Demos/Systems/12_thermal_orbital_evolution.ipynb` follows a melting mantle through a coupled thermal-orbital evolution.

## References

- Fischer, H.-J., and Spohn, T. (1990). Thermal-orbital histories of viscoelastic models of Io. *Icarus*, 83(1), 39-65.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Monteux, J., Andrault, D., and Samuel, H. (2016). On the cooling of a deep terrestrial magma ocean. *Earth and Planetary Science Letters*, 448, 140-149.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98.
- Simon, F., and Glatzel, G. (1929). Bemerkungen zur Schmelzdruckkurve. *Zeitschrift für anorganische und allgemeine Chemie*, 178(1), 309-316.
