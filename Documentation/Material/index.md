# Materials (`Material`)

_Updated: 2026-10-01_

`TidalPy.Material` contains functionality to calculate the state of planet-relevant materials. A material maps a pressure \[Pa\], temperature \[K\], and radius \[m\] onto its density \[kg m$^{-3}$\], isothermal and adiabatic bulk moduli \[Pa\], static shear modulus \[Pa\], shear and bulk viscosities \[Pa s\], thermal expansivity \[K$^{-1}$\], heat capacity \[J kg$^{-1}$ K$^{-1}$\], thermal conductivity \[W m$^{-1}$ K$^{-1}$\], and melt fraction. The whole-planet solve integrates the density from the center outward to find the body's radial structure (gravity, pressure, mass, and moment of inertia), and every later calculation reads the material's other properties from the same evaluation.

A material is built from small laws, each a model of its own family: an equation-of-state law and a shear-modulus law (this module), viscosity laws ([`Viscosity`](../Viscosity/index.md)), default rheologies ([`Rheology`](../Rheology/index.md)), and melting curves, melt weakening, and bulk mixing ([`PartialMelt`](../PartialMelt/index.md)). A `Phase` combines the first four with its thermal constants, and a `Material` combines a solid phase, a liquid phase, or both with the melting laws. MatPack ships named materials built this way.

| Page | Covers |
|---|---|
| [Equation-of-State and Shear-Modulus Laws](material_eos.md) | The equation-of-state laws, their thermal terms and pressure inversion, the shear-modulus laws, the factories, serialization, the C++ surface, and how to add a law. |
| [Phases and Materials](materials.md) | `Phase` and `Material`, the state a material reports, how it melts, the physics switches a layer passes it, changing a material, and serialization. |
| [MatPack](matpack.md) | The named materials TidalPy ships (simplified, rocky, icy, and giant-planet), loading them, overriding them, and adding new ones. |

```{toctree}
:maxdepth: 1

Equation-of-State and Shear-Modulus Laws <material_eos.md>
Phases and Materials <materials.md>
MatPack <matpack.md>
```

## Where Materials are Used

A layer holds one material: `Layer(..., material=...)` or `layer.material = ...` takes a `Material`, a MatPack name, or a material table. A world built from a TOML file gets each layer's material from its `material` key, a MatPack name or a `[layers.<name>.material]` table, and a layer that names none takes `[layers] material` of `TidalPy_Configs.toml` (`simple_rock`). See the [TOML schema](../Structures/config/toml_schema.md).

`BaseWorld.solve_eos()` integrates the planet's radial structure from the center to the surface, reading each layer's material at the local pressure, temperature, and radius with that layer's physics switches. The integration runs over radius with pressure as a state variable, so a pressure-dependent density law is evaluated at each step with the current pressure. The solver's outer loop adjusts the central pressure until the integrated surface pressure matches the requested boundary value. Where a layer can change state (`use_melting` on, `state = "auto"`, and a material that melts), the solve splits it into solid and liquid zones at the radius where the material's post-melt shear modulus crosses the liquid threshold, and the radial solver treats each zone as a layer. A thermal solve (`solve_temperature=True`) reads the expansivity, heat capacity, and conductivity of the same materials along the temperature profile its [cooling models](../Cooling/index.md) build. The solve, its zones, and its results are documented with the world class. See [Worlds](../Structures/worlds/worlds.md).

## Examples

`Demos/Physics/10_thermal_eos.ipynb` compares a constant-density interior with a compressible one and sweeps a melting mantle through its melting range, `Demos/Physics/15_thermal_interior.ipynb` solves a world's temperature profile, and `Benchmarks/EOS/EOS_vs_BurnMan.ipynb` checks the Birch-Murnaghan solve against BurnMan.

## References

This is not a comprehensive list but should give you a good starting point. Each MatPack file lists the sources of its own values.

- Birch, F. (1947). Finite elastic strain of cubic crystals. *Physical Review*, 71(11), 809-824. The finite-strain equation of state.
- Vinet, P., Ferrante, J., Rose, J. H., and Smith, J. R. (1987). Compressibility of solids. *Journal of Geophysical Research*, 92(B9), 9319-9325. The universal equation of state.
- Seager, S., Kuchner, M., Hier-Majumder, C. A., and Militzer, B. (2007). Mass-radius relationships for solid exoplanets. *The Astrophysical Journal*, 669, 1279-1297. The modified polytrope.
- Anderson, O. L. (1995). *Equations of State of Solids for Geophysics and Ceramic Science*. Oxford University Press. Thermal pressure and the adiabatic bulk modulus.
- Dziewonski, A. M., and Anderson, D. L. (1981). Preliminary reference Earth model. *Physics of the Earth and Planetary Interiors*, 25(4), 297-356. The tabulated profile shipped with TidalPy.
