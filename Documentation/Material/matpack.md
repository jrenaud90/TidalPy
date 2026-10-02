# MatPack (`Material.matpack`)

_Updated: 2026-10-02_

`TidalPy.Material.matpack` contains the named materials TidalPy ships and the functions that load them. Each MatPack material is a complete `Material`: a solid phase, a liquid phase, or both. Each phase has an equation of state and a thermal conductivity and heat capacity; a solid phase adds a shear-modulus law, a viscosity law, and a default tidal rheology, and most liquid phases a viscosity. A material that melts adds its melting curves, melt weakening, and latent heat. A material is loaded by name, optionally with overrides, and returns its state at a pressure [Pa], temperature [K], and radius [m] through `Material.calc_state` (see [Phases and Materials](materials.md)).

The values come from the literature, with each file listing its sources and the confidence of the less certain values in its comments. Several values are derived (fits through published data, mineral-physics assemblages computed with BurnMan) or are unsourced estimates; the file says which.

## Materials

**Simplified**

Constant density, moduli, and viscosity, and no melting. They are quick to evaluate and a good start for a new world.

| Material | Description |
|---|---|
| `simple_rock` | Uniform silicate rock; 3300 kg m$^{-3}$, $\mu$ = 60 GPa, $\eta$ = 10$^{21}$ Pa s |
| `simple_ice` | Uniform water ice; 920 kg m$^{-3}$, $\mu$ = 3.5 GPa, $\eta$ = 10$^{14}$ Pa s |
| `simple_iron_core` | Uniform solid iron-rich core; 8000 kg m$^{-3}$, $\mu$ = 80 GPa, $\eta$ = 10$^{20}$ Pa s |
| `simple_liquid_iron` | Uniform liquid iron-alloy core; always liquid |
| `simple_water` | Uniform liquid water; always liquid |
| `simple_gas` | Uniform gas-giant envelope placeholder; always fluid |

**Rocky**

| Material | Description |
|---|---|
| `peridotite` | Upper-mantle peridotite: Birch-Murnaghan, dry olivine diffusion and dislocation creep, Monteux et al. (2016) melting curves, and a silicate melt |
| `lower_mantle` | Bridgmanite-dominated lower mantle with the peridotite melting curves |
| `olivine` | Mantle olivine (Fo90) with a forsterite melting curve |
| `basalt` | Mafic crust (dense gabbro) with diabase creep |
| `felsic_crust` | Granitic continental crust with wet-quartzite creep |
| `chondrite` | Undifferentiated chondritic rock at zero porosity |
| `serpentinite` | Antigorite serpentinite; dehydration is not represented |
| `iron` | Solid iron with the Anzellini et al. (2013) melting curve, melting into liquid iron |
| `liquid_iron` | Liquid iron (Anderson and Ahrens 1994); always liquid |
| `iron_sulfide` | Fe-FeS core alloy at 16 wt% S with a eutectic solidus (low confidence) |

**Icy**

| Material | Description |
|---|---|
| `ice_ih` | Water ice Ih (IAPWS-06) with diffusion creep, melting into liquid water |
| `ice_iii`, `ice_v`, `ice_vi`, `ice_vii` | High-pressure ices, each with its IAPWS melting curve, melting into liquid water |
| `water` | Liquid water (IAPWS); always liquid |
| `brine` | Standard seawater (TEOS-10); always liquid |
| `ammonia_water` | Water ice with 10 wt% NH$_3$ and an ammonia-water melt (low confidence) |
| `methane_clathrate` | Structure I methane clathrate hydrate; dissociation is not represented |
| `nitrogen_ice` | Beta nitrogen ice, melting at 63.15 K |

**Giant planets**

| Material | Description |
|---|---|
| `h2_he_molecular` | Molecular H$_2$-He envelope on Jupiter's $n$ = 1 polytrope; always fluid |
| `h_he_metallic` | Metallic H-He envelope on the same polytrope; always fluid |
| `giant_core` | Rock core (Seager et al. 2007 MgSiO$_3$ modified polytrope); treated as fluid |

`available_materials()` lists the names, `available_materials("icy")` one category, and `material_info(name)` a material's description, category, and file.

### Behavior at the Limits

Each material is one phase assemblage at one composition:

- An ice polymorph is valid only in its own pressure band; a deep water layer is one layer per polymorph. Past the end of its band a melting curve can fall to zero (ice Ih's above 415 MPa), so the material reads as fully molten there when melting and pressure melting are on.
- Water's expansivity changes sign near 277 K at 1 bar; `water` uses a positive value, as in a pressurized ocean.
- Materials whose melt composition changes with temperature (`iron_sulfide`, `ammonia_water`) fix the bulk composition, and their melt is one composition.
- `serpentinite` and `methane_clathrate` decompose rather than melt, so they ship without melting.
- The two hydrogen-helium materials share one polytrope, so there is no density step between them.
- A material with a single melting temperature (the ices, `olivine`, `iron`, `nitrogen_ice`) melts as a step, and its latent heat does not enter the effective heat capacity: there is no melting range to spread it over. In a layer that can change state, the boundary between its solid and liquid zones carries the latent heat instead (the world's `layer_latent_capacity`).

## Python API

```python
from TidalPy.Material import load_material

# Upper-mantle peridotite
peridotite = load_material("peridotite")

# Its state at 3 GPa and 1600 K, melting on
state = peridotite.calc_state(
    3.0e9,
    1600.0,
    use_melting=True,
    use_pressure_melting=True)
print(state["density"], state["shear_viscosity"], state["solidus"])
```

A material is immutable. `with_parameters` and `replace` return changed copies.

### Overrides and Presets

Keyword overrides are merged over the material's table, one table at a time, so they name only what they change. A model table that names a different model than the material's replaces its table instead, since another model reads other keys, and `None` removes a slot. A material must keep a solid or a liquid phase.

```python
from TidalPy.Material import load_material

# Io's mantle after Renaud and Henning (2018): peridotite with constant melting temperatures
io_mantle = load_material(
    "peridotite",
    solid={"shear_rheology": {"model": "andrade", "alpha": 0.2}},
    melting={"solidus": {"model": "constant", "temperature_k": 1600.0},
             "liquidus": {"model": "constant", "temperature_k": 2000.0}})

# Peridotite that cannot melt
dry_rock = load_material(
    "peridotite",
    liquid=None,
    melting=None,
    latent_heat_j_kg=0.0)
```

The same overrides can be written as one table with a `preset` key, the form a TOML file holds. A `solid` or `liquid` table may name a preset too, and then starts from that material's phase in the same slot:

```python
from TidalPy.Material import load_material

# A freezing salty ocean: ice Ih over seawater, melting between the NaCl eutectic and the seawater freezing point
salty_ocean = load_material({
    "preset": "ice_ih",
    "liquid": {"preset": "brine"},
    "melting": {
        "solidus": {"model": "constant", "temperature_k": 251.9},
        "liquidus": {"model": "constant", "temperature_k": 271.2}}})
print(salty_ocean.calc_state(1.0e7, 260.0, use_melting=True)["melt_fraction"])
```

`material_config(source, overrides)` returns the resolved table without building the material, and `merge_material_tables(base, overrides)` applies the merge rules to any two tables.

### Convenience Functions

| Function | Returns |
|---|---|
| `load_material(source, **overrides)` | A `Material` from a name or a table |
| `material_config(source, overrides=None)` | The resolved material table |
| `available_materials(category=None)` | The sorted material names |
| `material_info(name)` | `name`, `description`, `category`, and the `path` it is read from |
| `install_matpack(force=False)` | Copies the packaged files into the data directory |

## Data Directory

The packaged files are copied into `<documents>/TidalPy/<version>/Materials` on first use (copy-if-absent), and a material is read from that copy, so an edit there changes the material everywhere it is named. A copy that differs from the packaged file is reported once per session, since it may be an edit or a copy left by an older install; `install_matpack(force=True)` replaces every copy, and `stale_matpack_copy = false` under `[warnings]` in `TidalPy_Configs.toml` silences the report. Without a writable data directory the packaged files are read directly. The WorldPack's bundled worlds work the same way.

## Adding a Material

A new file in the data directory is a new material, named by its file name. It holds the material table and three metadata keys, and may start from another material:

```toml
schema_version = "0.2.0"
description = "Io's mantle after Renaud and Henning (2018)."
category = "rocky"
preset = "peridotite"

[solid.shear_rheology]
model = "andrade"
alpha = 0.2

[melting.solidus]
model = "constant"
temperature_k = 1600.0

[melting.liquidus]
model = "constant"
temperature_k = 2000.0
```

A full material table has a `solid` table, a `liquid` table, or both, each with an `eos` table and optional `shear_modulus`, `shear_viscosity`, `bulk_viscosity`, `shear_rheology`, and `bulk_rheology` tables and the thermal parameters; a `melting` table with `solidus`, `liquidus`, and optional `weakening`, `bulk_modulus_mixing`, and `bulk_viscosity_mixing` tables; and `latent_heat_j_kg`. `Material.get_config_dict()` returns this form for any material. To add a material to TidalPy itself, add its file to `TidalPy/MatPack` with its references in the comments.

## References

Each MatPack file lists its own sources. The most used are:

- Stixrude, L., and Lithgow-Bertelloni, C. (2011). Thermodynamics of mantle minerals II. Phase equilibria. *Geophysical Journal International*, 184, 1180-1213.
- Monteux, J., Andrault, D., and Samuel, H. (2016). On the cooling of a deep terrestrial magma ocean. *Earth and Planetary Science Letters*, 448, 140-149.
- Hirth, G., and Kohlstedt, D. L. (2003). Rheology of the upper mantle and the mantle wedge: a view from the experimentalists. In *Inside the Subduction Factory*, Geophysical Monograph 138, AGU, 83-105.
- Feistel, R., and Wagner, W. (2006). A new equation of state for H2O ice Ih. *Journal of Physical and Chemical Reference Data*, 35, 1021-1047.
- Wagner, W., Riethmann, T., Feistel, R., and Harvey, A. H. (2011). New equations for the sublimation pressure and melting pressure of H2O ice Ih. *Journal of Physical and Chemical Reference Data*, 40, 043103.
- Journaux, B., et al. (2020). Holistic approach for studying planetary hydrospheres: Gibbs representation of ices thermodynamics, elasticity, and the water phase diagram to 2,300 MPa. *Journal of Geophysical Research: Planets*, 125, e2019JE006176.
- Goldsby, D. L., and Kohlstedt, D. L. (2001). Superplastic deformation of ice: experimental observations. *Journal of Geophysical Research*, 106, 11017-11030.
- Anderson, W. W., and Ahrens, T. J. (1994). An equation of state for liquid iron and implications for the Earth's core. *Journal of Geophysical Research*, 99, 4273-4284.
- Anzellini, S., Dewaele, A., Mezouar, M., Loubeyre, P., and Morard, G. (2013). Melting of iron at Earth's inner core boundary based on fast X-ray diffraction. *Science*, 340, 464-466.
- Seager, S., Kuchner, M., Hier-Majumder, C. A., and Militzer, B. (2007). Mass-radius relationships for solid exoplanets. *The Astrophysical Journal*, 669, 1279-1297.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857, 98.
