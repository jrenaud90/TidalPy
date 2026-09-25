# Rheology (`rheology_x`)

_Updated: 2026-09-16_

`TidalPy.rheology_x` utilizes a planet's material static shear (or bulk) modulus and viscosity to determine how it will respond to tidal (or loading) forcing. The result is a complex modulus $\mu^*(\omega)$ whose real part describes the elastic energy stored in the material and whose imaginary part quantifies the energy lost as frictional heat. Every dissipation quantity TidalPy produces, from a Love number to a heating rate to an orbital decay timescale, depends on that imaginary part.

"Static" here does not mean constant. It means purely real: the modulus and viscosity a material would show under a steady load, before any frequency dependence is applied. Those values are free to vary with radius, temperature, pressure, and time; what makes them static is that they carry no phase lag of their own.

| Page | Covers |
|---|---|
| [Rheology Models](rheology_models.md) | The seven constitutive models, their parameters and formulas, the factory, vectorized evaluation, serialization, and how to add a model. |

```{toctree}
:maxdepth: 1

Rheology Models <rheology_models.md>
```

## Where a Rheology is Used

A rheology model is the bridge between a material description and a response calculation, and it usually appears in three places.

The first is a layer. `PhysicsLayer` and its subclasses hold up to two rheology models, one for the shear response and one for the bulk response, attached with `set_shear_rheology` and `set_bulk_rheology`. The layer's `calc_complex_shear_modulus(frequency)` then returns $\mu^*(\omega)$ instead of the purely real static modulus. See [PhysicsLayer](../structures_x/layers/physics_layer.md).

The second is the radial solver. Solving for Love numbers requires a complex modulus at every radius, which is what these models supply. A world with elastic layers returns a real $k_2$ and no dissipation; the same world with a Maxwell or Andrade rheology (for example) returns a complex $k_2$ whose imaginary part sets the heating. See [Calculating Love Numbers](../RadialSolver_x/calculating_love_numbers.md).

The third is the tidal solve itself, which evaluates the rheology once per forcing frequency in the Kaula mode list. Because a single tidal solve can touch hundreds of modes, and each mode can need the modulus at every radial slice, the models are written in C++ and exposed with vectorized entry points. See [Global Tides](../Tides_x/global_tides.md).

## Rheology and Viscosity

In the literature, _rheology_ generally includes how viscosity changes with temperature, pressure, and melt fraction. TidalPy separates the two: a [`viscosity_x`](../viscosity_x/viscosity_models.md) model maps temperature and pressure onto a viscosity, and a rheology model maps that viscosity and the static modulus onto the complex modulus. A layer with a rheology but no viscosity model falls back to its static viscosity.

## Examples

`Demos_x/Physics/06_rheology_io.ipynb` builds models, sweeps them across frequency, and saves and reloads them. `Demos_x/Physics/05_tidal_basics.ipynb` shows the same models driving a tidal solve.

## References

- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Efroimsky, M. (2012). Tidal dissipation compared to seismic dissipation: In small bodies, Earths, and super-Earths. *The Astrophysical Journal*, 746(2), 150.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models. *The Astrophysical Journal*, 857(2), 98.
