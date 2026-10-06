# Rheology (`Rheology`)

_Updated: 2026-10-06_

`TidalPy.Rheology` utilizes a planet's material static shear (or bulk) modulus and viscosity to determine how it will respond to tidal (or loading) forcing. The result is a complex modulus $\mu^*(\omega)$ whose real part describes the elastic energy stored in the material while the imaginary portion quantifies the energy lost as frictional heat.

"Static" here means purely real, it does not mean constant. These are free to vary with radius, temperature, pressure, and time; what makes them static is that they carry no phase lag of their own.

| Page | Covers |
|---|---|
| [Rheology Models](rheology_models.md) | The constitutive models, their parameters and formulas, the factory, vectorized evaluation, serialization, and how to add a model. |

```{toctree}
:maxdepth: 1

Rheology Models <rheology_models.md>
```

## Where a Rheology is Used

A rheology model is the bridge between a material description and a response calculation, and it usually appears in three places.

The first is a layer. Every layer uses up to two rheology models, one for the shear response and one for the bulk response. A material's phases can carry default rheologies (`Phase(shear_rheology=..., bulk_rheology=...)`, as every MatPack solid does), and a layer's own `shear_rheology` and `bulk_rheology` override them. The layer's `calc_complex_shear_modulus(frequency)` then returns $\mu^*(\omega)$ instead of the purely real static modulus. A layer with neither responds elastically. See [Layer](../Structures/layers/layer.md) and [Phases and Materials](../Material/materials.md).

The second is the radial solver. Solving for Love numbers requires a complex modulus at every radius, which is what these models supply. A world with elastic layers returns a real $k_2$ and no dissipation; the same world with a Maxwell or Andrade rheology (for example) returns a complex $k_2$ whose imaginary part sets the heating. See [Calculating Love Numbers](../RadialSolver/calculating_love_numbers.md).

The third is the tidal solve itself, which evaluates the rheology once per forcing frequency in the Kaula mode list. See [Global Tides](../Tides/global_tides.md).

## Rheology and Viscosity

In the literature, the term _rheology_ generally includes how viscosity changes with temperature, pressure, and melt fraction. TidalPy separates the two: a [`Viscosity`](../Viscosity/viscosity_models.md) model maps temperature and pressure onto a viscosity, the material's [melt weakening](../PartialMelt/partial_melt_models.md#melt-weakening) lowers it where melt is present, and a rheology model maps that viscosity and the static modulus onto the complex modulus. A phase with no viscosity law reports a NaN viscosity, and a viscoelastic rheology then returns a NaN modulus.

## Examples

`Demos/Physics/P03_rheology_io.ipynb` maps the tidal heating of a homogeneous Io across shear modulus and viscosity for the Maxwell, Burgers, Andrade, and Sundberg-Cooper rheologies, and plots each rheology's heating against eccentricity. `Demos/Physics/P02_tidal_basics.ipynb` uses a fixed-Q tide model, which needs no rheology, for comparison: it computes Io's heating and sweeps it against eccentricity and $Q$.

## References

- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Efroimsky, M. (2012). Tidal dissipation compared to seismic dissipation: In small bodies, Earths, and super-Earths. *The Astrophysical Journal*, 746(2), 150.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models. *The Astrophysical Journal*, 857(2), 98.
