"""
Tests for TidalPy.radiogenics_x.radiogenics.

Covers construction, model names and aliases, radiogenic-heating physics
(cross-checked against independent reference formulas and the legacy
``TidalPy/radiogenics/radiogenic_models.py``), the ``make_radiogenics`` factory,
vectorization, config dicts, binary round-trip, TOML config save, and the
isinstance chain.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import math
import os
import tempfile

import numpy as np
import pytest


# =====================================================================================================================
# Imports
# =====================================================================================================================
def _import_radiogenics():
    try:
        from TidalPy import radiogenics_x as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.radiogenics_x not compiled — run uv pip install first."
        )


# Representative MKS inputs. Times/half-lives share units (seconds here); the
# physics is unit-agnostic provided the inputs are consistent.
_MASS = 1.0e22          # radiogenic mass [kg]
_REF  = 0.0             # reference time  [s]

# A small synthetic isotope set (values chosen for an easy independent check).
_HPR    = [9.48e-5, 2.69e-5]   # heat production [W/kg]
_HALF   = [4.47e17, 1.40e18]   # half lives      [s]
_FRAC   = [0.9928, 0.9998]     # isotope mass fractions
_CONC   = [0.012e-6, 0.04e-6]  # element concentrations


# =====================================================================================================================
# Reference heating formulas (independent re-implementation of the legacy math).
# =====================================================================================================================
def _ref_isotope(time, mass, hpr, half, frac, conc, ref):
    specific = 0.0
    for h, hl, mf, c in zip(hpr, half, frac, conc):
        gamma = math.log(0.5) / hl
        specific += mf * c * h * math.exp(gamma * (time - ref))
    return specific * mass


def _ref_fixed(time, mass, rate, half, ref):
    if half <= 0.0:
        return mass * rate
    gamma = math.log(0.5) / half
    return mass * rate * math.exp(gamma * (time - ref))


# =====================================================================================================================
# Construction & names
# =====================================================================================================================
def test_model_names():
    """Each model reports its canonical lower-case model name."""
    mod = _import_radiogenics()
    assert mod.OffRadiogenics().model_name == "off"
    assert mod.IsotopeRadiogenics().model_name == "isotope"
    assert mod.FixedRadiogenics().model_name == "fixed"


def test_abstract_base_not_instantiable():
    """RadiogenicsBase is abstract and cannot be constructed directly."""
    mod = _import_radiogenics()
    with pytest.raises(TypeError):
        mod.RadiogenicsBase()


# =====================================================================================================================
# Heating physics (cross-check vs reference formulas)
# =====================================================================================================================
def test_off_heating_is_zero():
    """The Off model produces zero heating for any time and mass."""
    mod = _import_radiogenics()
    o = mod.OffRadiogenics()
    assert o.calc_heating(0.0, _MASS) == 0.0
    assert o.calc_heating(1.0e17, 5.0e21) == 0.0


def test_isotope_heating_matches_reference():
    """Isotope heating matches the independent reference at several times."""
    mod = _import_radiogenics()
    iso = mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC, _REF)
    for t in (0.0, 1.0e17, 5.0e17):
        got = iso.calc_heating(t, _MASS)
        expected = _ref_isotope(t, _MASS, _HPR, _HALF, _FRAC, _CONC, _REF)
        assert got == pytest.approx(expected, rel=1e-12)


def test_isotope_matches_legacy_model():
    """Isotope heating matches the legacy radiogenic_models.isotope function."""
    mod = _import_radiogenics()
    from TidalPy.radiogenics import radiogenic_models as legacy
    iso = mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC, _REF)
    t = 3.0e17
    got = iso.calc_heating(t, _MASS)
    # Legacy signature: isotope(time, mass, massfracs, concentrations, halflives, hpr, ref_time)
    expected = legacy.isotope(
        np.array([t]), _MASS,
        tuple(_FRAC), tuple(_CONC), tuple(_HALF), tuple(_HPR), _REF)[0]
    assert got == pytest.approx(expected, rel=1e-12)


def test_isotope_decays_over_time():
    """Isotope heating decreases monotonically for times after the reference."""
    mod = _import_radiogenics()
    iso = mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC, _REF)
    h0 = iso.calc_heating(0.0, _MASS)
    h1 = iso.calc_heating(5.0e17, _MASS)
    h2 = iso.calc_heating(1.0e18, _MASS)
    assert h0 > h1 > h2 > 0.0


def test_fixed_no_decay_is_constant():
    """The Fixed model with no decay scales linearly with mass and ignores time."""
    mod = _import_radiogenics()
    f = mod.FixedRadiogenics(1.0e-11, 0.0, 0.0)
    assert f.calc_heating(0.0, _MASS) == pytest.approx(1.0e-11 * _MASS)
    assert f.calc_heating(1.0e18, _MASS) == pytest.approx(1.0e-11 * _MASS)
    assert f.calc_heating(0.0, 2.0 * _MASS) == pytest.approx(2.0e-11 * _MASS)


def test_fixed_with_decay_matches_reference():
    """The Fixed model with decay matches the independent reference."""
    mod = _import_radiogenics()
    f = mod.FixedRadiogenics(1.0e-11, 4.47e17, 0.0)
    for t in (0.0, 4.47e17, 1.0e18):
        got = f.calc_heating(t, _MASS)
        expected = _ref_fixed(t, _MASS, 1.0e-11, 4.47e17, 0.0)
        assert got == pytest.approx(expected, rel=1e-12)


def test_fixed_half_life_halves_heating():
    """At one half life past the reference, decayed heating is half the initial."""
    mod = _import_radiogenics()
    half = 4.47e17
    f = mod.FixedRadiogenics(1.0e-11, half, 0.0)
    h0 = f.calc_heating(0.0, _MASS)
    h_half = f.calc_heating(half, _MASS)
    assert h_half == pytest.approx(0.5 * h0, rel=1e-12)


# =====================================================================================================================
# Parameter getters
# =====================================================================================================================
def test_isotope_parameters():
    mod = _import_radiogenics()
    iso = mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC, 1.0e17)
    assert iso.num_isotopes == 2
    assert iso.ref_time == pytest.approx(1.0e17)
    assert list(iso.heat_production) == pytest.approx(_HPR)
    assert list(iso.half_lives) == pytest.approx(_HALF)
    assert list(iso.mass_fracs) == pytest.approx(_FRAC)
    assert list(iso.concentrations) == pytest.approx(_CONC)
    # Names are auto-generated when not supplied.
    assert list(iso.isotope_names) == ["isotope_0", "isotope_1"]


def test_isotope_explicit_names():
    """Explicit isotope names are preserved."""
    mod = _import_radiogenics()
    iso = mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC, names=["U238", "Th232"])
    assert list(iso.isotope_names) == ["U238", "Th232"]


def test_isotope_array_length_mismatch_raises():
    """Mismatched isotope array lengths raise ValueError."""
    mod = _import_radiogenics()
    with pytest.raises(ValueError):
        mod.IsotopeRadiogenics([1.0e-4, 2.0e-4], [1.0e17], [0.9], [1.0e-8])


def test_fixed_parameters():
    mod = _import_radiogenics()
    f = mod.FixedRadiogenics(2.5e-11, 1.0e18, 3.0e17)
    assert f.fixed_heat_production == pytest.approx(2.5e-11)
    assert f.average_half_life == pytest.approx(1.0e18)
    assert f.ref_time == pytest.approx(3.0e17)


# =====================================================================================================================
# Factory
# =====================================================================================================================
@pytest.mark.parametrize("alias,canonical", [
    ("off",       "off"),
    ("none",      "off"),
    ("OFF",       "off"),
    ("isotope",   "isotope"),
    ("isotopes",  "isotope"),
    ("Isotope",   "isotope"),
    ("fixed",     "fixed"),
    ("constant",  "fixed"),
    ("FIXED",     "fixed"),
])
def test_make_radiogenics_aliases(alias, canonical):
    """The factory resolves aliases (case-insensitive) to canonical models."""
    mod = _import_radiogenics()
    inst = mod.make_radiogenics(alias)
    assert inst.model_name == canonical


def test_make_radiogenics_unknown_raises():
    """An unknown model name raises ValueError (from the C++ enum mapping)."""
    mod = _import_radiogenics()
    with pytest.raises(ValueError):
        mod.make_radiogenics("not_a_real_model")


def test_make_radiogenics_fixed_config():
    """The factory forwards fixed-model config parameters."""
    mod = _import_radiogenics()
    f = mod.make_radiogenics("fixed", {
        "fixed_heat_production_w_kg": 3.0e-11,
        "average_half_life_s": 1.0e18,
        "ref_time_s": 2.0e17,
    })
    assert f.fixed_heat_production == pytest.approx(3.0e-11)
    assert f.average_half_life == pytest.approx(1.0e18)
    assert f.ref_time == pytest.approx(2.0e17)


def test_make_radiogenics_isotope_explicit():
    """The factory builds an isotope model from explicit MKS arrays."""
    mod = _import_radiogenics()
    iso = mod.make_radiogenics("isotope", {
        "heat_production_w_kg": _HPR,
        "half_lives_s": _HALF,
        "mass_fracs": _FRAC,
        "concentrations": _CONC,
        "ref_time_s": 0.0,
    })
    assert iso.num_isotopes == 2
    expected = _ref_isotope(0.0, _MASS, _HPR, _HALF, _FRAC, _CONC, 0.0)
    assert iso.calc_heating(0.0, _MASS) == pytest.approx(expected, rel=1e-12)


def test_make_radiogenics_isotope_named_dataset():
    """The factory looks up a named isotope dataset and converts Myr -> s."""
    mod = _import_radiogenics()
    iso = mod.make_radiogenics("isotope", {"isotopes": "modern_day_chondritic"})
    assert iso.num_isotopes == 4
    # Dataset ref_time is 4600 Myr; converted to seconds it must be > 1e17 s.
    assert iso.ref_time > 1.0e17
    # Half lives are stored in seconds (Myr-scale -> ~1e16-1e17 s).
    assert all(hl > 1.0e16 for hl in iso.half_lives)


def test_make_radiogenics_inline_dataset():
    """An inline isotope dataset dict (Myr units) is converted to seconds."""
    mod = _import_radiogenics()
    inline = {
        "ref_time": 4600.0,
        "U238": {"hpr": 9.48e-5, "half_life": 4470.0,
                 "iso_mass_fraction": 0.9928, "element_concentration": 0.012e-6},
    }
    iso = mod.make_radiogenics("isotope", {"isotopes": inline})
    assert iso.num_isotopes == 1
    assert iso.ref_time == pytest.approx(4600.0 * 1.0e6 * 365.25 * 24.0 * 3600.0)
    assert iso.half_lives[0] == pytest.approx(4470.0 * 1.0e6 * 365.25 * 24.0 * 3600.0)


# =====================================================================================================================
# Built-in literature isotope datasets
# =====================================================================================================================
def test_available_isotope_datasets():
    """The three built-in datasets are advertised."""
    mod = _import_radiogenics()
    names = mod.available_isotope_datasets()
    assert set(names) == {"modern_day_chondritic", "llri_and_slri", "bulk_silicate_earth"}


_MYR_S = 1.0e6 * 365.25 * 24.0 * 3600.0


@pytest.mark.parametrize("name,n_isotopes,ref_time_s", [
    ("modern_day_chondritic", 4, 4600.0 * _MYR_S),
    ("llri_and_slri", 7, 0.0),
    ("bulk_silicate_earth", 4, 4600.0 * _MYR_S),
])
def test_isotope_dataset_contents(name, n_isotopes, ref_time_s):
    """Each built-in dataset returns an MKS dict with the expected isotope count.

    The present-epoch datasets quote concentrations 4600 Myr after formation; the
    LLRI+SLRI dataset quotes formation (CAI) abundances, so its reference time is 0.
    """
    mod = _import_radiogenics()
    ds = mod.isotope_dataset(name)
    assert set(ds) == {
        "heat_production_w_kg", "half_lives_s", "mass_fracs",
        "concentrations", "isotope_names", "ref_time_s"}
    assert len(ds["isotope_names"]) == n_isotopes
    assert ds["ref_time_s"] == pytest.approx(ref_time_s, abs=1.0)
    assert all(hl > 0.0 for hl in ds["half_lives_s"])


def test_llri_heating_finite_at_formation():
    """The LLRI+SLRI dataset gives finite, decaying heating from formation onward.

    The short-lived isotopes (Al26 with a 0.72 Myr half life, Fe60, Mn53) contribute
    roughly half the formation-epoch heating and are gone within a few tens of Myr;
    the long-lived isotopes persist, so the total decays substantially (the SLRI
    share dies off) but not by orders of magnitude.
    """
    import math
    mod = _import_radiogenics()
    model = mod.IsotopeRadiogenics.from_dataset("llri_and_slri")
    heating_formation = model.calc_heating(0.0, _MASS)
    heating_10myr = model.calc_heating(10.0 * _MYR_S, _MASS)
    heating_100myr = model.calc_heating(100.0 * _MYR_S, _MASS)
    assert math.isfinite(heating_formation)
    assert heating_formation > 0.0
    # Monotonic decay, with the SLRI share (about half the formation total) gone by 100 Myr.
    assert heating_10myr < heating_formation
    assert heating_100myr < heating_10myr
    assert heating_100myr < 0.6 * heating_formation


def test_isotope_dataset_unknown_raises():
    """An unknown built-in dataset name raises ValueError."""
    mod = _import_radiogenics()
    with pytest.raises(ValueError):
        mod.isotope_dataset("not_a_dataset")


def test_llri_includes_short_lived_isotopes():
    """The LLRI+SLRI dataset includes the short-lived isotopes Al26, Fe60, Mn53."""
    mod = _import_radiogenics()
    ds = mod.isotope_dataset("llri_and_slri")
    assert {"Al26", "Fe60", "Mn53"}.issubset(set(ds["isotope_names"]))


def test_from_dataset_matches_factory():
    """IsotopeRadiogenics.from_dataset matches make_radiogenics for the same name."""
    mod = _import_radiogenics()
    a = mod.IsotopeRadiogenics.from_dataset("modern_day_chondritic")
    b = mod.make_radiogenics("isotope", {"isotopes": "modern_day_chondritic"})
    assert isinstance(a, mod.IsotopeRadiogenics)
    assert list(a.isotope_names) == list(b.isotope_names)
    assert a.calc_heating(a.ref_time, _MASS) == pytest.approx(
        b.calc_heating(b.ref_time, _MASS), rel=1e-12)


def test_dataset_heating_cross_check():
    """Built-in dataset heating matches an independent reference (via dataset dict)."""
    mod = _import_radiogenics()
    ds = mod.isotope_dataset("modern_day_chondritic")
    iso = mod.IsotopeRadiogenics.from_dataset("modern_day_chondritic")
    t = ds["ref_time_s"]
    expected = _ref_isotope(
        t, _MASS, ds["heat_production_w_kg"], ds["half_lives_s"],
        ds["mass_fracs"], ds["concentrations"], ds["ref_time_s"])
    assert iso.calc_heating(t, _MASS) == pytest.approx(expected, rel=1e-12)


def test_make_radiogenics_returns_rich_subclass():
    """The factory returns the correct concrete Python subclass (not the base)."""
    mod = _import_radiogenics()
    assert type(mod.make_radiogenics("off")).__name__ == "OffRadiogenics"
    assert type(mod.make_radiogenics("isotope")).__name__ == "IsotopeRadiogenics"
    assert type(mod.make_radiogenics("fixed")).__name__ == "FixedRadiogenics"


def test_make_radiogenics_adopted_object_is_usable():
    """An object built by the C++ enum factory is fully usable and round-trips."""
    mod = _import_radiogenics()
    f = mod.make_radiogenics("fixed", {
        "fixed_heat_production_w_kg": 1.5e-11, "average_half_life_s": 5.0e17})
    assert f.fixed_heat_production == pytest.approx(1.5e-11)
    expected = _ref_fixed(1.0e17, _MASS, 1.5e-11, 5.0e17, 0.0)
    assert f.calc_heating(1.0e17, _MASS) == pytest.approx(expected, rel=1e-12)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "adopted.tpyb")
        f.save_binary(path)
        restored = mod.FixedRadiogenics()
        restored.load_binary(path)
    assert restored.get_config_dict() == f.get_config_dict()


# =====================================================================================================================
# Config dict
# =====================================================================================================================
def test_config_dict_keys():
    """Each model's config dict carries the expected keys."""
    mod = _import_radiogenics()
    assert set(mod.OffRadiogenics().get_config_dict()) == {"model_name"}
    assert set(mod.FixedRadiogenics().get_config_dict()) == {
        "model_name", "fixed_heat_production_w_kg", "average_half_life_s", "ref_time_s"}
    assert set(mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC).get_config_dict()) == {
        "model_name", "heat_production_w_kg", "half_lives_s",
        "mass_fracs", "concentrations", "isotope_names", "ref_time_s"}


def test_save_config_writes_toml():
    """save_config writes a TOML file that round-trips through toml.load."""
    import toml
    mod = _import_radiogenics()
    f = mod.FixedRadiogenics(2.5e-11, 1.0e18, 3.0e17)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "radio.toml")
        f.save_config(path)
        loaded = toml.load(path)
    assert loaded["model_name"] == "fixed"
    assert loaded["fixed_heat_production_w_kg"] == pytest.approx(2.5e-11)
    assert loaded["average_half_life_s"] == pytest.approx(1.0e18)


# =====================================================================================================================
# Vectorized methods and direct convenience functions
# =====================================================================================================================
def test_convenience_scalar_matches_class():
    """The lower-case convenience function matches the class for scalar input."""
    mod = _import_radiogenics()
    got = mod.fixed(1.0e17, _MASS, 1.0e-11, 5.0e17, 0.0)
    inst = mod.FixedRadiogenics(1.0e-11, 5.0e17, 0.0)
    expected = inst.calc_heating(1.0e17, _MASS)
    assert isinstance(got, float)
    assert got == pytest.approx(expected)


def test_convenience_vectorize_time():
    """Array time at constant mass returns a matching ndarray."""
    mod = _import_radiogenics()
    inst = mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC, 0.0)
    times = np.array([0.0, 1.0e17, 5.0e17, 1.0e18])
    got = mod.isotope(times, _MASS, _HPR, _HALF, _FRAC, _CONC, ref_time_s=0.0)
    assert isinstance(got, np.ndarray)
    assert got.shape == (4,)
    assert got.dtype == np.float64
    for i in range(4):
        assert got[i] == pytest.approx(inst.calc_heating(times[i], _MASS))


def test_convenience_vectorize_mass():
    """Array mass at constant time returns a matching ndarray."""
    mod = _import_radiogenics()
    inst = mod.FixedRadiogenics(1.0e-11, 0.0, 0.0)
    masses = np.array([1.0e21, 5.0e21, 1.0e22])
    got = mod.fixed(0.0, masses, 1.0e-11, 0.0, 0.0)
    assert got.shape == (3,)
    for i in range(3):
        assert got[i] == pytest.approx(inst.calc_heating(0.0, masses[i]))


def test_convenience_vectorize_all_and_broadcast():
    """All-array and mixed (broadcast) inputs return the broadcast shape."""
    mod = _import_radiogenics()
    times  = np.array([0.0, 1.0e17, 5.0e17])
    masses = np.array([1.0e21, 5.0e21, 1.0e22])
    got_all = mod.fixed(times, masses, 1.0e-11, 5.0e17, 0.0)
    assert got_all.shape == (3,)
    # Mixed: time array, mass scalar -> broadcast.
    got_mixed = mod.fixed(times, _MASS, 1.0e-11, 5.0e17, 0.0)
    assert got_mixed.shape == (3,)
    inst = mod.FixedRadiogenics(1.0e-11, 5.0e17, 0.0)
    for i in range(3):
        assert got_mixed[i] == pytest.approx(inst.calc_heating(times[i], _MASS))


def test_convenience_preserves_2d_shape():
    """A 2-D input shape is preserved in the output array."""
    mod = _import_radiogenics()
    times = np.array([[0.0, 1.0e17], [2.0e17, 3.0e17]])
    got = mod.fixed(times, _MASS, 1.0e-11, 5.0e17, 0.0)
    assert got.shape == (2, 2)


def test_class_vectorize_methods_match_scalar():
    """The class-level vectorized methods match element-wise scalar calls."""
    mod = _import_radiogenics()
    iso = mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC, 0.0)
    times  = np.array([0.0, 1.0e17, 5.0e17])
    masses = np.array([1.0e21, 5.0e21, 1.0e22])

    by_time = iso.calc_heating_vectorize_time(times, _MASS)
    by_mass = iso.calc_heating_vectorize_mass(0.0, masses)
    by_all  = iso.calc_heating_vectorize_all(times, masses)
    for i in range(3):
        assert by_time[i] == pytest.approx(iso.calc_heating(times[i], _MASS))
        assert by_mass[i] == pytest.approx(iso.calc_heating(0.0, masses[i]))
        assert by_all[i]  == pytest.approx(iso.calc_heating(times[i], masses[i]))


def test_class_vectorize_size_mismatch_raises():
    """Mismatched input lengths raise ValueError from the C++ layer."""
    mod = _import_radiogenics()
    f = mod.FixedRadiogenics(1.0e-11, 0.0, 0.0)
    with pytest.raises(ValueError):
        f.calc_heating_vectorize_all(np.ones(3), np.ones(2))


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
def test_binary_round_trip_off():
    mod = _import_radiogenics()
    original = mod.OffRadiogenics()
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "off.tpyb")
        original.save_binary(path)
        restored = mod.OffRadiogenics()
        restored.load_binary(path)
    assert restored.model_name == "off"
    assert restored.get_config_dict() == original.get_config_dict()


def test_binary_round_trip_fixed():
    mod = _import_radiogenics()
    original = mod.FixedRadiogenics(2.5e-11, 1.0e18, 3.0e17)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "fixed.tpyb")
        original.save_binary(path)
        restored = mod.FixedRadiogenics()
        restored.load_binary(path)
    assert restored.model_name == "fixed"
    assert restored.get_config_dict() == original.get_config_dict()


def test_binary_round_trip_isotope():
    """Isotope binary round-trip preserves the variable-length isotope list."""
    mod = _import_radiogenics()
    original = mod.IsotopeRadiogenics(_HPR, _HALF, _FRAC, _CONC, 1.0e17, names=["U238", "Th232"])
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "isotope.tpyb")
        original.save_binary(path)
        restored = mod.IsotopeRadiogenics()
        restored.load_binary(path)
    assert restored.model_name == "isotope"
    assert restored.num_isotopes == 2
    assert list(restored.isotope_names) == ["U238", "Th232"]
    assert list(restored.half_lives) == pytest.approx(_HALF)
    assert restored.calc_heating(0.0, _MASS) == pytest.approx(
        original.calc_heating(0.0, _MASS), rel=1e-12)
    assert restored.get_config_dict() == original.get_config_dict()


def test_binary_round_trip_dataset():
    """A dataset-built isotope model round-trips through binary."""
    mod = _import_radiogenics()
    original = mod.IsotopeRadiogenics.from_dataset("llri_and_slri")
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "llri.tpyb")
        original.save_binary(path)
        restored = mod.IsotopeRadiogenics()
        restored.load_binary(path)
    assert restored.num_isotopes == 7
    assert list(restored.isotope_names) == list(original.isotope_names)
    assert restored.get_config_dict() == original.get_config_dict()


def test_load_binary_missing_file_raises():
    """Loading a non-existent binary raises FileNotFoundError."""
    mod = _import_radiogenics()
    with pytest.raises(FileNotFoundError):
        mod.FixedRadiogenics().load_binary("does_not_exist_12345.tpyb")


# =====================================================================================================================
# isinstance chain
# =====================================================================================================================
@pytest.mark.parametrize("cls", ["OffRadiogenics", "IsotopeRadiogenics", "FixedRadiogenics"])
def test_isinstance_chain(cls):
    """Every model is a RadiogenicsBase, PhysicsBase and TidalPyBaseClass."""
    mod = _import_radiogenics()
    from TidalPy.Utilities_x.classes_x.classes import PhysicsBase, TidalPyBaseClass
    inst = getattr(mod, cls)()
    assert isinstance(inst, mod.RadiogenicsBase)
    assert isinstance(inst, PhysicsBase)
    assert isinstance(inst, TidalPyBaseClass)
