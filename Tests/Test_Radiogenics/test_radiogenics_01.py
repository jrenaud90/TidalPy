"""Radiogenic heating models (off, isotope, fixed): physics, datasets, factory, vectorization, and I/O.

Legacy values in ``frozen/test_radiogenics_01.npz`` come from the classic backend (0.8.0 snapshot 8b8e0b12).
"""

import math
from pathlib import Path

import numpy as np
import pytest

import TidalPy
from TidalPy import Radiogenics
from TidalPy.constants import seconds_per_myr
from TidalPy.Utilities.classes.classes import PhysicsBase, TidalPyBaseClass

_FROZEN_PATH = Path(__file__).parent / "frozen" / "test_radiogenics_01.npz"
with np.load(_FROZEN_PATH, allow_pickle=False) as _frozen_file:
    _LEGACY_REFERENCE = {key: _frozen_file[key] for key in _frozen_file.files}

_MASS = 1.0e22          # [kg]
_MYR_S = 1.0e6 * 365.25 * 24.0 * 3600.0

# A small synthetic isotope set: heat production [W/kg], half lives [s], mass fractions, concentrations.
_HPR = [9.48e-5, 2.69e-5]
_HALF = [4.47e17, 1.40e18]
_FRAC = [0.9928, 0.9998]
_CONC = [0.012e-6, 0.04e-6]
_ISOTOPES = (_HPR, _HALF, _FRAC, _CONC)

# Fixed-model parameters: (heat production [W/kg], half life [s], reference time [s]).
_STEADY_FIXED = (1.0e-11, 0.0, 0.0)
_DECAYING_FIXED = (1.0e-11, 5.0e17, 0.0)

_MODELS = [("OffRadiogenics", "off"), ("IsotopeRadiogenics", "isotope"), ("FixedRadiogenics", "fixed")]


def _isotope_model(ref_time=0.0, **kwargs):
    return Radiogenics.IsotopeRadiogenics(*_ISOTOPES, ref_time=ref_time, **kwargs)


def _reference_isotope_heating(time, mass, isotopes, ref_time):
    """Independent re-implementation of the legacy isotope heating."""
    specific = 0.0
    for heat_production, half_life, mass_frac, concentration in zip(*isotopes):
        gamma = math.log(0.5) / half_life
        specific += mass_frac * concentration * heat_production * math.exp(gamma * (time - ref_time))
    return specific * mass


def _reference_fixed_heating(time, mass, rate, half_life):
    return mass * rate * math.exp(math.log(0.5) / half_life * time)


def _binary_copy(model, blank, tmp_path):
    """Save a model to a binary file and load it into a blank model."""
    path = str(tmp_path / "radiogenics.tpyb")
    model.save_binary(path)
    blank.load_binary(path)
    return blank


# =====================================================================================================================
# Model identity
# =====================================================================================================================
@pytest.mark.parametrize("cls_name,model_name", _MODELS)
def test_model_name(cls_name, model_name):
    assert getattr(Radiogenics, cls_name)().model_name == model_name


@pytest.mark.parametrize("cls_name,model_name", _MODELS)
def test_factory_returns_rich_subclass(cls_name, model_name):
    assert type(Radiogenics.make_radiogenics(model_name)).__name__ == cls_name


@pytest.mark.parametrize("cls_name,model_name", _MODELS)
def test_isinstance_chain(cls_name, model_name):
    model = getattr(Radiogenics, cls_name)()
    assert isinstance(model, Radiogenics.RadiogenicsBase)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)


@pytest.mark.parametrize("call,error", [
    (lambda: Radiogenics.RadiogenicsBase(), TypeError),
    (lambda: Radiogenics.make_radiogenics("not_a_real_model"), ValueError),
    (lambda: Radiogenics.IsotopeRadiogenics([1.0e-4, 2.0e-4], [1.0e17], [0.9], [1.0e-8]), ValueError),
    (lambda: Radiogenics.isotope_dataset("not_a_dataset"), ValueError),
    (lambda: Radiogenics.FixedRadiogenics(*_STEADY_FIXED).calc_heating_vectorize_all(np.ones(3), np.ones(2)),
     ValueError),
    (lambda: Radiogenics.FixedRadiogenics().load_binary("does_not_exist_12345.tpyb"), FileNotFoundError),
], ids=["abstract_base", "unknown_model", "isotope_length_mismatch", "unknown_dataset", "vectorize_size_mismatch",
        "missing_binary"])
def test_invalid_use_raises(call, error):
    with pytest.raises(error):
        call()


# =====================================================================================================================
# Heating physics
# =====================================================================================================================
def test_off_heating_is_zero():
    off = Radiogenics.OffRadiogenics()
    assert off.calc_heating(0.0, _MASS) == 0.0
    assert off.calc_heating(1.0e17, 5.0e21) == 0.0


@pytest.mark.parametrize("time", [0.0, 1.0e17, 5.0e17])
def test_isotope_heating_matches_reference(time):
    expected = _reference_isotope_heating(time, _MASS, _ISOTOPES, 0.0)
    assert _isotope_model().calc_heating(time, _MASS) == pytest.approx(expected, rel=1e-12)


def test_isotope_matches_legacy_model():
    # Frozen from the legacy isotope(np.array([t]), _MASS, _FRAC, _CONC, _HALF, _HPR, 0.0)[0].
    expected = float(_LEGACY_REFERENCE["isotope_legacy_model"])
    assert _isotope_model().calc_heating(3.0e17, _MASS) == pytest.approx(expected, rel=1e-12)


def test_isotope_decays_over_time():
    isotope = _isotope_model()
    heating = [isotope.calc_heating(time, _MASS) for time in (0.0, 5.0e17, 1.0e18)]
    assert heating[0] > heating[1] > heating[2] > 0.0


def test_fixed_no_decay_is_constant():
    """With no decay the fixed model scales with mass and ignores time."""
    fixed = Radiogenics.FixedRadiogenics(*_STEADY_FIXED)
    assert fixed.calc_heating(0.0, _MASS) == pytest.approx(1.0e-11 * _MASS)
    assert fixed.calc_heating(1.0e18, _MASS) == pytest.approx(1.0e-11 * _MASS)
    assert fixed.calc_heating(0.0, 2.0 * _MASS) == pytest.approx(2.0e-11 * _MASS)


@pytest.mark.parametrize("time", [0.0, 4.47e17, 1.0e18])
def test_fixed_with_decay_matches_reference(time):
    fixed = Radiogenics.FixedRadiogenics(1.0e-11, 4.47e17, 0.0)
    expected = _reference_fixed_heating(time, _MASS, 1.0e-11, 4.47e17)
    assert fixed.calc_heating(time, _MASS) == pytest.approx(expected, rel=1e-12)


def test_fixed_half_life_halves_heating():
    fixed = Radiogenics.FixedRadiogenics(1.0e-11, 4.47e17, 0.0)
    assert fixed.calc_heating(4.47e17, _MASS) == pytest.approx(0.5 * fixed.calc_heating(0.0, _MASS), rel=1e-12)


@pytest.mark.parametrize("model_kind", ["fixed", "isotope"])
def test_overflow_before_reference_time_is_nan(model_kind):
    """Far before the reference time the decay exponential overflows; the result is NaN, not inf."""
    half_life = 1.0          # [s] tiny, so a modest time offset overflows exp()
    ref_time = 1.0e6         # [s]
    if model_kind == "fixed":
        model = Radiogenics.FixedRadiogenics(1.0e-11, half_life, ref_time)
        one_shot = Radiogenics.fixed(
            0.0,
            _MASS,
            fixed_heat_production=1.0e-11,
            average_half_life=half_life,
            ref_time=ref_time)
    else:
        isotopes = ([1.0e-5], [half_life], [1.0], [1.0e-6])
        model = Radiogenics.IsotopeRadiogenics(*isotopes, ref_time=ref_time)
        one_shot = Radiogenics.isotope(0.0, _MASS, *isotopes, ref_time=ref_time)
    assert math.isnan(model.calc_heating(0.0, _MASS))
    assert math.isnan(one_shot)
    assert np.all(np.isnan(model.calc_heating_vectorize_time(np.array([0.0, 1.0]), _MASS)))


# =====================================================================================================================
# Parameter getters
# =====================================================================================================================
def test_isotope_parameters():
    isotope = _isotope_model(ref_time=1.0e17)
    assert isotope.num_isotopes == 2
    assert isotope.ref_time == pytest.approx(1.0e17)
    assert list(isotope.heat_production) == pytest.approx(_HPR)
    assert list(isotope.half_lives) == pytest.approx(_HALF)
    assert list(isotope.mass_fracs) == pytest.approx(_FRAC)
    assert list(isotope.concentrations) == pytest.approx(_CONC)
    # Names are generated when not supplied.
    assert list(isotope.isotope_names) == ["isotope_0", "isotope_1"]


def test_isotope_explicit_names():
    assert list(_isotope_model(names=["U238", "Th232"]).isotope_names) == ["U238", "Th232"]


def test_fixed_parameters():
    fixed = Radiogenics.FixedRadiogenics(2.5e-11, 1.0e18, 3.0e17)
    assert fixed.fixed_heat_production == pytest.approx(2.5e-11)
    assert fixed.average_half_life == pytest.approx(1.0e18)
    assert fixed.ref_time == pytest.approx(3.0e17)


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
    assert Radiogenics.make_radiogenics(alias).model_name == canonical


def test_make_radiogenics_fixed_config():
    fixed = Radiogenics.make_radiogenics("fixed", {
        "fixed_heat_production_w_kg": 3.0e-11,
        "average_half_life_s": 1.0e18,
        "ref_time_s": 2.0e17,
    })
    assert fixed.fixed_heat_production == pytest.approx(3.0e-11)
    assert fixed.average_half_life == pytest.approx(1.0e18)
    assert fixed.ref_time == pytest.approx(2.0e17)


def test_make_radiogenics_isotope_explicit():
    isotope = Radiogenics.make_radiogenics("isotope", {
        "heat_production_w_kg": _HPR,
        "half_lives_s": _HALF,
        "mass_fracs": _FRAC,
        "concentrations": _CONC,
        "ref_time_s": 0.0,
    })
    assert isotope.num_isotopes == 2
    expected = _reference_isotope_heating(0.0, _MASS, _ISOTOPES, 0.0)
    assert isotope.calc_heating(0.0, _MASS) == pytest.approx(expected, rel=1e-12)


def test_make_radiogenics_isotope_named_dataset():
    """A named dataset is looked up and its Myr values converted to seconds."""
    isotope = Radiogenics.make_radiogenics("isotope", {"isotopes": "modern_day_chondritic"})
    assert isotope.num_isotopes == 4
    assert isotope.ref_time > 1.0e17
    assert all(half_life > 1.0e16 for half_life in isotope.half_lives)


def test_make_radiogenics_user_dataset_from_config():
    """A dataset added to the configuration's known_isotope_data is found and converted to seconds."""
    known = TidalPy.config["radiogenics"]["known_isotope_data"]
    known["user_test_dataset"] = {
        "ref_time": 100.0,
        "U238": {"hpr": 9.48e-5, "half_life": 4470.0, "iso_mass_fraction": 0.9928, "element_concentration": 0.012e-6},
    }
    try:
        isotope = Radiogenics.make_radiogenics("isotope", {"isotopes": "user_test_dataset"})
        assert isotope.num_isotopes == 1
        assert isotope.ref_time == pytest.approx(100.0 * seconds_per_myr)
        assert isotope.half_lives[0] == pytest.approx(4470.0 * seconds_per_myr)
        with pytest.raises(ValueError, match="unknown isotope dataset"):
            Radiogenics.make_radiogenics("isotope", {"isotopes": "not_a_dataset"})
    finally:
        del known["user_test_dataset"]


def test_make_radiogenics_inline_dataset():
    """An inline dataset dict in Myr is converted to seconds."""
    inline = {
        "ref_time": 4600.0,
        "U238": {"hpr": 9.48e-5, "half_life": 4470.0,
                 "iso_mass_fraction": 0.9928, "element_concentration": 0.012e-6},
    }
    isotope = Radiogenics.make_radiogenics("isotope", {"isotopes": inline})
    assert isotope.num_isotopes == 1
    assert isotope.ref_time == pytest.approx(4600.0 * _MYR_S)
    assert isotope.half_lives[0] == pytest.approx(4470.0 * _MYR_S)


def test_explicit_arrays_win_over_a_dataset_name():
    explicit = {"heat_production_w_kg": [1.0e-5], "half_lives_s": [1.0e30], "mass_fracs": [1.0],
                "concentrations": [1.0]}
    alone = Radiogenics.make_radiogenics("isotope", dict(explicit))
    with_name = Radiogenics.make_radiogenics("isotope", {**explicit, "isotopes": "modern_day_chondritic"})
    assert with_name.calc_heating(0.0, 1.0) == pytest.approx(alone.calc_heating(0.0, 1.0), rel=1e-12)
    chondritic = Radiogenics.make_radiogenics("isotope", {"isotopes": "modern_day_chondritic"})
    assert with_name.calc_heating(0.0, 1.0) != pytest.approx(chondritic.calc_heating(0.0, 1.0), rel=1e-3)


def test_reference_time_applies_to_a_dataset():
    """The dataset shifted 100 Myr later heats at 100 Myr as the unshifted one does at its own epoch."""
    default = Radiogenics.make_radiogenics("isotope", {"isotopes": "llri_and_slri"})
    shifted = Radiogenics.make_radiogenics("isotope", {"isotopes": "llri_and_slri", "ref_time_s": 100.0 * _MYR_S})
    assert shifted.calc_heating(100.0 * _MYR_S, 1.0) == pytest.approx(default.calc_heating(0.0, 1.0), rel=1e-9)


def test_make_radiogenics_adopted_object_is_usable(tmp_path):
    """An object built by the C++ enum factory is usable and round trips."""
    fixed = Radiogenics.make_radiogenics("fixed", {
        "fixed_heat_production_w_kg": 1.5e-11, "average_half_life_s": 5.0e17})
    assert fixed.fixed_heat_production == pytest.approx(1.5e-11)
    expected = _reference_fixed_heating(1.0e17, _MASS, 1.5e-11, 5.0e17)
    assert fixed.calc_heating(1.0e17, _MASS) == pytest.approx(expected, rel=1e-12)
    restored = _binary_copy(fixed, Radiogenics.FixedRadiogenics(), tmp_path)
    assert restored.get_config_dict() == fixed.get_config_dict()


# =====================================================================================================================
# Built-in literature isotope datasets
# =====================================================================================================================
def test_available_isotope_datasets():
    assert set(Radiogenics.available_isotope_datasets()) == {
        "modern_day_chondritic", "llri", "slri", "llri_and_slri", "bulk_silicate_earth"}


@pytest.mark.parametrize("name,n_isotopes,ref_time", [
    ("modern_day_chondritic", 4, 4600.0 * _MYR_S),
    ("llri", 4, 0.0),
    ("slri", 3, 0.0),
    ("llri_and_slri", 7, 0.0),
    ("bulk_silicate_earth", 4, 4600.0 * _MYR_S),
])
def test_isotope_dataset_contents(name, n_isotopes, ref_time):
    """Each dataset is an MKS dict with the expected isotope count and reference epoch."""
    # Present-epoch datasets quote abundances 4600 Myr after formation; the Castillo-Rogez sets quote formation.
    dataset = Radiogenics.isotope_dataset(name)
    assert set(dataset) == {
        "heat_production_w_kg", "half_lives_s", "mass_fracs",
        "concentrations", "isotope_names", "ref_time_s"}
    assert len(dataset["isotope_names"]) == n_isotopes
    assert dataset["ref_time_s"] == pytest.approx(ref_time, abs=1.0)
    assert all(half_life > 0.0 for half_life in dataset["half_lives_s"])


# Castillo-Rogez et al. (2007) Table 3: each isotope's concentration in ordinary chondrites at CAI formation [kg/kg],
# with 60Fe at the 60Fe/56Fe = 1e-6 end of its 22.5 to 225 ppb range.
_CASTILLO_ROGEZ_TABLE_3 = {
    "Al26": 600.0e-9, "Fe60": 225.0e-9, "Mn53": 25.7e-9,
    "K40": 1104.0e-9, "Th232": 53.8e-9, "U235": 8.2e-9, "U238": 26.2e-9}
_LONG_LIVED = {"U238", "U235", "Th232", "K40"}
_SHORT_LIVED = {"Al26", "Fe60", "Mn53"}


@pytest.mark.parametrize("name,expected_isotopes", [
    ("llri", _LONG_LIVED),
    ("slri", _SHORT_LIVED),
    ("llri_and_slri", _LONG_LIVED | _SHORT_LIVED),
])
def test_castillo_rogez_isotope_concentrations_reproduce_table_3(name, expected_isotopes):
    """Each isotope's mass fraction times concentration is its Table 3 concentration at formation."""
    # Table 3 already folds in the isotopic abundance, so applying the Table 4 or 5 abundance again would double it.
    dataset = Radiogenics.isotope_dataset(name)
    assert set(dataset["isotope_names"]) == expected_isotopes
    for isotope, mass_frac, concentration in zip(
            dataset["isotope_names"], dataset["mass_fracs"], dataset["concentrations"]):
        assert mass_frac * concentration == pytest.approx(_CASTILLO_ROGEZ_TABLE_3[isotope], rel=1e-12)
    # Long-lived isotopes carry their own concentration; short-lived ones their element's times the Table 5 ratio.
    fractions = dict(zip(dataset["isotope_names"], dataset["mass_fracs"]))
    for isotope in expected_isotopes & _LONG_LIVED:
        assert fractions[isotope] == 1.0
    for isotope, ratio in (("Al26", 5.0e-5), ("Fe60", 1.0e-6), ("Mn53", 1.0e-5)):
        if isotope in expected_isotopes:
            assert fractions[isotope] == ratio


@pytest.mark.parametrize("time_myr", [0.0, 0.5, 3.0, 10.0, 100.0, 4568.0])
def test_llri_and_slri_is_the_sum_of_llri_and_slri(time_myr):
    combined = Radiogenics.IsotopeRadiogenics.from_dataset("llri_and_slri")
    long_lived = Radiogenics.IsotopeRadiogenics.from_dataset("llri")
    short_lived = Radiogenics.IsotopeRadiogenics.from_dataset("slri")
    time = time_myr * _MYR_S
    assert combined.calc_heating(time, _MASS) == pytest.approx(
        long_lived.calc_heating(time, _MASS) + short_lived.calc_heating(time, _MASS), rel=1e-12)


def test_llri_heating_finite_at_formation():
    """Formation-epoch heating is the Table 3 sum (about 2.21e-7 W/kg, mostly Al26) and decays away within 100 Myr."""
    model = Radiogenics.IsotopeRadiogenics.from_dataset("llri_and_slri")
    dataset = Radiogenics.isotope_dataset("llri_and_slri")
    heating_formation = model.calc_heating(0.0, _MASS)
    heating_10myr = model.calc_heating(10.0 * _MYR_S, _MASS)
    heating_100myr = model.calc_heating(100.0 * _MYR_S, _MASS)
    assert math.isfinite(heating_formation)
    expected = sum(
        heat_production * _CASTILLO_ROGEZ_TABLE_3[isotope]
        for isotope, heat_production in zip(dataset["isotope_names"], dataset["heat_production_w_kg"]))
    assert heating_formation == pytest.approx(expected * _MASS, rel=1e-12)
    assert heating_formation / _MASS == pytest.approx(2.2128e-7, rel=1e-3)
    assert heating_10myr < heating_formation
    assert heating_100myr < heating_10myr
    assert heating_100myr < 1.0e-3 * heating_formation


def test_llri_heating_today_is_chondritic():
    """Decayed to the present, the long-lived set heats as ordinary chondrites do today, about 5e-12 W/kg."""
    model = Radiogenics.IsotopeRadiogenics.from_dataset("llri")
    heating_today = model.calc_heating(4568.0 * _MYR_S, _MASS) / _MASS
    assert 4.0e-12 < heating_today < 6.0e-12


@pytest.mark.parametrize("time_myr", [0.0, 3.0, 4568.0])
def test_legacy_config_llri_and_slri_matches(time_myr):
    """The legacy LLRI_and_SLRI dataset through the legacy formula heats as the built-in set does."""
    # The legacy dataset is in Myr, so it was evaluated at the time in Myr.
    legacy_heating = float(_LEGACY_REFERENCE[f"legacy_config_llri_and_slri__time_myr_{time_myr:g}"])
    heating = Radiogenics.IsotopeRadiogenics.from_dataset("llri_and_slri").calc_heating(time_myr * _MYR_S, _MASS)
    assert legacy_heating == pytest.approx(heating, rel=1e-12)


def test_from_dataset_matches_factory():
    from_dataset = Radiogenics.IsotopeRadiogenics.from_dataset("modern_day_chondritic")
    from_factory = Radiogenics.make_radiogenics("isotope", {"isotopes": "modern_day_chondritic"})
    assert isinstance(from_dataset, Radiogenics.IsotopeRadiogenics)
    assert list(from_dataset.isotope_names) == list(from_factory.isotope_names)
    assert from_dataset.calc_heating(from_dataset.ref_time, _MASS) == pytest.approx(
        from_factory.calc_heating(from_factory.ref_time, _MASS), rel=1e-12)


def test_dataset_heating_cross_check():
    dataset = Radiogenics.isotope_dataset("modern_day_chondritic")
    model = Radiogenics.IsotopeRadiogenics.from_dataset("modern_day_chondritic")
    time = dataset["ref_time_s"]
    isotopes = (dataset["heat_production_w_kg"], dataset["half_lives_s"], dataset["mass_fracs"],
                dataset["concentrations"])
    expected = _reference_isotope_heating(time, _MASS, isotopes, dataset["ref_time_s"])
    assert model.calc_heating(time, _MASS) == pytest.approx(expected, rel=1e-12)


# =====================================================================================================================
# Config dict and binary round trip
# =====================================================================================================================
@pytest.mark.parametrize("factory,keys", [
    (Radiogenics.OffRadiogenics, {"model"}),
    (Radiogenics.FixedRadiogenics, {"model", "fixed_heat_production_w_kg", "average_half_life_s", "ref_time_s"}),
    (_isotope_model, {"model", "heat_production_w_kg", "half_lives_s", "mass_fracs", "concentrations",
                      "isotope_names", "ref_time_s"}),
], ids=["off", "fixed", "isotope"])
def test_config_dict_keys(factory, keys):
    assert set(factory().get_config_dict()) == keys


def test_save_config_writes_toml(tmp_path):
    import toml
    path = str(tmp_path / "radio.toml")
    Radiogenics.FixedRadiogenics(2.5e-11, 1.0e18, 3.0e17).save_config(path)
    loaded = toml.load(path)
    assert loaded["model"] == "fixed"
    assert loaded["fixed_heat_production_w_kg"] == pytest.approx(2.5e-11)
    assert loaded["average_half_life_s"] == pytest.approx(1.0e18)


@pytest.mark.parametrize("factory,model_name,num_isotopes", [
    (Radiogenics.OffRadiogenics, "off", None),
    (lambda: Radiogenics.FixedRadiogenics(2.5e-11, 1.0e18, 3.0e17), "fixed", None),
    (lambda: _isotope_model(ref_time=1.0e17, names=["U238", "Th232"]), "isotope", 2),
    (lambda: Radiogenics.IsotopeRadiogenics.from_dataset("llri_and_slri"), "isotope", 7),
], ids=["off", "fixed", "isotope", "dataset"])
def test_binary_round_trip(factory, model_name, num_isotopes, tmp_path):
    """A binary round trip keeps the model, its parameters, and its heating."""
    original = factory()
    restored = _binary_copy(original, type(original)(), tmp_path)
    assert restored.model_name == model_name
    assert restored.get_config_dict() == original.get_config_dict()
    assert restored.calc_heating(0.0, _MASS) == pytest.approx(original.calc_heating(0.0, _MASS), rel=1e-12)
    if num_isotopes is not None:
        assert restored.num_isotopes == num_isotopes
        assert list(restored.isotope_names) == list(original.isotope_names)
        assert list(restored.half_lives) == pytest.approx(list(original.half_lives))


def test_isotope_binary_payload_size_mismatch_raises(tmp_path):
    """An isotope record whose header claims more payload than its isotope list is refused, unread."""
    path = str(tmp_path / "isotope.tpyb")
    _isotope_model(ref_time=1.0e17, names=["U238", "Th232"]).save_binary(path)
    with open(path, "rb") as binary_file:
        record = bytearray(binary_file.read())
    # The header's payload size is the little-endian uint64 at byte 12.
    payload_size = int.from_bytes(record[12:20], "little")
    record[12:20] = (payload_size + 8).to_bytes(8, "little")
    record += bytes(8)
    with open(path, "wb") as binary_file:
        binary_file.write(record)
    restored = Radiogenics.IsotopeRadiogenics()
    with pytest.raises(IOError, match="payload bytes"):
        restored.load_binary(path)
    assert restored.num_isotopes == 0


# =====================================================================================================================
# Vectorized methods and direct convenience functions
# =====================================================================================================================
def test_convenience_scalar_matches_class():
    got = Radiogenics.fixed(1.0e17, _MASS, *_DECAYING_FIXED)
    assert isinstance(got, float)
    assert got == pytest.approx(Radiogenics.FixedRadiogenics(*_DECAYING_FIXED).calc_heating(1.0e17, _MASS))


def test_convenience_vectorize_time():
    times = np.array([0.0, 1.0e17, 5.0e17, 1.0e18])
    got = Radiogenics.isotope(times, _MASS, *_ISOTOPES, ref_time=0.0)
    assert isinstance(got, np.ndarray)
    assert got.shape == (4,)
    assert got.dtype == np.float64
    model = _isotope_model()
    assert got == pytest.approx(np.array([model.calc_heating(time, _MASS) for time in times]))


def test_convenience_vectorize_mass():
    masses = np.array([1.0e21, 5.0e21, 1.0e22])
    got = Radiogenics.fixed(0.0, masses, *_STEADY_FIXED)
    assert got.shape == (3,)
    model = Radiogenics.FixedRadiogenics(*_STEADY_FIXED)
    assert got == pytest.approx(np.array([model.calc_heating(0.0, mass) for mass in masses]))


def test_convenience_vectorize_all_and_broadcast():
    """All-array, mixed, and 2-D broadcast inputs match element-wise scalar calls."""
    times = np.array([0.0, 1.0e17, 5.0e17])
    masses = np.array([1.0e21, 5.0e21, 1.0e22])
    fixed = Radiogenics.FixedRadiogenics(*_DECAYING_FIXED)
    got_all = Radiogenics.fixed(times, masses, *_DECAYING_FIXED)
    got_mixed = Radiogenics.fixed(times, _MASS, *_DECAYING_FIXED)
    assert got_all.shape == (3,)
    assert got_mixed.shape == (3,)
    assert got_all == pytest.approx(np.array([fixed.calc_heating(time, mass) for time, mass in zip(times, masses)]))
    assert got_mixed == pytest.approx(np.array([fixed.calc_heating(time, _MASS) for time in times]))

    isotope = _isotope_model()
    got_2d = Radiogenics.isotope(times[:, None], masses[None, :], *_ISOTOPES, ref_time=0.0)
    assert got_2d.shape == (3, 3)
    assert got_2d == pytest.approx(np.array([[isotope.calc_heating(time, mass) for mass in masses] for time in times]))


def test_convenience_preserves_2d_shape():
    times = np.array([[0.0, 1.0e17], [2.0e17, 3.0e17]])
    assert Radiogenics.fixed(times, _MASS, *_DECAYING_FIXED).shape == (2, 2)


def test_class_vectorize_methods_match_scalar():
    isotope = _isotope_model()
    times = np.array([0.0, 1.0e17, 5.0e17])
    masses = np.array([1.0e21, 5.0e21, 1.0e22])
    assert isotope.calc_heating_vectorize_time(times, _MASS) == pytest.approx(
        np.array([isotope.calc_heating(time, _MASS) for time in times]))
    assert isotope.calc_heating_vectorize_mass(0.0, masses) == pytest.approx(
        np.array([isotope.calc_heating(0.0, mass) for mass in masses]))
    assert isotope.calc_heating_vectorize_all(times, masses) == pytest.approx(
        np.array([isotope.calc_heating(time, mass) for time, mass in zip(times, masses)]))
