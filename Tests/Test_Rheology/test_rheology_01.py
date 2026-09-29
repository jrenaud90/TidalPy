"""Rheology models: complex moduli against reference compliances, factory, vectorization, and I/O."""

import math

import numpy as np
import pytest

from TidalPy import Rheology
from TidalPy.Utilities.classes.classes import PhysicsBase, TidalPyBaseClass

_ETA = 1.0e20      # viscosity [Pa s]
_MU = 5.0e10       # shear modulus [Pa]
_OMEGA = 1.0e-5    # forcing frequency [rad/s]

_MODEL_CLASSES = [
    ("elastic",  "Elastic"),
    ("viscous",  "Viscous"),
    ("voigt",    "Voigt"),
    ("maxwell",  "Maxwell"),
    ("burgers",  "Burgers"),
    ("andrade",  "Andrade"),
    ("sundberg", "Sundberg"),
    ("zener",    "Zener"),
]
_MODEL_NAMES = [name for name, _ in _MODEL_CLASSES]


# =====================================================================================================================
# Reference compliances at (_ETA, _MU, _OMEGA), re-implementing the TidalPy 0.7 compliance_models math.
# =====================================================================================================================
def _maxwell_compliance():
    return complex(1.0 / _MU, -1.0 / (_ETA * _OMEGA))


def _voigt_compliance(modulus_frac=5.0, viscosity_frac=0.02):
    voigt_compliance = (1.0 / _MU) / modulus_frac
    voigt_viscosity = viscosity_frac * _ETA
    denominator = (voigt_compliance * voigt_viscosity * _OMEGA) ** 2 + 1.0
    return complex(voigt_compliance / denominator, -(voigt_compliance ** 2) * voigt_viscosity * _OMEGA / denominator)


def _andrade_compliance(alpha=0.3, zeta=1.0):
    compliance = 1.0 / _MU
    andrade_term = compliance * ((compliance * _ETA * _OMEGA * zeta) ** (-alpha)) * math.gamma(1.0 + alpha)
    andrade = complex(math.cos(alpha * math.pi / 2.0) * andrade_term, -math.sin(alpha * math.pi / 2.0) * andrade_term)
    return _maxwell_compliance() + andrade


def _zener_modulus(relaxed_frac=0.5, frequency=_OMEGA):
    """Standard linear solid: relaxed spring in parallel with a Maxwell arm (Nowick and Berry 1972)."""
    arm = (1.0 - relaxed_frac) * _MU
    iwt = 1j * frequency * _ETA / arm
    return relaxed_frac * _MU + arm * iwt / (1.0 + iwt)


_REF_COMPLIANCE = {
    "elastic":  lambda: complex(1.0 / _MU, 0.0),
    "viscous":  lambda: complex(0.0, -1.0 / (_ETA * _OMEGA)),
    "voigt":    _voigt_compliance,
    "maxwell":  _maxwell_compliance,
    "burgers":  lambda: _maxwell_compliance() + _voigt_compliance(),
    "andrade":  _andrade_compliance,
    "sundberg": lambda: _andrade_compliance() + _voigt_compliance(),
    "zener":    lambda: 1.0 / _zener_modulus(),
}


def _make(name):
    """Construct a model by canonical name with default parameters."""
    return getattr(Rheology, dict(_MODEL_CLASSES)[name])()


def _assert_complex_close(got, expected):
    assert got.real == pytest.approx(expected.real, rel=1e-12)
    assert got.imag == pytest.approx(expected.imag, rel=1e-12)


def _binary_copy(model, blank, tmp_path):
    """Save a model to a binary file and load it into a blank model."""
    path = str(tmp_path / "rheology.tpyb")
    model.save_binary(path)
    blank.load_binary(path)
    return blank


# =====================================================================================================================
# Model identity
# =====================================================================================================================
@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_model_name(name, cls):
    assert _make(name).model_name == name


@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_make_rheology_returns_rich_subclass(name, cls):
    model = Rheology.make_rheology(name)
    assert type(model).__name__ == cls
    assert isinstance(model, getattr(Rheology, cls))
    assert isinstance(model, Rheology.RheologyBase)


@pytest.mark.parametrize("name", _MODEL_NAMES)
def test_isinstance_chain(name):
    model = _make(name)
    assert isinstance(model, Rheology.RheologyBase)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)


def test_model_name_is_read_only():
    model = Rheology.Maxwell()
    with pytest.raises(AttributeError):
        model.model_name = "andrade"
    assert model.model_name == "maxwell"


@pytest.mark.parametrize("call,error", [
    (lambda: Rheology.RheologyBase(), TypeError),
    (lambda: Rheology.make_rheology("not_a_real_model"), ValueError),
    (lambda: Rheology.Maxwell().calc_complex_modulus_vectorize_all(np.ones(3), np.ones(2), np.ones(3)), ValueError),
    (lambda: Rheology.Maxwell().calc_complex_modulus_vectorize_modulus(np.ones(3), np.ones(2), _OMEGA), ValueError),
    (lambda: Rheology.Maxwell().load_binary("does_not_exist_12345.tpyb"), FileNotFoundError),
], ids=["abstract_base", "unknown_model", "vectorize_all_size_mismatch", "vectorize_modulus_size_mismatch",
        "missing_binary"])
def test_invalid_use_raises(call, error):
    with pytest.raises(error):
        call()


# =====================================================================================================================
# Complex modulus physics
# =====================================================================================================================
@pytest.mark.parametrize("name", _MODEL_NAMES)
def test_complex_modulus_matches_reference(name):
    """The complex modulus is 1 / reference compliance."""
    _assert_complex_close(_make(name).calc_complex_modulus(_MU, _ETA, _OMEGA), 1.0 / _REF_COMPLIANCE[name]())


@pytest.mark.parametrize("name", _MODEL_NAMES)
def test_loss_modulus_nonnegative(name):
    # Elastic is about zero; every other model is positive.
    assert _make(name).calc_complex_modulus(_MU, _ETA, _OMEGA).imag >= -1.0e-6


@pytest.mark.parametrize("name", _MODEL_NAMES)
def test_zero_frequency_finite(name):
    modulus = _make(name).calc_complex_modulus(_MU, _ETA, 0.0)
    assert math.isfinite(modulus.real)
    assert math.isfinite(modulus.imag)


def test_elastic_is_real_lossless_and_frequency_independent():
    elastic = _make("elastic")
    modulus = elastic.calc_complex_modulus(_MU, _ETA, _OMEGA)
    assert modulus.real == pytest.approx(_MU)
    assert modulus.imag == pytest.approx(0.0)
    assert elastic.calc_complex_modulus(_MU, _ETA, 1.0e-7) == elastic.calc_complex_modulus(_MU, _ETA, 1.0e-3)


def test_viscous_modulus_purely_dissipative():
    """The viscous modulus is i * viscosity * frequency."""
    modulus = _make("viscous").calc_complex_modulus(_MU, _ETA, _OMEGA)
    assert modulus.real == pytest.approx(0.0)
    assert modulus.imag == pytest.approx(_ETA * _OMEGA)


@pytest.mark.parametrize("frequency", [_OMEGA, 0.0])
@pytest.mark.parametrize("name", ["maxwell", "burgers", "andrade", "sundberg", "zener"])
def test_infinite_viscosity_is_elastic(name, frequency):
    """An infinite viscosity (the cold limit) locks every dashpot, so the series models are elastic."""
    # Also at zero frequency, where the viscous term would otherwise be inf * 0.
    modulus = _make(name).calc_complex_modulus(_MU, math.inf, frequency)
    assert modulus.real == pytest.approx(_MU, rel=1.0e-12)
    assert modulus.imag == 0.0


def test_infinite_viscosity_voigt_is_rigid():
    """A Voigt solid with infinite viscosity is infinitely stiff under forcing and a spring at rest."""
    voigt = _make("voigt")
    forced = voigt.calc_complex_modulus(_MU, math.inf, _OMEGA)
    assert forced.real == pytest.approx(5.0 * _MU, rel=1.0e-12)
    assert math.isinf(forced.imag) and forced.imag > 0.0
    static = voigt.calc_complex_modulus(_MU, math.inf, 0.0)
    assert static == pytest.approx(complex(5.0 * _MU, 0.0), rel=1.0e-12)


# =====================================================================================================================
# Parameters and factory
# =====================================================================================================================
@pytest.mark.parametrize("cls_name,params", [
    ("Voigt", {"voigt_modulus_frac": 0.3, "voigt_viscosity_frac": 0.05}),
    ("Andrade", {"alpha": 0.25, "zeta": 2.0}),
    ("Sundberg", {"alpha": 0.4, "zeta": 3.0, "voigt_modulus_frac": 0.15, "voigt_viscosity_frac": 0.03}),
    ("Zener", {"relaxed_modulus_frac": 0.8}),
])
def test_positional_parameters(cls_name, params):
    model = getattr(Rheology, cls_name)(*params.values())
    for name, value in params.items():
        assert getattr(model, name) == pytest.approx(value), name


def test_andrade_parameters_affect_modulus():
    modulus_default = Rheology.Andrade(0.3, 1.0).calc_complex_modulus(_MU, _ETA, _OMEGA)
    modulus_changed = Rheology.Andrade(0.5, 1.0).calc_complex_modulus(_MU, _ETA, _OMEGA)
    assert modulus_default != modulus_changed
    _assert_complex_close(modulus_changed, 1.0 / _andrade_compliance(alpha=0.5, zeta=1.0))


@pytest.mark.parametrize("alias,canonical", [
    ("off",             "elastic"),
    ("Elastic",         "elastic"),
    ("newton",          "viscous"),
    ("VISCOUS",         "viscous"),
    ("voigt-kelvin",    "voigt"),
    ("voigt_kelvin",    "voigt"),
    ("Maxwell",         "maxwell"),
    ("burgers",         "burgers"),
    ("andrade",         "andrade"),
    ("Sundberg-Cooper", "sundberg"),
    ("sundberg_cooper", "sundberg"),
    ("Zener",           "zener"),
    ("sls",             "zener"),
    ("Standard_Linear_Solid", "zener"),
])
def test_make_rheology_aliases(alias, canonical):
    assert Rheology.make_rheology(alias).model_name == canonical


def test_make_rheology_with_config():
    sundberg = Rheology.make_rheology("sundberg", {
        "alpha": 0.45, "zeta": 2.5,
        "voigt_modulus_frac": 0.18, "voigt_viscosity_frac": 0.04,
    })
    assert sundberg.alpha == pytest.approx(0.45)
    assert sundberg.zeta == pytest.approx(2.5)
    assert sundberg.voigt_modulus_frac == pytest.approx(0.18)
    assert sundberg.voigt_viscosity_frac == pytest.approx(0.04)


def test_make_rheology_adopted_object_is_usable(tmp_path):
    """An object built by the C++ enum factory keeps its parameters, computes, and round trips."""
    andrade = Rheology.make_rheology("andrade", {"alpha": 0.42, "zeta": 1.7})
    assert andrade.alpha == pytest.approx(0.42)
    assert andrade.zeta == pytest.approx(1.7)
    _assert_complex_close(andrade.calc_complex_modulus(_MU, _ETA, _OMEGA),
                          1.0 / _andrade_compliance(alpha=0.42, zeta=1.7))
    restored = _binary_copy(andrade, Rheology.Andrade(), tmp_path)
    assert restored.get_config_dict() == andrade.get_config_dict()


# =====================================================================================================================
# Config dict and binary round trip
# =====================================================================================================================
@pytest.mark.parametrize("cls_name,keys", [
    ("Elastic", {"model"}),
    ("Maxwell", {"model"}),
    ("Voigt", {"model", "voigt_modulus_frac", "voigt_viscosity_frac"}),
    ("Andrade", {"model", "alpha", "zeta"}),
    ("Sundberg", {"model", "alpha", "zeta", "voigt_modulus_frac", "voigt_viscosity_frac"}),
    ("Zener", {"model", "relaxed_modulus_frac"}),
])
def test_config_dict_keys(cls_name, keys):
    assert set(getattr(Rheology, cls_name)().get_config_dict()) == keys


def test_save_config_writes_toml(tmp_path):
    import toml
    path = str(tmp_path / "rheo.toml")
    Rheology.Sundberg(0.4, 2.0, 0.15, 0.03).save_config(path)
    loaded = toml.load(path)
    assert loaded["model"] == "sundberg"
    assert loaded["alpha"] == pytest.approx(0.4)
    assert loaded["voigt_modulus_frac"] == pytest.approx(0.15)


# Non-default constructor arguments for the models that take any.
_NON_DEFAULT_ARGS = {"voigt": (0.33, 0.07), "burgers": (0.33, 0.07), "andrade": (0.42, 1.7),
                     "sundberg": (0.42, 1.7, 0.33, 0.07), "zener": (0.23,)}


@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_binary_round_trip(name, cls, tmp_path):
    model_class = getattr(Rheology, cls)
    original = model_class(*_NON_DEFAULT_ARGS.get(name, ()))
    restored = _binary_copy(original, model_class(), tmp_path)
    assert restored.model_name == name
    assert restored.get_config_dict() == original.get_config_dict()


# =====================================================================================================================
# Vectorized methods and direct convenience functions
# =====================================================================================================================
@pytest.mark.parametrize("name", _MODEL_NAMES)
def test_convenience_scalar_matches_class(name):
    # The convenience function takes the class method's order: (modulus, viscosity, frequency).
    got = getattr(Rheology, name)(_MU, _ETA, _OMEGA)
    assert isinstance(got, complex)
    assert got == pytest.approx(_make(name).calc_complex_modulus(_MU, _ETA, _OMEGA))


@pytest.mark.parametrize("name", _MODEL_NAMES)
def test_convenience_vectorize_modulus(name):
    moduli = np.array([1.0e10, 5.0e10, 1.0e11])
    viscosities = np.array([1.0e19, 1.0e20, 1.0e21])
    got = getattr(Rheology, name)(moduli, viscosities, _OMEGA)
    assert isinstance(got, np.ndarray)
    assert got.shape == (3,)
    assert got.dtype == np.complex128
    model = _make(name)
    assert got == pytest.approx(np.array(
        [model.calc_complex_modulus(modulus, viscosity, _OMEGA) for modulus, viscosity in zip(moduli, viscosities)]))


def test_convenience_vectorize_frequency_with_params():
    """Array frequency works and model parameters reach the stack-allocated C++ model."""
    frequencies = np.array([1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4])
    got = Rheology.andrade(
        _MU,
        _ETA,
        frequencies,
        alpha=0.3,
        zeta=1.0)
    assert got.shape == (4,)
    model = Rheology.Andrade(0.3, 1.0)
    assert got == pytest.approx(np.array([model.calc_complex_modulus(_MU, _ETA, freq) for freq in frequencies]))
    tuned_got = Rheology.andrade(
        _MU,
        _ETA,
        _OMEGA,
        alpha=0.42,
        zeta=1.7)
    assert tuned_got == pytest.approx(Rheology.Andrade(0.42, 1.7).calc_complex_modulus(_MU, _ETA, _OMEGA))


def test_convenience_vectorize_all_and_broadcast():
    """All-array, mixed (broadcast), and 2-D inputs return the broadcast shape."""
    frequencies = np.array([1.0e-6, 1.0e-5, 1.0e-4])
    moduli = np.array([1.0e10, 5.0e10, 1.0e11])
    viscosities = np.array([1.0e19, 1.0e20, 1.0e21])
    assert Rheology.maxwell(moduli, viscosities, frequencies).shape == (3,)
    got_mixed = Rheology.maxwell(moduli, _ETA, frequencies)
    assert got_mixed.shape == (3,)
    maxwell = Rheology.Maxwell()
    assert got_mixed == pytest.approx(np.array(
        [maxwell.calc_complex_modulus(modulus, _ETA, freq) for modulus, freq in zip(moduli, frequencies)]))
    moduli_2d = np.array([[1.0e10, 2.0e10], [3.0e10, 4.0e10]])
    assert Rheology.maxwell(moduli_2d, np.full((2, 2), 1.0e20), _OMEGA).shape == (2, 2)


def test_class_vectorize_methods_match_scalar():
    sundberg = Rheology.Sundberg(0.3, 1.0, 0.2, 0.02)
    moduli = np.array([1.0e10, 5.0e10, 1.0e11])
    viscosities = np.array([1.0e19, 1.0e20, 1.0e21])
    frequencies = np.array([1.0e-6, 1.0e-5, 1.0e-4])
    assert sundberg.calc_complex_modulus_vectorize_modulus(moduli, viscosities, _OMEGA) == pytest.approx(np.array(
        [sundberg.calc_complex_modulus(mod, visc, _OMEGA) for mod, visc in zip(moduli, viscosities)]))
    assert sundberg.calc_complex_modulus_vectorize_frequency(_MU, _ETA, frequencies) == pytest.approx(np.array(
        [sundberg.calc_complex_modulus(_MU, _ETA, freq) for freq in frequencies]))
    assert sundberg.calc_complex_modulus_vectorize_all(moduli, viscosities, frequencies) == pytest.approx(np.array(
        [sundberg.calc_complex_modulus(*args) for args in zip(moduli, viscosities, frequencies)]))


def test_vectorized_keeps_each_complex_value():
    """Vectorized results equal the scalar ones exactly, including a locked Voigt dashpot's finite storage modulus."""
    viscosities = np.array([1.0e20, np.inf])
    got = Rheology.voigt(np.full(2, _MU), viscosities, _OMEGA)
    for index in range(2):
        expected = Rheology.Voigt().calc_complex_modulus(_MU, viscosities[index], _OMEGA)
        assert got[index].real == expected.real
        assert got[index].imag == expected.imag
    assert np.isfinite(got[1].real)


def test_read_only_inputs_are_accepted():
    def read_only(values):
        array = np.array(values, dtype=np.float64)
        array.setflags(write=False)
        return array

    modulus = read_only([5.0e10, 6.0e10])
    viscosity = read_only([1.0e19, 1.0e20])
    frequency = read_only([1.0e-5, 2.0e-5])
    direct = Rheology.maxwell(modulus, viscosity, frequency)
    method = Rheology.Maxwell().calc_complex_modulus_vectorize_all(modulus, viscosity, frequency)
    assert np.allclose(direct, method)


# =====================================================================================================================
# Zener (standard linear solid)
# =====================================================================================================================
@pytest.mark.parametrize("frequency", [1.0e-9, _OMEGA, 1.0e-2])
def test_zener_with_no_relaxed_spring_is_maxwell(frequency):
    _assert_complex_close(Rheology.Zener(0.0).calc_complex_modulus(_MU, _ETA, frequency),
                          Rheology.Maxwell().calc_complex_modulus(_MU, _ETA, frequency))


def test_zener_with_a_fully_relaxed_spring_is_elastic():
    assert Rheology.Zener(1.0).calc_complex_modulus(_MU, _ETA, _OMEGA) == complex(_MU, 0.0)


@pytest.mark.parametrize("relaxed_frac", [0.0, 0.3, 0.9])
def test_zener_limits_and_loss_peak(relaxed_frac):
    """It relaxes to r M at zero frequency, stays at M at high frequency, and peaks at omega tau = 1."""
    zener = Rheology.Zener(relaxed_frac)
    assert zener.calc_complex_modulus(_MU, _ETA, 0.0) == complex(relaxed_frac * _MU, 0.0)
    fast = zener.calc_complex_modulus(_MU, _ETA, 1.0e10)
    assert fast.real == pytest.approx(_MU, rel=1.0e-12)
    peak_frequency = (1.0 - relaxed_frac) * _MU / _ETA
    peak = zener.calc_complex_modulus(_MU, _ETA, peak_frequency)
    assert peak.real == pytest.approx(0.5 * (1.0 + relaxed_frac) * _MU, rel=1.0e-12)
    assert peak.imag == pytest.approx(0.5 * (1.0 - relaxed_frac) * _MU, rel=1.0e-12)
    for factor in (0.5, 2.0):
        assert zener.calc_complex_modulus(_MU, _ETA, factor * peak_frequency).imag < peak.imag


@pytest.mark.parametrize("frequency", [1.0e-300, 1.0e300])
def test_zener_extreme_frequencies_are_finite(frequency):
    modulus = Rheology.Zener(0.4).calc_complex_modulus(_MU, _ETA, frequency)
    assert math.isfinite(modulus.real) and math.isfinite(modulus.imag)
    assert modulus.imag >= 0.0


def test_zener_matches_reference_at_other_fractions():
    _assert_complex_close(Rheology.Zener(0.85).calc_complex_modulus(_MU, _ETA, _OMEGA), _zener_modulus(0.85))
    _assert_complex_close(Rheology.zener(_MU, _ETA, _OMEGA, relaxed_modulus_frac=0.85), _zener_modulus(0.85))


@pytest.mark.parametrize("relaxed_frac", [-0.1, 1.5, math.nan])
def test_zener_rejects_a_fraction_outside_zero_to_one(relaxed_frac):
    with pytest.raises(ValueError, match="relaxed_modulus_frac"):
        Rheology.Zener(relaxed_frac)
    with pytest.raises(ValueError, match="relaxed_modulus_frac"):
        Rheology.make_rheology("zener", {"relaxed_modulus_frac": relaxed_frac})
