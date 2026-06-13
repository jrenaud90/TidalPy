"""
Tests for TidalPy.rheology_x.rheology.

Covers construction, model names and aliases, complex-compliance physics
(cross-checked against reference formulas), complex-modulus inversion, the
``make_rheology`` factory, config dicts, binary round-trip, TOML config save,
and the isinstance chain.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import math
import os
import tempfile

import pytest


# =====================================================================================================================
# Imports
# =====================================================================================================================
def _import_rheology():
    try:
        from TidalPy import rheology_x as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.rheology_x not compiled — run uv pip install first."
        )


# Representative MKS inputs (silicate-mantle-ish, low forcing frequency)
_ETA   = 1.0e20    # viscosity   [Pa·s]
_MU    = 5.0e10    # shear modulus [Pa]
_OMEGA = 1.0e-5    # forcing frequency [rad/s]


# =====================================================================================================================
# Reference complex-modulus formulas.
#
# The reference element compliances independently re-implement the legacy
# TidalPy/rheology/complex_compliance/compliance_models.py math; the reference
# complex modulus is then mu* = 1 / J* (the quantity the public API returns).
# =====================================================================================================================
def _ref_compliance_elastic(eta, mu, omega):
    return complex(1.0 / mu, 0.0)


def _ref_compliance_viscous(eta, mu, omega):
    return complex(0.0, -1.0 / (eta * omega))


def _ref_compliance_maxwell(eta, mu, omega):
    return complex(1.0 / mu, -1.0 / (eta * omega))


def _ref_compliance_voigt(eta, mu, omega, mf=5.0, vf=0.02):
    vc = (1.0 / mu) / mf
    vv = vf * eta
    den = (vc * vv * omega) ** 2 + 1.0
    return complex(vc / den, -(vc ** 2) * vv * omega / den)


def _ref_compliance_burgers(eta, mu, omega, mf=5.0, vf=0.02):
    return _ref_compliance_maxwell(eta, mu, omega) + _ref_compliance_voigt(eta, mu, omega, mf, vf)


def _ref_compliance_andrade(eta, mu, omega, alpha=0.3, zeta=1.0):
    j = 1.0 / mu
    at = j * eta * omega * zeta
    ct = j * (at ** (-alpha)) * math.gamma(1.0 + alpha)
    andrade = complex(math.cos(alpha * math.pi / 2.0) * ct,
                      -math.sin(alpha * math.pi / 2.0) * ct)
    return _ref_compliance_maxwell(eta, mu, omega) + andrade


def _ref_compliance_sundberg(eta, mu, omega, alpha=0.3, zeta=1.0, mf=5.0, vf=0.02):
    return (_ref_compliance_andrade(eta, mu, omega, alpha, zeta)
            + _ref_compliance_voigt(eta, mu, omega, mf, vf))


_REF_COMPLIANCE = {
    "elastic":  _ref_compliance_elastic,
    "viscous":  _ref_compliance_viscous,
    "voigt":    _ref_compliance_voigt,
    "maxwell":  _ref_compliance_maxwell,
    "burgers":  _ref_compliance_burgers,
    "andrade":  _ref_compliance_andrade,
    "sundberg": _ref_compliance_sundberg,
}


def _ref_modulus(name, eta, mu, omega):
    """Reference complex modulus = 1 / reference compliance."""
    return 1.0 / _REF_COMPLIANCE[name](eta, mu, omega)

_MODEL_CLASSES = [
    ("elastic",  "Elastic"),
    ("viscous",  "Viscous"),
    ("voigt",    "Voigt"),
    ("maxwell",  "Maxwell"),
    ("burgers",  "Burgers"),
    ("andrade",  "Andrade"),
    ("sundberg", "Sundberg"),
]


def _make(model):
    """Construct a model by canonical name with default parameters."""
    mod = _import_rheology()
    return getattr(mod, dict(_MODEL_CLASSES)[model])()


# =====================================================================================================================
# Construction & names
# =====================================================================================================================
@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_model_name(name, cls):
    """Each model reports its canonical lower-case model name."""
    inst = _make(name)
    assert inst.model_name == name


def test_abstract_base_not_instantiable():
    """RheologyBase is abstract and cannot be constructed directly."""
    mod = _import_rheology()
    with pytest.raises(TypeError):
        mod.RheologyBase()


# =====================================================================================================================
# Complex modulus physics (cross-check vs reference formulas)
# =====================================================================================================================
@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_complex_modulus_matches_reference(name, cls):
    """Complex modulus matches the independent reference (1 / reference compliance)."""
    inst = _make(name)
    got = inst.calc_complex_modulus(_MU, _ETA, _OMEGA)
    expected = _ref_modulus(name, _ETA, _MU, _OMEGA)
    assert got.real == pytest.approx(expected.real, rel=1e-12)
    assert got.imag == pytest.approx(expected.imag, rel=1e-12)


def test_elastic_modulus_real_and_lossless():
    """Elastic modulus equals the static modulus with no loss component."""
    e = _make("elastic")
    mu = e.calc_complex_modulus(_MU, _ETA, _OMEGA)
    assert mu.real == pytest.approx(_MU)
    assert mu.imag == pytest.approx(0.0)


def test_elastic_frequency_independent():
    """Elastic response does not depend on forcing frequency."""
    e = _make("elastic")
    g1 = e.calc_complex_modulus(_MU, _ETA, 1.0e-7)
    g2 = e.calc_complex_modulus(_MU, _ETA, 1.0e-3)
    assert g1 == g2


def test_viscous_modulus_purely_dissipative():
    """Viscous modulus is purely imaginary: mu* = i * viscosity * frequency."""
    v = _make("viscous")
    mu = v.calc_complex_modulus(_MU, _ETA, _OMEGA)
    assert mu.real == pytest.approx(0.0)
    assert mu.imag == pytest.approx(_ETA * _OMEGA)


def test_loss_modulus_nonnegative():
    """Every dissipative model has a non-negative loss (imaginary) modulus."""
    for name, _cls in _MODEL_CLASSES:
        inst = _make(name)
        mu = inst.calc_complex_modulus(_MU, _ETA, _OMEGA)
        assert mu.imag >= -1.0e-6  # elastic is ~0; all others positive


def test_zero_frequency_finite():
    """At zero forcing frequency every model returns a finite complex modulus."""
    for name, _cls in _MODEL_CLASSES:
        inst = _make(name)
        mu = inst.calc_complex_modulus(_MU, _ETA, 0.0)
        assert math.isfinite(mu.real)
        assert math.isfinite(mu.imag)


# =====================================================================================================================
# Parameter getters
# =====================================================================================================================
def test_voigt_parameters():
    mod = _import_rheology()
    v = mod.Voigt(0.3, 0.05)
    assert v.voigt_modulus_frac == pytest.approx(0.3)
    assert v.voigt_viscosity_frac == pytest.approx(0.05)


def test_andrade_parameters():
    mod = _import_rheology()
    a = mod.Andrade(0.25, 2.0)
    assert a.alpha == pytest.approx(0.25)
    assert a.zeta == pytest.approx(2.0)


def test_sundberg_parameters():
    mod = _import_rheology()
    s = mod.Sundberg(0.4, 3.0, 0.15, 0.03)
    assert s.alpha == pytest.approx(0.4)
    assert s.zeta == pytest.approx(3.0)
    assert s.voigt_modulus_frac == pytest.approx(0.15)
    assert s.voigt_viscosity_frac == pytest.approx(0.03)


def test_andrade_parameters_affect_modulus():
    """Different Andrade parameters give a different complex modulus."""
    mod = _import_rheology()
    a1 = mod.Andrade(0.3, 1.0).calc_complex_modulus(_MU, _ETA, _OMEGA)
    a2 = mod.Andrade(0.5, 1.0).calc_complex_modulus(_MU, _ETA, _OMEGA)
    assert a1 != a2
    expected = 1.0 / _ref_compliance_andrade(_ETA, _MU, _OMEGA, alpha=0.5, zeta=1.0)
    assert a2.real == pytest.approx(expected.real, rel=1e-12)
    assert a2.imag == pytest.approx(expected.imag, rel=1e-12)


# =====================================================================================================================
# Factory
# =====================================================================================================================
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
])
def test_make_rheology_aliases(alias, canonical):
    """The factory resolves aliases (case-insensitive) to canonical models."""
    mod = _import_rheology()
    inst = mod.make_rheology(alias)
    assert inst.model_name == canonical


def test_make_rheology_with_config():
    """The factory forwards config parameters to the model."""
    mod = _import_rheology()
    s = mod.make_rheology("sundberg", {
        "alpha": 0.45, "zeta": 2.5,
        "voigt_modulus_frac": 0.18, "voigt_viscosity_frac": 0.04,
    })
    assert s.alpha == pytest.approx(0.45)
    assert s.zeta == pytest.approx(2.5)
    assert s.voigt_modulus_frac == pytest.approx(0.18)
    assert s.voigt_viscosity_frac == pytest.approx(0.04)


def test_make_rheology_unknown_raises():
    """An unknown model name raises ValueError (from the C++ enum mapping)."""
    mod = _import_rheology()
    with pytest.raises(ValueError):
        mod.make_rheology("not_a_real_model")


@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_make_rheology_returns_rich_subclass(name, cls):
    """The factory returns the correct concrete Python subclass (not the base)."""
    mod = _import_rheology()
    inst = mod.make_rheology(name)
    assert type(inst).__name__ == cls
    assert isinstance(inst, getattr(mod, cls))
    assert isinstance(inst, mod.RheologyBase)


def test_make_rheology_adopted_object_is_usable():
    """An object built by the C++ enum factory is fully usable and outlives the call."""
    mod = _import_rheology()
    a = mod.make_rheology("andrade", {"alpha": 0.42, "zeta": 1.7})
    # Parameters survived the unique_ptr adoption.
    assert a.alpha == pytest.approx(0.42)
    assert a.zeta == pytest.approx(1.7)
    # Calculation dispatches to the correct concrete model.
    got = a.calc_complex_modulus(_MU, _ETA, _OMEGA)
    expected = 1.0 / _ref_compliance_andrade(_ETA, _MU, _OMEGA, alpha=0.42, zeta=1.7)
    assert got.real == pytest.approx(expected.real, rel=1e-12)
    assert got.imag == pytest.approx(expected.imag, rel=1e-12)
    # Binary round-trip still works on a factory-built object.
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "adopted.tpyb")
        a.save_binary(path)
        restored = mod.Andrade()
        restored.load_binary(path)
    assert restored.get_config_dict() == a.get_config_dict()


# =====================================================================================================================
# Config dict
# =====================================================================================================================
def test_config_dict_keys():
    """Each model's config dict carries the expected keys."""
    mod = _import_rheology()
    assert set(mod.Elastic().get_config_dict()) == {"model_name"}
    assert set(mod.Maxwell().get_config_dict()) == {"model_name"}
    assert set(mod.Voigt().get_config_dict()) == {
        "model_name", "voigt_modulus_frac", "voigt_viscosity_frac"}
    assert set(mod.Andrade().get_config_dict()) == {"model_name", "alpha", "zeta"}
    assert set(mod.Sundberg().get_config_dict()) == {
        "model_name", "alpha", "zeta", "voigt_modulus_frac", "voigt_viscosity_frac"}


def test_save_config_writes_toml():
    """save_config writes a TOML file that round-trips through toml.load."""
    import toml
    mod = _import_rheology()
    s = mod.Sundberg(0.4, 2.0, 0.15, 0.03)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "rheo.toml")
        s.save_config(path)
        loaded = toml.load(path)
    assert loaded["model_name"] == "sundberg"
    assert loaded["alpha"] == pytest.approx(0.4)
    assert loaded["voigt_modulus_frac"] == pytest.approx(0.15)


# =====================================================================================================================
# Vectorized methods and direct convenience functions
# =====================================================================================================================
import numpy as np  # noqa: E402


@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_convenience_scalar_matches_class(name, cls):
    """The lower-case convenience function matches the class for scalar input."""
    mod = _import_rheology()
    func = getattr(mod, name)            # e.g. mod.maxwell
    inst = _make(name)                   # corresponding class instance, default params
    got = func(_OMEGA, _MU, _ETA)        # convenience order: (frequency, modulus, viscosity)
    expected = inst.calc_complex_modulus(_MU, _ETA, _OMEGA)
    assert isinstance(got, complex)
    assert got == pytest.approx(expected)


@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_convenience_vectorize_modulus(name, cls):
    """Array modulus+viscosity at constant frequency returns matching ndarray."""
    mod = _import_rheology()
    func = getattr(mod, name)
    inst = _make(name)
    mods  = np.array([1.0e10, 5.0e10, 1.0e11])
    viscs = np.array([1.0e19, 1.0e20, 1.0e21])
    got = func(_OMEGA, mods, viscs)
    assert isinstance(got, np.ndarray)
    assert got.shape == (3,)
    assert got.dtype == np.complex128
    for i in range(3):
        assert got[i] == pytest.approx(
            inst.calc_complex_modulus(mods[i], viscs[i], _OMEGA))


def test_convenience_vectorize_frequency():
    """Array frequency at constant modulus/viscosity returns matching ndarray."""
    mod = _import_rheology()
    inst = mod.Andrade(0.3, 1.0)
    freqs = np.array([1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4])
    got = mod.andrade(freqs, _MU, _ETA, alpha=0.3, zeta=1.0)
    assert got.shape == (4,)
    for i in range(4):
        assert got[i] == pytest.approx(
            inst.calc_complex_modulus(_MU, _ETA, freqs[i]))


def test_convenience_vectorize_all_and_broadcast():
    """All-array and mixed (broadcast) inputs return the broadcast shape."""
    mod = _import_rheology()
    freqs = np.array([1.0e-6, 1.0e-5, 1.0e-4])
    mods  = np.array([1.0e10, 5.0e10, 1.0e11])
    viscs = np.array([1.0e19, 1.0e20, 1.0e21])
    got_all = mod.maxwell(freqs, mods, viscs)
    assert got_all.shape == (3,)
    # Mixed: frequency array, modulus array, viscosity scalar -> broadcast.
    got_mixed = mod.maxwell(freqs, mods, _ETA)
    assert got_mixed.shape == (3,)
    m = mod.Maxwell()
    for i in range(3):
        assert got_mixed[i] == pytest.approx(
            m.calc_complex_modulus(mods[i], _ETA, freqs[i]))


def test_convenience_preserves_2d_shape():
    """A 2-D input shape is preserved in the output array."""
    mod = _import_rheology()
    mods  = np.array([[1.0e10, 2.0e10], [3.0e10, 4.0e10]])
    viscs = np.full((2, 2), 1.0e20)
    got = mod.maxwell(_OMEGA, mods, viscs)
    assert got.shape == (2, 2)


def test_convenience_params_forwarded():
    """Model parameters reach the stack-allocated C++ model."""
    mod = _import_rheology()
    inst = mod.Andrade(0.42, 1.7)
    got = mod.andrade(_OMEGA, _MU, _ETA, alpha=0.42, zeta=1.7)
    assert got == pytest.approx(inst.calc_complex_modulus(_MU, _ETA, _OMEGA))


def test_class_vectorize_methods_match_scalar():
    """The class-level vectorized methods match element-wise scalar calls."""
    mod = _import_rheology()
    s = mod.Sundberg(0.3, 1.0, 0.2, 0.02)
    mods  = np.array([1.0e10, 5.0e10, 1.0e11])
    viscs = np.array([1.0e19, 1.0e20, 1.0e21])
    freqs = np.array([1.0e-6, 1.0e-5, 1.0e-4])

    by_mod = s.calc_complex_modulus_vectorize_modulus(mods, viscs, _OMEGA)
    by_frq = s.calc_complex_modulus_vectorize_frequency(_MU, _ETA, freqs)
    by_all = s.calc_complex_modulus_vectorize_all(mods, viscs, freqs)
    for i in range(3):
        assert by_mod[i] == pytest.approx(s.calc_complex_modulus(mods[i], viscs[i], _OMEGA))
        assert by_frq[i] == pytest.approx(s.calc_complex_modulus(_MU, _ETA, freqs[i]))
        assert by_all[i] == pytest.approx(s.calc_complex_modulus(mods[i], viscs[i], freqs[i]))


def test_class_vectorize_size_mismatch_raises():
    """Mismatched input lengths raise ValueError from the C++ layer."""
    mod = _import_rheology()
    m = mod.Maxwell()
    with pytest.raises(ValueError):
        m.calc_complex_modulus_vectorize_all(np.ones(3), np.ones(2), np.ones(3))
    with pytest.raises(ValueError):
        m.calc_complex_modulus_vectorize_modulus(np.ones(3), np.ones(2), _OMEGA)


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_binary_round_trip(name, cls):
    """save_binary / load_binary preserves model name and parameters."""
    mod = _import_rheology()
    cls_obj = getattr(mod, cls)
    # Build with non-default params where applicable.
    if name in ("voigt", "burgers"):
        original = cls_obj(0.33, 0.07)
    elif name == "andrade":
        original = cls_obj(0.42, 1.7)
    elif name == "sundberg":
        original = cls_obj(0.42, 1.7, 0.33, 0.07)
    else:
        original = cls_obj()

    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, f"{name}.tpyb")
        original.save_binary(path)
        restored = cls_obj()
        restored.load_binary(path)

    assert restored.model_name == name
    assert restored.get_config_dict() == original.get_config_dict()


def test_load_binary_missing_file_raises():
    """Loading a non-existent binary raises FileNotFoundError."""
    mod = _import_rheology()
    with pytest.raises(FileNotFoundError):
        mod.Maxwell().load_binary("does_not_exist_12345.tpyb")


# =====================================================================================================================
# isinstance chain
# =====================================================================================================================
@pytest.mark.parametrize("name,cls", _MODEL_CLASSES)
def test_isinstance_chain(name, cls):
    """Every model is a RheologyBase, PhysicsBase and TidalPyBaseClass."""
    mod = _import_rheology()
    from TidalPy.Utilities_x.classes_x.classes import PhysicsBase, TidalPyBaseClass
    inst = _make(name)
    assert isinstance(inst, mod.RheologyBase)
    assert isinstance(inst, PhysicsBase)
    assert isinstance(inst, TidalPyBaseClass)
