"""
Tests for TidalPy.Tides_x.love.love — LoveNumbers (Phase 2).

Covers construction, property access, equality, iteration/tuple unpacking,
to_dict, repr, and integration with PhysicsLayer.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import pytest


# =====================================================================================================================
# Helpers
# =====================================================================================================================

def _import_love():
    try:
        from TidalPy.Tides_x.love import love as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.Tides_x.love.love not compiled — run uv pip install first."
        )


def _make(k=0+0j, h=0+0j, l=0+0j):
    mod = _import_love()
    return mod.LoveNumbers(k=k, h=h, l=l)


# =====================================================================================================================
# Construction and property access
# =====================================================================================================================
def test_love_numbers_defaults_zero():
    """LoveNumbers defaults to k=h=l=0+0j."""
    ln = _make()
    assert ln.k == pytest.approx(0.0 + 0.0j)
    assert ln.h == pytest.approx(0.0 + 0.0j)
    assert ln.l == pytest.approx(0.0 + 0.0j)


def test_love_numbers_real_values():
    """LoveNumbers stores real-valued k, h, l correctly."""
    ln = _make(k=0.3, h=0.6, l=0.1)
    assert ln.k == pytest.approx(0.3 + 0.0j)
    assert ln.h == pytest.approx(0.6 + 0.0j)
    assert ln.l == pytest.approx(0.1 + 0.0j)


def test_love_numbers_complex_values():
    """LoveNumbers stores fully complex k, h, l correctly."""
    ln = _make(k=0.3 - 0.01j, h=0.6 - 0.02j, l=0.1 - 0.005j)
    assert ln.k == pytest.approx(0.3  - 0.01j)
    assert ln.h == pytest.approx(0.6  - 0.02j)
    assert ln.l == pytest.approx(0.1  - 0.005j)


def test_love_numbers_properties_are_complex():
    """k, h, l properties return Python complex objects."""
    ln = _make(k=1.0 + 2.0j, h=3.0 + 4.0j, l=5.0 + 6.0j)
    assert isinstance(ln.k, complex)
    assert isinstance(ln.h, complex)
    assert isinstance(ln.l, complex)


# =====================================================================================================================
# Equality
# =====================================================================================================================
def test_love_numbers_equal_to_itself():
    """A LoveNumbers is equal to itself."""
    ln = _make(k=0.3 - 0.01j, h=0.6 - 0.02j, l=0.1 - 0.005j)
    assert ln == ln


def test_love_numbers_equal_to_copy():
    """Two LoveNumbers with identical values are equal."""
    ln1 = _make(k=0.3 - 0.01j, h=0.6 - 0.02j, l=0.1 - 0.005j)
    ln2 = _make(k=0.3 - 0.01j, h=0.6 - 0.02j, l=0.1 - 0.005j)
    assert ln1 == ln2


def test_love_numbers_not_equal_when_different():
    """LoveNumbers with different values are not equal."""
    ln1 = _make(k=0.3, h=0.6, l=0.1)
    ln2 = _make(k=0.4, h=0.6, l=0.1)
    assert ln1 != ln2


def test_love_numbers_eq_not_implemented_for_non_love():
    """Equality against a non-LoveNumbers returns NotImplemented."""
    ln = _make()
    result = ln.__eq__(42)
    assert result is NotImplemented


# =====================================================================================================================
# Iteration / tuple unpacking
# =====================================================================================================================
def test_love_numbers_iter():
    """Iterating LoveNumbers yields (k, h, l) in order."""
    ln   = _make(k=0.3 - 0.01j, h=0.6 - 0.02j, l=0.1 - 0.005j)
    vals = list(ln)
    assert len(vals) == 3
    assert vals[0] == pytest.approx(0.3  - 0.01j)
    assert vals[1] == pytest.approx(0.6  - 0.02j)
    assert vals[2] == pytest.approx(0.1  - 0.005j)


def test_love_numbers_tuple_unpack():
    """LoveNumbers can be unpacked as (k, h, l)."""
    ln = _make(k=1.0j, h=2.0j, l=3.0j)
    k, h, l = ln
    assert k == pytest.approx(1.0j)
    assert h == pytest.approx(2.0j)
    assert l == pytest.approx(3.0j)


# =====================================================================================================================
# to_dict
# =====================================================================================================================
def test_love_numbers_to_dict_keys():
    """to_dict returns all six expected keys."""
    ln  = _make()
    d   = ln.to_dict()
    for key in ("love_number_k_re", "love_number_k_im",
                "love_number_h_re", "love_number_h_im",
                "love_number_l_re", "love_number_l_im"):
        assert key in d, f"Missing key: {key}"


def test_love_numbers_to_dict_values():
    """to_dict values match the stored love numbers."""
    ln = _make(k=0.3 - 0.01j, h=0.6 - 0.02j, l=0.1 - 0.005j)
    d  = ln.to_dict()
    assert d["love_number_k_re"] == pytest.approx(0.3)
    assert d["love_number_k_im"] == pytest.approx(-0.01)
    assert d["love_number_h_re"] == pytest.approx(0.6)
    assert d["love_number_h_im"] == pytest.approx(-0.02)
    assert d["love_number_l_re"] == pytest.approx(0.1)
    assert d["love_number_l_im"] == pytest.approx(-0.005)


def test_love_numbers_to_dict_zeros():
    """to_dict returns zeros for a default-constructed LoveNumbers."""
    d = _make().to_dict()
    for key in d:
        assert d[key] == pytest.approx(0.0), f"{key} should be 0.0"


# =====================================================================================================================
# repr
# =====================================================================================================================
def test_love_numbers_repr():
    """repr includes k, h, l values."""
    ln = _make(k=0.3 - 0.01j, h=0.6, l=0.0)
    r  = repr(ln)
    assert "LoveNumbers" in r
    assert "k=" in r
    assert "h=" in r
    assert "l=" in r


# =====================================================================================================================
# Integration: LoveNumbers from PhysicsLayer
# =====================================================================================================================
def test_love_numbers_from_physics_layer():
    """PhysicsLayer.love_numbers returns a LoveNumbers with correct values."""
    try:
        from TidalPy.structures_x.layers.physics import PhysicsLayer
        from TidalPy.Tides_x.love.love import LoveNumbers
    except ImportError:
        pytest.skip("PhysicsLayer or LoveNumbers not compiled.")

    pl = PhysicsLayer("test", 0, 0.0, 1e6, 1e20,
                      love_number_k=0.3 - 0.01j,
                      love_number_h=0.6 - 0.02j,
                      love_number_l=0.1 - 0.005j)
    ln = pl.love_numbers
    assert isinstance(ln, LoveNumbers)
    assert ln.k == pytest.approx(0.3  - 0.01j)
    assert ln.h == pytest.approx(0.6  - 0.02j)
    assert ln.l == pytest.approx(0.1  - 0.005j)


def test_love_numbers_from_physics_layer_to_dict():
    """PhysicsLayer.love_numbers.to_dict matches get_config_dict love entries."""
    try:
        from TidalPy.structures_x.layers.physics import PhysicsLayer
    except ImportError:
        pytest.skip("PhysicsLayer not compiled.")

    pl  = PhysicsLayer("test", 0, 0.0, 1e6, 1e20,
                       love_number_k=0.3 - 0.01j,
                       love_number_h=0.6 - 0.02j,
                       love_number_l=0.1 - 0.005j)
    ld  = pl.love_numbers.to_dict()
    cfg = pl.get_config_dict()
    for key in ("love_number_k_re", "love_number_k_im",
                "love_number_h_re", "love_number_h_im",
                "love_number_l_re", "love_number_l_im"):
        assert ld[key] == pytest.approx(cfg[key]), f"Mismatch for {key}"
