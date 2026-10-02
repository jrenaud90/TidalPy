"""Shared fixtures for the radial solver tests."""
import pytest


def _starting_conditions_supported(layer_type, is_static, is_incompressible, use_kamata):
    """Whether the shooting method has starting conditions for a start layer of this kind
    (RadialSolver/starting/driver_.hpp): Saito's for any static liquid, and otherwise every combination except an
    incompressible static solid with Kamata's and any incompressible layer with Takeuchi's."""
    if (layer_type == "liquid") and is_static:
        return True
    if use_kamata and (layer_type == "solid") and is_static and is_incompressible:
        return False
    return use_kamata or not is_incompressible


@pytest.fixture
def starting_conditions_supported():
    """The predicate :func:`_starting_conditions_supported`, so a test asserts that an unsupported combination raises
    rather than skipping on whatever raises."""
    return _starting_conditions_supported
