import pytest
import gcmotion as gcm

# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

# ========================== BASE OBJECTS SETUP ===============================

params = gcm.PhysicalParameters(
    species=species,
    mu=Q(1e-5, "NUMagnetic_moment"),
    Pzeta=Q(-0.025, "NUMagnetic_flux"),
)


def test_units():
    """Check that all params attributes have the correct units"""

    # Input Parameters
    assert isinstance(params.species, str)
    assert params.mu.units == ureg.Magnetic_moment
    assert params.Pzeta.units == ureg.Magnetic_flux

    # Quantities
    assert isinstance(params.species_name, str)
    assert params.muNU.units == ureg.NUMagnetic_moment
    assert params.PzetaNU.units == ureg.NUMagnetic_flux


def test_str_repr_functionality():
    """Tests that __repr__() and __str__() can return with no errors"""
    _ = params.__repr__()
    _ = params.__str__()
