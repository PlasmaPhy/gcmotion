import pytest
import gcmotion as gcm

# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

# ========================== BASE OBJECTS SETUP ===============================

R = Q(Rnum, "meters")
a = Q(anum, "meters")
B0 = Q(B0num, "Tesla")
i = Q(0, "NUPlasma_current")
g = Q(1, "NUPlasma_current")


tokamak = gcm.Tokamak(
    R=R,
    a=a,
    qfactor=gcm.qfactor.Unity(),
    bfield=gcm.bfield.LAR(B0=B0, i=i, g=g),
    efield=None,
)


def test_units():
    """Check that all tokamak attributes have the correct units"""

    # Input Parameters
    assert tokamak.R.units == ureg.meter
    assert tokamak.a.units == ureg.meter

    # Quantities
    assert tokamak.B0.units == ureg.Tesla
    assert tokamak.RNU.units == ureg.NUmeter
    assert tokamak.aNU.units == ureg.NUmeter
    assert tokamak.B0NU.units == ureg.NUTesla

    # Last closed surfaces
    assert tokamak.psi_wall.units == ureg.Magnetic_flux
    assert tokamak.psip_wall.units == ureg.Magnetic_flux
    assert tokamak.psi_wallNU.units == ureg.NUMagnetic_flux
    assert tokamak.psip_wallNU.units == ureg.NUMagnetic_flux


def test_str_repr_functionality():
    """Tests that __repr__() and __str__() can return with no errors"""
    _ = tokamak.__repr__()
    _ = tokamak.__str__()
