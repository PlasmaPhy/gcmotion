import gcmotion as gcm
import pint
import numpy as np
from pint import get_application_registry
import pytest


Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
Q = gcm.QuantityConstructor(R=Rnum, a=anum, B0=B0num, species=species)
ureg = get_application_registry()


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
    efield=gcm.efield.Nofield(),
)

params = gcm.PhysicalParameters(
    species=species,
    mu=Q(1e-5, "NUMagnetic_moment"),
    Pzeta0=Q(-0.027, "NUmagnetic_flux"),
)

profile = gcm.Profile(tokamak, params)


# =========================== TESTS ================================


class TestProfileInitialization:

    def test_attribute_inheritance(self):
        """Test that all tokamak's and params' attributes are also profile's
        attributes (by name only)."""

        tokamakvars = set(vars(tokamak))
        paramsvars = set(vars(params))
        profilevars = set(vars(profile))

        assert profilevars.issuperset(tokamakvars)
        assert profilevars.issuperset(paramsvars)

    def test_str_repr_functionality(self):
        """Tests that __repr__() and __str__() can return with no errors"""
        _ = tokamak.__repr__()
        _ = tokamak.__str__()


class TestProfileFindPtheta:

    def test_units(self):

        psi = Q(0.5, "Magnetic_flux")
        psiNU = Q(0.5, "NUMagnetic_flux")
        psiT = Q(0.5, "Tesla*meter^2")
        psi_wall = Q(0.5, "psi_wall")
        psi_float = 0.5

        assert profile.findPtheta(psi).units == psi.units
        assert profile.findPtheta(psiNU).units == psiNU.units
        assert profile.findPtheta(psiT).units == psiT.units
        assert profile.findPtheta(psi_wall).units == psi_wall.units
        with pytest.raises(AttributeError):
            profile.findPtheta(psi_float)

    def test_return_types(self):

        psi = Q(0.5, "NUMagnetic_flux")
        psilist = Q([0, 1], "psi_wall")
        psiarray = Q(np.linspace(0, 5, 10), "NUMagnetic_flux")
        psigrid = Q(
            np.meshgrid(np.linspace(0, 5, 10), np.linspace(0, 5, 10)),
            "Magnetic_flux"
        )

        Ptheta = profile.findPtheta(psi)
        Pthetalist = profile.findPtheta(psilist)
        Pthetaarray = profile.findPtheta(psiarray)
        Pthetagrid = profile.findPtheta(psigrid)

        assert isinstance(Ptheta, pint.Quantity)
        assert isinstance(Pthetalist, pint.Quantity)
        assert isinstance(Pthetaarray, pint.Quantity)
        assert isinstance(Pthetagrid, pint.Quantity)

        assert isinstance(Ptheta.m, float)
        assert len(Pthetalist) == len(psilist)
        assert Pthetaarray.shape == psiarray.shape
        assert Pthetagrid.shape == psigrid.shape


class TestProfileFindEnergy:

    def test_units(self):

        theta = 0
        psi = Q(0.5, "Magnetic_flux")
        psiNU = Q(0.5, "NUMagnetic_flux")
        psiT = Q(0.5, "Tesla*meter^2")
        psi_wall = Q(0.5, "psi_wall")

        assert profile.findEnergy(psi, theta, "keV").units == ureg.keV
        assert profile.findEnergy(psiNU, theta, "keV").units == ureg.keV
        assert profile.findEnergy(psiT, theta, "keV").units == ureg.keV
        assert profile.findEnergy(psi_wall, theta, "keV").units == ureg.keV

        assert profile.findEnergy(psi, theta, "keV").units == ureg.keV
        assert profile.findEnergy(psi, theta, "eV").units == ureg.eV
        assert profile.findEnergy(psi, theta, "Joule").units == ureg.Joule
        assert profile.findEnergy(psi, theta, "NUJoule").units == ureg.NUJoule

    def test_return_types(self):

        theta = 0
        psi = Q(0.5, "NUMagnetic_flux")
        psilist = Q([0, 1], "psi_wall")
        psiarray = Q(np.linspace(0, 5, 10), "NUMagnetic_flux")
        psigrid = Q(
            np.meshgrid(np.linspace(0, 5, 10), np.linspace(0, 5, 10)),
            "Magnetic_flux"
        )

        Energy = profile.findEnergy(psi, theta, "keV")
        Energylist = profile.findEnergy(psilist, theta, "keV")
        Energyarray = profile.findEnergy(psiarray, theta, "keV")
        Energygrid = profile.findEnergy(psigrid, theta, "keV")

        assert isinstance(Energy, pint.Quantity)
        assert isinstance(Energylist, pint.Quantity)
        assert isinstance(Energyarray, pint.Quantity)
        assert isinstance(Energygrid, pint.Quantity)

        assert isinstance(Energy.m, float)
        assert len(Energylist) == len(psilist)
        assert Energyarray.shape == psiarray.shape
        assert Energygrid.shape == psigrid.shape
