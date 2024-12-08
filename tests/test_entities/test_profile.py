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

# Construct a Tokamak
tokamak = gcm.Tokamak(
    R=R,
    a=a,
    qfactor=gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2),
    bfield=gcm.bfield.LAR(B0=B0, i=i, g=g),
    efield=gcm.efield.Nofield(),
    # efield=gcm.efield.Radial(a, Ea, B0, peak=0.9, rw=1 / 10),
)

# Construct a PhysicalParameters set
paramsE = gcm.PhysicalParameters(
    species=species,
    mu=Q(1e-5, "NUMagnetic_moment"),
    Pzeta=Q(-0.015, "NUMagnetic_flux"),
)
paramsPzeta = gcm.PhysicalParameters(
    species=species,
    mu=Q(1e-5, "NUMagnetic_moment"),
    E=Q(4, "keV"),
)
paramsmu = gcm.PhysicalParameters(
    species=species,
    Pzeta=Q(-0.015, "NUMagnetic_flux"),
    E=Q(4, "keV"),
)

# Create a Profile
profileE = gcm.Profile(tokamak, paramsE)
profilePzeta = gcm.Profile(tokamak, paramsPzeta)
profilemu = gcm.Profile(tokamak, paramsmu)

# =========================== TESTS ================================


class TestProfileInitialization:

    def test_attribute_inheritance(self):
        """Test that all tokamak's and params' attributes are also profile's
        attributes (by name only)."""

        tokamakvars = set(vars(tokamak))
        paramsvars = set(vars(paramsE))
        profilevars = set(vars(profileE))

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

        assert profileE.findPtheta(psi).units == psi.units
        assert profileE.findPtheta(psiNU).units == psiNU.units
        assert profileE.findPtheta(psiT).units == psiT.units
        assert profileE.findPtheta(psi_wall).units == psi_wall.units
        with pytest.raises(AttributeError):
            profileE.findPtheta(psi_float)

    def test_return_types(self):

        psi = Q(0.5, "NUMagnetic_flux")
        psilist = Q([0, 1], "psi_wall")
        psiarray = Q(np.linspace(0, 5, 10), "NUMagnetic_flux")
        psigrid = Q(
            np.meshgrid(np.linspace(0, 5, 10), np.linspace(0, 5, 10)),
            "Magnetic_flux",
        )

        Ptheta = profileE.findPtheta(psi)
        Pthetalist = profileE.findPtheta(psilist)
        Pthetaarray = profileE.findPtheta(psiarray)
        Pthetagrid = profileE.findPtheta(psigrid)

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

        assert profileE.findEnergy(psi, theta, "keV").units == ureg.keV
        assert profileE.findEnergy(psiNU, theta, "keV").units == ureg.keV
        assert profileE.findEnergy(psiT, theta, "keV").units == ureg.keV
        assert profileE.findEnergy(psi_wall, theta, "keV").units == ureg.keV

        assert profileE.findEnergy(psi, theta, "keV").units == ureg.keV
        assert profileE.findEnergy(psi, theta, "eV").units == ureg.eV
        assert profileE.findEnergy(psi, theta, "Joule").units == ureg.Joule
        assert profileE.findEnergy(psi, theta, "NUJoule").units == ureg.NUJoule

    def test_return_types(self):

        theta = 0
        psi = Q(0.5, "NUMagnetic_flux")
        psilist = Q([0, 1], "psi_wall")
        psiarray = Q(np.linspace(0, 5, 10), "NUMagnetic_flux")
        psigrid = Q(
            np.meshgrid(np.linspace(0, 5, 10), np.linspace(0, 5, 10)),
            "Magnetic_flux",
        )

        Energy = profileE.findEnergy(psi, theta, "keV")
        Energylist = profileE.findEnergy(psilist, theta, "keV")
        Energyarray = profileE.findEnergy(psiarray, theta, "keV")
        Energygrid = profileE.findEnergy(psigrid, theta, "keV")

        assert isinstance(Energy, pint.Quantity)
        assert isinstance(Energylist, pint.Quantity)
        assert isinstance(Energyarray, pint.Quantity)
        assert isinstance(Energygrid, pint.Quantity)

        assert isinstance(Energy.m, float)
        assert len(Energylist) == len(psilist)
        assert Energyarray.shape == psiarray.shape
        assert Energygrid.shape == psigrid.shape


class TestProfileFindPzeta:

    def test_units(self):

        theta = 0
        psi = Q(0.5, "psi_wall").to("Magnetic_flux")
        psiNU = Q(0.5, "psi_wall").to("NUMagnetic_flux")
        psiT = Q(0.5, "psi_wall").to("Tesla*meter^2")
        psi_wall = Q(0.5, "psi_wall")

        assert (
            profilePzeta.findPzeta(psi, theta, "Magnetic_flux").units
            == ureg.Magnetic_flux
        )
        assert (
            profilePzeta.findPzeta(psiNU, theta, "Magnetic_flux").units
            == ureg.Magnetic_flux
        )
        assert (
            profilePzeta.findPzeta(psiT, theta, "Magnetic_flux").units
            == ureg.Magnetic_flux
        )
        assert (
            profilePzeta.findPzeta(psi_wall, theta, "Magnetic_flux").units
            == ureg.Magnetic_flux
        )

        assert (
            profilePzeta.findPzeta(psi, theta, "Magnetic_flux").units
            == ureg.Magnetic_flux
        )
        assert (
            profilePzeta.findPzeta(psi, theta, "NUMagnetic_flux").units
            == ureg.NUMagnetic_flux
        )
        assert (
            profilePzeta.findPzeta(psi, theta, "psi_wall").units
            == ureg.psi_wall
        )
        assert (
            profilePzeta.findPzeta(psi, theta, "NUpsi_wall").units
            == ureg.NUpsi_wall
        )

    def test_return_types(self):

        theta = 0
        psi = Q(0.5, "psi_wall").to("NUMagnetic_flux")
        psilist = Q([0, 0.1], "psi_wall")
        psiarray = Q(np.linspace(0, 0.5, 10), "psi_wall").to("NUMagnetic_flux")
        psigrid = Q(
            np.meshgrid(np.linspace(0, 0.5, 10), np.linspace(0, 0.5, 10)),
            "psi_wall",
        ).to("Magnetic_flux")

        Pzeta = profilePzeta.findPzeta(psi, theta, "Magnetic_flux")
        Pzetalist = profilePzeta.findPzeta(psilist, theta, "Magnetic_flux")
        Pzetaarray = profilePzeta.findPzeta(psiarray, theta, "Magnetic_flux")
        Pzetagrid = profilePzeta.findPzeta(psigrid, theta, "Magnetic_flux")

        assert isinstance(Pzeta, pint.Quantity)
        assert isinstance(Pzetalist, pint.Quantity)
        assert isinstance(Pzetaarray, pint.Quantity)
        assert isinstance(Pzetagrid, pint.Quantity)

        assert isinstance(Pzeta.m, float)
        assert len(Pzetalist) == len(psilist)
        assert Pzetaarray.shape == psiarray.shape
        assert Pzetagrid.shape == psigrid.shape


class TestProfileFindmu:

    def test_units(self):

        theta = 0
        psi = Q(0.5, "psi_wall").to("Magnetic_flux")
        psiNU = Q(0.5, "psi_wall").to("NUMagnetic_flux")
        psiT = Q(0.5, "psi_wall").to("Tesla*meter^2")
        psi_wall = Q(0.5, "psi_wall")

        assert (
            profilemu.findmu(psi, theta, "Magnetic_moment").units
            == ureg.Magnetic_moment
        )
        assert (
            profilemu.findmu(psiNU, theta, "Magnetic_moment").units
            == ureg.Magnetic_moment
        )
        assert (
            profilemu.findmu(psiT, theta, "Magnetic_moment").units
            == ureg.Magnetic_moment
        )
        assert (
            profilemu.findmu(psi_wall, theta, "Magnetic_moment").units
            == ureg.Magnetic_moment
        )

        assert (
            profilemu.findmu(psi, theta, "Magnetic_moment").units
            == ureg.Magnetic_moment
        )
        assert (
            profilemu.findmu(psi, theta, "Magnetic_moment").units
            == ureg.Magnetic_moment
        )
        assert (
            profilemu.findmu(psi, theta, "NUMagnetic_moment").units
            == ureg.NUMagnetic_moment
        )

    def test_return_types(self):

        theta = 0
        psi = Q(0.5, "psi_wall").to("NUMagnetic_flux")
        psilist = Q([0, 0.1], "psi_wall")
        psiarray = Q(np.linspace(0, 0.5, 10), "psi_wall").to("NUMagnetic_flux")
        psigrid = Q(
            np.meshgrid(np.linspace(0, 0.5, 10), np.linspace(0, 0.5, 10)),
            "psi_wall",
        ).to("Magnetic_flux")

        mu = profilemu.findmu(psi, theta, "Magnetic_moment")
        mulist = profilemu.findmu(psilist, theta, "Magnetic_moment")
        muarray = profilemu.findmu(psiarray, theta, "Magnetic_moment")
        mugrid = profilemu.findmu(psigrid, theta, "Magnetic_moment")

        assert isinstance(mu, pint.Quantity)
        assert isinstance(mulist, pint.Quantity)
        assert isinstance(muarray, pint.Quantity)
        assert isinstance(mugrid, pint.Quantity)

        assert isinstance(mu.m, float)
        assert len(mulist) == len(psilist)
        assert muarray.shape == psiarray.shape
        assert mugrid.shape == psigrid.shape
