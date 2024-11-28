from pint import get_application_registry
import gcmotion as gcm
import pytest


Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
Q = gcm.QuantityConstructor(R=Rnum, a=anum, B0=B0num, species=species)
ureg = get_application_registry()


@pytest.fixture(scope="module")
def tokamak():
    R = Q(Rnum, "meters")
    a = Q(anum, "meters")
    B0 = Q(B0num, "Tesla")
    i = Q(0, "NUPlasma_current")
    g = Q(1, "NUPlasma_current")

    yield gcm.Tokamak(
        R=R,
        a=a,
        qfactor=gcm.qfactor.Unity(),
        bfield=gcm.bfield.LAR(B0=B0, i=i, g=g),
        efield=gcm.efield.Nofield(),
    )


class TestTokamakInitialization:

    @pytest.fixture(autouse=True)
    def _set_tokamak(self, tokamak):
        self.tokamak = tokamak

    def test_units(self):
        """Check that all tokamak attributes have the correct units"""

        # Input Parameters
        # assert tokamak.R.units == ureg.meter
        assert self.tokamak.a.units == ureg.meter

        # Quantities
        assert self.tokamak.B0.units == ureg.Tesla
        assert self.tokamak.RNU.units == ureg.NUmeter
        assert self.tokamak.aNU.units == ureg.NUmeter
        assert self.tokamak.B0NU.units == ureg.NUTesla

        # Last closed surfaces
        assert self.tokamak.psi_wall.units == ureg.Magnetic_flux
        assert self.tokamak.psip_wall.units == ureg.Magnetic_flux
        assert self.tokamak.psi_wallNU.units == ureg.NUMagnetic_flux
        assert self.tokamak.psip_wallNU.units == ureg.NUMagnetic_flux

    def test_str_repr_functionality(self):
        """Tests that __repr__() and __str__() can return with no errors"""
        _ = self.tokamak.__repr__()
        _ = self.tokamak.__str__()
