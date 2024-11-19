import gcmotion as gcm
import numpy as np

# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

a = Q(anum, "meters")
B0 = Q(B0num, "Tesla")
Ea = Q(73500, "Volts/meter")

# ================================ EFIELDS =================================

efields = [
    gcm.efield.Nofield(),
]

efields += [
    gcm.efield.Radial(a, Ea, B0, peak=peak, rw=rw)
    for peak in np.linspace(0.1, 1, 10)
    for rw in np.linspace(0.01, 0.5, 10)
]

# ================================== TESTS =================================


def test_PhiNU_return_types():
    """Test that the PhiNU() I/O is as annotated:
    (psi: float | ndarray, theta: float | ndarray) → float | ndarray"""

    for efield in efields:
        # floats
        for psi in np.linspace(0.01, 2, 10):
            for theta in np.linspace(0, 2 * np.pi, 10):
                Phi = efield.PhiNU(psi, theta)
                assert isinstance(Phi, (int, float))

        # 1d arrays
        psi = np.random.rand(100)
        theta = np.random.rand(100)
        Phi = efield.PhiNU(psi, theta)
        assert isinstance(Phi, np.ndarray)
        assert psi.shape == Phi.shape

        # 2d arrays
        psi = np.random.rand(100, 100)
        theta = np.random.rand(100, 100)
        Phi = efield.PhiNU(psi, theta)
        assert isinstance(Phi, np.ndarray)
        assert psi.shape == Phi.shape


def test_Er_return_types():
    """Test that the Er() I/O is as annotated:
    (psi: ndarray) → ndarray."""

    for efield in efields:
        # 1d arrays
        psi = np.random.rand(100)
        Er = efield.Er(psi)
        assert isinstance(Er, np.ndarray)
        assert psi.shape == Er.shape

        # 2d arrays
        psi = np.random.rand(100, 100)
        Er = efield.Er(psi)
        assert isinstance(Er, np.ndarray)
        assert psi.shape == Er.shape
