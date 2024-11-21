import gcmotion as gcm
import numpy as np

# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

B0 = Q(B0num, "Tesla")
i = Q(0, "NUPlasma_current")
g = Q(1, "NUPlasma_current")

# ================================ Q_FACTORS================================

bfields = [gcm.bfield.LAR(B0=B0, i=i, g=g)]

# ================================== TESTS =================================


def test_bigNU_return_types():
    """Test that the bigNU() I/O is as annotated:
    (phi: float | ndarray, theta: float | ndarray) → float | ndarray."""

    for bfield in bfields:
        # floats
        for psi in np.linspace(0.01, 2, 10):
            for theta in np.linspace(0, 2 * np.pi, 10):
                b, i, g = bfield.bigNU(psi, theta)
                assert isinstance(b, (int, float))
                assert isinstance(g, (int, float))
                assert isinstance(g, (int, float))

        # 1d arrays
        psi = np.random.rand(100)
        theta = np.random.rand(100)
        b, i, g = bfield.bigNU(psi, theta)
        assert isinstance(b, np.ndarray)
        assert isinstance(i, np.ndarray)
        assert isinstance(g, np.ndarray)
        assert psi.shape == b.shape
        assert psi.shape == i.shape
        assert psi.shape == g.shape

        # 2d arrays
        psi = np.random.rand(100, 100)
        theta = np.random.rand(100, 100)
        b, i, g = bfield.bigNU(psi, theta)
        assert isinstance(b, np.ndarray)
        assert isinstance(i, np.ndarray)
        assert isinstance(g, np.ndarray)
        assert psi.shape == b.shape
        assert psi.shape == i.shape
        assert psi.shape == g.shape


def test_solverbNU_return_types():
    """Test that the solverbNU() I/O is as annotated:
    (psi: float, theta: float) → tuple[float, float, float]."""

    for bfield in bfields:
        for psi in np.linspace(0.01, 2, 10):
            for theta in np.linspace(0, 2 * np.pi, 10):
                b, b_der, currents, currents_der = bfield.solverbNU(psi, theta)
                b_der_psi, b_der_theta = b_der
                i, g = currents
                i_der, g_der = currents_der
                assert isinstance(b, (int, float))
                assert isinstance(b_der_psi, (int, float))
                assert isinstance(b_der_theta, (int, float))
                assert isinstance(i, (int, float))
                assert isinstance(g, (int, float))
                assert isinstance(i_der, (int, float))
                assert isinstance(g_der, (int, float))
