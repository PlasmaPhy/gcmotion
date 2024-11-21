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

# ================================ Q_FACTORS================================

qfactors = [
    gcm.qfactor.Unity(),
]

qfactors += [
    gcm.qfactor.Parabolic(a, B0, q0=q0, q_wall=q_wall)
    for q0 in np.linspace(1, 4, 5)
    for q_wall in np.linspace(q0 + 0.1, 5, 5)  # q_wall =/= q0
]

qfactors += [
    gcm.qfactor.Hypergeometric(a, B0, q0=q0, q_wall=q_wall, n=n)
    for q0 in np.linspace(1, 4, 5)
    for q_wall in np.linspace(q0, 5, 5)
    for n in [1, 2, 4]
]

# ================================== TESTS =================================


def test_solverqNU_return_types():
    """Test that the solverqNU() I/O is as annotated (float -> float)."""

    for qfactor in qfactors:
        for psi in np.linspace(0.01, 2, 100):
            q = qfactor.solverqNU(psi)
            assert isinstance(q, (int, float))


def test_psipNU_return_types():
    """Test that the psipNU() I/O is as annotated (float | np.ndarray -> float | np.ndarray)."""

    for qfactor in qfactors:
        # float
        for psi in np.linspace(0.01, 2, 100):
            psip = qfactor.psipNU(psi)
            assert isinstance(psip, (int, float))

        # 1d array
        psi = np.random.rand(100)
        psip = qfactor.psipNU(psi)
        assert isinstance(psip, np.ndarray)
        assert psi.shape == psip.shape

        # 2d array
        psi = np.random.rand(100, 100)
        psip = qfactor.psipNU(psi)
        assert isinstance(psip, np.ndarray)
        assert psi.shape == psip.shape


def test_nondecreasing():
    """Test that q(ψ) and ψ_p(ψ) are nondecreasing."""

    psis = np.linspace(0, 2, 1000)

    q = np.zeros(psis.shape)
    for qfactor in qfactors:

        for i in range(len(psis)):
            q[i] = qfactor.solverqNU(psis[i])
        psip = qfactor.psipNU(psis)

        assert np.all(np.diff(q) >= 0)
        assert np.all(np.diff(psip) >= 0)
