import numba

from scipy.special import hyp2f1

n = 2
q0 = 1.1
q_wall = 4
psi_wallNU = 0.4


def hypergeometric_psipNU(psi):

    a = b = 1 / n
    c = 1 + 1 / n
    z = (1 - (q_wall / q0) ** n) * (psi / psi_wallNU) ** n
    return psi / q0 * hyp2f1(a, b, c, z)


@numba.jit
def _psipNU(psi, n, q_wall, q0, psi_wall):

    a = b = 1 / n
    c = 1 + 1 / n
    z = (1 - (q_wall / q0) ** n) * (psi / psi_wall) ** n
    return a, b, c, z


def hypergeometric_psipNU_jit(psi):
    a, b, c, z = _psipNU(
        psi,
        n=n,
        q0=q0,
        q_wall=q_wall,
        psi_wall=psi_wallNU,
    )
    return psi / q0 * hyp2f1(a, b, c, z)


# Call it once to compile
hypergeometric_psipNU_jit(psi=0.1)


psi = 0.1


def test_hypergeometric_psipNU(benchmark):
    benchmark(hypergeometric_psipNU, psi)


def test_hypergeometric_psipNU_jit(benchmark):
    benchmark(hypergeometric_psipNU_jit, psi)
