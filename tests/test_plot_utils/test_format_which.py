from gcmotion.plot.particle_evolution import format_which


def test_format_which():

    which1 = ""
    which2 = "theta"
    which3 = "theta, zeta"
    which4 = "zeta, theta"
    which5 = "zeta theta"
    which6 = "all"
    which7 = "Pzeta  ,   Ptheta theta rho,"

    assert format_which(which1) == []
    assert format_which(which2) == ["theta"]
    assert format_which(which3) == ["theta", "zeta"]
    assert format_which(which4) == ["theta", "zeta"]
    assert format_which(which5) == ["theta", "zeta"]
    assert format_which(which6) == ["theta", "zeta", "psi", "psip", "rho",
                                    "Ptheta", "Pzeta"]
    assert format_which(which7) == ["theta", "rho", "Ptheta", "Pzeta"]
