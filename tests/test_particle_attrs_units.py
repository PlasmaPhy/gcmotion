import gcmotion as gcm
import numpy as np

# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

# ============================= PARTICLE SETUP =============================

R = Q(Rnum * 1000, "millimeters")
a = Q(anum / 1000, "kilometers")
B0 = Q(B0num, "Tesla")
i = Q(0, "NUPlasma_current")
g = Q(1, "NUPlasma_current")
Ea = Q(73500, "Volts/meter")

tokamak = {
    "R": R,
    "a": a,
    "B0": B0,
    "qfactor": gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2),
    "bfield": gcm.bfield.LAR(B0=B0, i=i, g=g),
    "efield": gcm.efield.Radial(a, Ea, B0, peak=0.98, rw=1 / 50),
}
parameters = {
    "species": species,
    "mu/muB": Q(5, "keV"),
    "theta0": 0,
    "zeta0": 0,
    "psi0": Q(0.78, "psi_wall"),
    "Pzeta0": Q(-0.0272, "NUMagnetic_flux"),
    "t_eval": Q(np.linspace(0, 1e-3, 100000), "seconds"),
}

cwp = gcm.Particle(tokamak, parameters)
events = [gcm.events.when_theta(parameters["theta0"], 8)]
cwp.run(events=events)

# ================================== TESTS ==================================


def test_attrs_units():
    # in the order of definition
    assert cwp.mi.units == ureg.kilogram
    assert cwp.qi.units == ureg.Coulomb

    assert cwp.R.units == ureg.meters
    assert cwp.a.units == ureg.meters
    assert cwp.B0.units == ureg.Tesla
    assert cwp.psi_wall.units == ureg.Magnetic_flux
    assert cwp.psip_wall.units == ureg.Magnetic_flux

    assert cwp.theta0.units == ureg.radians
    assert cwp.zeta0.units == ureg.radians
    assert cwp.psi0.units == ureg.Magnetic_flux
    assert cwp.Pzeta0.units == ureg.Magnetic_flux
    assert cwp.t_eval.units == ureg.seconds
    assert cwp.tfinal.units == ureg.seconds
    assert cwp.mu.units == ureg.keV / ureg.Tesla
    assert cwp.muB.units == ureg.keV
    assert cwp.psip0.units == ureg.Magnetic_flux
    assert cwp.rho0.units == ureg.meters
    assert cwp.Ptheta0.units == ureg.Magnetic_flux
    assert cwp.E.units == ureg.Joule
    assert cwp.EkeV.units == ureg.keV

    # NU counterparts
    assert cwp.miNU.units == ureg.Proton_mass
    assert cwp.qiNU.units == ureg.Proton_charge

    assert cwp.RNU.units == ureg.NUmeters
    assert cwp.aNU.units == ureg.NUmeters
    assert cwp.B0NU.units == ureg.NUTesla
    assert cwp.psi_wallNU.units == ureg.NUMagnetic_flux
    assert cwp.psip_wallNU.units == ureg.NUMagnetic_flux

    assert cwp.psi0NU.units == ureg.NUMagnetic_flux
    assert cwp.Pzeta0NU.units == ureg.NUMagnetic_flux
    assert cwp.t_evalNU.units == ureg.NUseconds
    assert cwp.muNU.units == ureg.NUMagnetic_moment
    assert cwp.muBNU.units == ureg.NUJoule
    assert cwp.psip0NU.units == ureg.NUMagnetic_flux
    assert cwp.rho0NU.units == ureg.NUmeters
    assert cwp.Ptheta0NU.units == ureg.NUMagnetic_flux
    assert cwp.ENU.units == ureg.NUJoule

    # Orbit arrays
    assert cwp.theta.units == ureg.radians
    assert cwp.zeta.units == ureg.radians
    assert cwp.psi.units == ureg.Magnetic_flux
    assert cwp.psip.units == ureg.Magnetic_flux
    assert cwp.rho.units == ureg.meters
    assert cwp.Ptheta.units == ureg.Magnetic_flux
    assert cwp.Pzeta.units == ureg.Magnetic_flux

    assert cwp.psiNU.units == ureg.NUMagnetic_flux
    assert cwp.psipNU.units == ureg.NUMagnetic_flux
    assert cwp.rhoNU.units == ureg.NUmeters
    assert cwp.PthetaNU.units == ureg.NUMagnetic_flux
    assert cwp.PzetaNU.units == ureg.NUMagnetic_flux
