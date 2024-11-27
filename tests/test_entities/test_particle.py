import gcmotion as gcm
import numpy as np
from pint import get_application_registry


# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
Q = gcm.QuantityConstructor(R=Rnum, a=anum, B0=B0num, species=species)
ureg = get_application_registry()

# ========================== BASE OBJECTS SETUP ===============================

R = Q(Rnum, "meters")
a = Q(anum, "meters")
B0 = Q(B0num, "Tesla")
i = Q(0, "NUPlasma_current")
g = Q(1, "NUPlasma_current")
Ea = Q(73500, "Volts/meter")

tokamak = gcm.Tokamak(
    R=R,
    a=a,
    qfactor=gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2),
    bfield=gcm.bfield.LAR(B0=B0, i=i, g=g),
    efield=gcm.efield.Radial(a, Ea, B0, peak=0.98, rw=1 / 50),
)

params = gcm.PhysicalParameters(
    species=species,
    mu=Q(1e-5, "NUMagnetic_moment"),
    Pzeta=Q(-0.0272, "NUMagnetic_flux"),
)

init = gcm.InitialConditions(
    muB=Q(0.5, "keV"),
    theta0=0,
    zeta0=0,
    psi0=Q(0.98, "psi_wall"),
    t_eval=Q(np.linspace(0, 1e-3, 1000), "seconds"),
)
# ====================== "COMPOSITE" OBJECT SETUP ============================

profile = gcm.Profile(tokamak, params)

particle = gcm.Particle(profile, init)
particle.run()
# ================================== TESTS ==================================


def test_attrs_units():
    # in the order of definition
    assert particle.mi.units == ureg.kilogram
    assert particle.qi.units == ureg.Coulomb

    assert particle.R.units == ureg.meters
    assert particle.a.units == ureg.meters
    assert particle.B0.units == ureg.Tesla
    assert particle.psi_wall.units == ureg.Magnetic_flux
    assert particle.psip_wall.units == ureg.Magnetic_flux

    assert isinstance(particle.theta0, (int, float))
    assert isinstance(particle.zeta0, (int, float))
    assert particle.psi0.units == ureg.Magnetic_flux
    assert particle.Pzeta.units == ureg.Magnetic_flux
    assert particle.t_eval.units == ureg.seconds
    assert particle.mu.units == ureg.keV / ureg.Tesla
    assert particle.muB.units == ureg.keV
    assert particle.psip0.units == ureg.Magnetic_flux
    assert particle.rho0.units == ureg.meters
    assert particle.Ptheta0.units == ureg.Magnetic_flux
    assert particle.E.units == ureg.Joule
    assert particle.EkeV.units == ureg.keV

    # NU counterparts
    assert particle.miNU.units == ureg.Proton_mass
    assert particle.qiNU.units == ureg.Proton_charge

    assert particle.RNU.units == ureg.NUmeters
    assert particle.aNU.units == ureg.NUmeters
    assert particle.B0NU.units == ureg.NUTesla
    assert particle.psi_wallNU.units == ureg.NUMagnetic_flux
    assert particle.psip_wallNU.units == ureg.NUMagnetic_flux

    assert particle.psi0NU.units == ureg.NUMagnetic_flux
    assert particle.PzetaNU.units == ureg.NUMagnetic_flux
    assert particle.t_evalNU.units == ureg.NUseconds
    assert particle.muNU.units == ureg.NUMagnetic_moment
    assert particle.muBNU.units == ureg.NUJoule
    assert particle.psip0NU.units == ureg.NUMagnetic_flux
    assert particle.rho0NU.units == ureg.NUmeters
    assert particle.Ptheta0NU.units == ureg.NUMagnetic_flux
    assert particle.ENU.units == ureg.NUJoule

    # Orbit arrays
    assert particle.theta.units == ureg.radians
    assert particle.zeta.units == ureg.radians
    assert particle.psi.units == ureg.Magnetic_flux
    assert particle.psip.units == ureg.Magnetic_flux
    assert particle.rho.units == ureg.meters
    assert particle.Ptheta.units == ureg.Magnetic_flux
    assert particle.Pzeta.units == ureg.Magnetic_flux

    assert particle.psiNU.units == ureg.NUMagnetic_flux
    assert particle.psipNU.units == ureg.NUMagnetic_flux
    assert particle.rhoNU.units == ureg.NUmeters
    assert particle.PthetaNU.units == ureg.NUMagnetic_flux
    assert particle.PzetaNU.units == ureg.NUMagnetic_flux


def test_quantitiesfunc():
    """Make sure it doesn't print any errors"""
