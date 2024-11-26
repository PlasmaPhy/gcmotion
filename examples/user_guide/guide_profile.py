import gcmotion as gcm
import gcmotion.plotters as gplt
import numpy as np

# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

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
    muB=Q(5, "keV"),
    theta0=0,
    zeta0=0,
    psi0=Q(0.78, "psi_wall"),
    t_eval=Q(np.linspace(0, 1e-3, 10000), "seconds"),
)

# ====================== "COMPOSITE" OBJECT SETUP ============================

profile = gcm.Profile(tokamak, params)

gplt.profile_contour(profile)
