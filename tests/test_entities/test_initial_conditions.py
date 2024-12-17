import gcmotion as gcm
import numpy as np

# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
Q = gcm.QuantityConstructor(R=Rnum, a=anum, B0=B0num, species=species)

# ========================== BASE OBJECTS SETUP ===============================

R = Q(Rnum, "meters")
a = Q(anum, "meters")
B0 = Q(B0num, "Tesla")
i = Q(0, "NUPlasma_current")
g = Q(1, "NUPlasma_current")
Ea = Q(73500, "Volts/meter")

R = Q(Rnum, "meters")
a = Q(anum, "meters")
B0 = Q(B0num, "Tesla")
i = Q(0, "NUPlasma_current")
g = Q(1, "NUPlasma_current")

tokamak = gcm.Tokamak(
    R=R,
    a=a,
    qfactor=gcm.qfactor.Unity(),
    bfield=gcm.bfield.LAR(B0=B0, i=i, g=g),
    efield=gcm.efield.Nofield(),
)


def test_init_str_repr_functionality():
    """Tests that __repr__() and __str__() can return with no errors"""
    _ = init.__repr__()
    _ = init.__str__()


init = gcm.InitialConditions(
    species=species,
    muB=Q(0.5, "keV"),
    theta0=0,
    zeta0=0,
    psi0=Q(0.6, "psi_wall"),
    Pzeta0=Q(-0.015, "NUMagnetic_flux"),
    t_eval=Q(np.linspace(0, 1e-3, 1000), "seconds"),
)

init_full = init._calculate_full_set(tokamak)


def test_init_full_str_repr_functionality():
    """Tests that __repr__() and __str__() can return with no errors"""
    _ = init_full.__repr__()
    _ = init_full.__str__()
