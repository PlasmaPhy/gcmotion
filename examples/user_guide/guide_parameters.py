import gcmotion as gcm
import numpy as np

N = 10
Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

# ==========================================================

R = Q(Rnum, "meters")
a = Q(anum, "meters")
B0 = Q(B0num, "Tesla")
i = Q(0, "NUPlasma_current")
g = Q(1, "NUPlasma_current")
Ea = Q(73500, "Volts/meter")

# fmt: off

params = {

    "R"           :   R,
    "a"           :   a,
    "B0"          :   B0,
    "qfactor"     :   gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2),
    "bfield"      :   gcm.bfield.LAR(B0=B0, i=i, g=g),
    "efield"      :   gcm.efield.Radial(a, Ea, B0, peak=0.98, rw=1/50),    

    "species"     :   species,
    "mu/muB"      :   Q(5, "keV"),
    "theta0"      :   np.repeat(0, 15),
    #"theta0"      :   np.concat((np.repeat(0, 5), np.repeat(np.pi, 3))),
    "zeta0"       :   0,
    "psi0"        :   Q(np.linspace(0.3, 1, 15), "psi_wall"),
    #"psi0"        :   Q(np.concat((np.linspace(0.6, 1.2, 5), np.linspace(0.95, 1.05,3))), "psi_wall"),
    "Pzeta0"      :   Q(-0.0272, "NUMagnetic_flux"),
    "t_eval"      :   Q(np.linspace(0, 1e-3, 100000), "seconds") # t0, tf, steps
    
}

# fmt: on
