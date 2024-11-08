import gcmotion as gcm
import numpy as np

Rnum = 1.65
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, B0=B0num, species=species)

R = Q(Rnum, "meters")
a = Q(0.5, "meters")
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
    # "mu": Q(1e-5, "NUmagnetic_moment"),
    "muB": Q(5, "keV"),
    "theta0": 0,
    "zeta0": 0,
    "psi0": 0.8,  # times psi_wall
    "Pzeta0": Q(-0.0272, "NUMagnetic_flux"),
    "t_eval": Q(np.linspace(0, 1e-3, 10000), "seconds"),
}

cwp = gcm.Particle(tokamak, parameters)
cwp.run()

gcm.time_evolution(cwp)

gcm.drift(cwp)

gcm.drifts(cwp)

gcm.energy_contour(cwp, units="SI")
