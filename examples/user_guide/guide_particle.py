import gcmotion as gcm
import numpy as np

# ========================== QUANTITY CONSTRUCTOR ==========================

Rnum = 1.65
anum = 0.5
B0num = 1
species = "p"
ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

# ============================= PARTICLE SETUP =============================

R = Q(Rnum, "meters")
a = Q(anum, "meters")
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

# ================================= PLOTS =================================

gcm.qfactor_profile(tokamak, Q, units="SI")

gcm.magnetic_profile(tokamak, Q)

gcm.electric_profile(tokamak, Q, zoom=[0.8, 1.1])

gcm.time_evolution(cwp)

gcm.drift(cwp)

gcm.drifts(cwp)

gcm.energy_contour(cwp, wall_shade=True)

gcm.poloidal_cut(cwp)

gcm.torus2d(cwp)

gcm.torus3d(cwp, bold="bold")
