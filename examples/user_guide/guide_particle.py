import gcmotion as gcm
import numpy as np

R, a = 12, 2  # Major/Minor Radius in [m]
qfactor = gcm.qfactor.Hypergeometric(R, a, q0=1.1, q_wall=3.5, n=2)
Bfield = gcm.bfield.LAR(i=0, g=1, B0=5)
Efield = gcm.efield.Radial(R, a, qfactor, Ea=75000, minimum=0.9, r_w=1 / 50)

species = "p"
mu = 1e-5
theta0 = 0
psi0 = 0.5  # times psi_wall
zeta0 = np.pi
Pzeta0 = -0.027
t_eval = np.linspace(0, 10000, 10000)  # t0, tf, steps

tokamak = {
    "R": R,
    "a": a,
    "qfactor": qfactor,
    "Bfield": Bfield,
    "Efield": Efield,
}
init_cond = {"theta0": theta0, "psi0": psi0, "zeta0": zeta0, "Pzeta0": Pzeta0}

particle1 = gcm.Particle(tokamak, t_eval, init_cond, mu, species)
cwp = particle1

cwp.run()

gcm.time_evolution(cwp, percentage=100, units="s")
gcm.tokamak_profile(cwp, zoom=[0, 1.1])
gcm.drift(cwp, angle="theta", lim=[-np.pi, np.pi])
gcm.drifts(cwp)
gcm.poloidal_cut(cwp)