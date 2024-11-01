import gcmotion as gcm
import numpy as np

# --------------------Particle Setup--------------------
R, a = 1.65, 0.5  # Major/Minor Radius in [m]
qfactor = gcm.qfactor.Hypergeometric(R, a, q0=1.1, q_wall=3.5, n=2)
Bfield = gcm.bfield.LAR(i=0, g=1, B0=1)
Efield = gcm.efield.Radial(R, a, qfactor, Ea=75000, minimum=0.9, r_w=1 / 50)

species = "d"
mu = 1e-5
theta0 = 0
psi0 = 0.533  # times psi_wall
zeta0 = np.pi
Pzeta0 = -0.022
t_eval = np.linspace(0, 1000000, 100000)  # t0, tf, steps

tokamak = {
    "R": R,
    "a": a,
    "qfactor": qfactor,
    "Bfield": Bfield,
    "Efield": Efield,
}
parameters = {
    "species": species,
    "mu": mu,
    "theta0": theta0,
    "psi0": psi0,
    "zeta0": zeta0,
    "Pzeta0": Pzeta0,
    "t_eval": t_eval,
}

particle1 = gcm.Particle(tokamak, parameters)
cwp = particle1

cwp.run()

# -------------------------Plots-------------------------

gcm.tokamak_profile(cwp, zoom=[0, 1.1])

gcm.time_evolution(cwp, percentage=20, units="s")

gcm.drift(cwp, angle="theta", lim=[-np.pi, np.pi], plot_initial=True)

gcm.drifts(cwp, theta_lim=[-np.pi, np.pi], plot_initial=True)

gcm.energy_contour(
    cwp,
    theta_lim=[-np.pi, np.pi],
    psi_lim="auto",
    plot_drift=True,
    contour_Phi=True,
    units="keV",
    levels=20,
)

gcm.parabolas(cwp)

gcm.poloidal_cut(cwp, wall_shade=True)

gcm.torus2d(cwp, percentage=20)

gcm.torus3d(cwp, bold="bold", truescale=True, percentage=50)
