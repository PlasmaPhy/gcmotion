import gcmotion as gcm
import numpy as np

# To be passed to objects as parameters for ease
R = 1.65
a = 0.5
q = gcm.qfactor.Hypergeometric(R, a, q0=1.1, q_wall=3.5, n=2)
N = 20  # Number of particles

# fmt: off

params = {

    "R"           :   R,
    "a"           :   a,
    "q"           :   q,
    "Bfield"      :   gcm.bfield.LAR(i = 0, g = 1, B0 = 1),
    "Efield"      :   gcm.efield.Radial(R, a, q, Ea = 73500, minimum = 0.98, r_w=1/50),
    "species"     :   "p",
    "mu"          :   0.000019173,
    "theta0"      :   np.repeat(0*np.pi, N),
    "psi0"        :   np.linspace(0.1, 1, N), #times psi_wall
    "z0"          :   0,
    "Pz0"         :   -0.015,
    "t_eval"      :   np.linspace(0, 100000, 1000000) # t0, tf, steps

    
}

# fmt: on
