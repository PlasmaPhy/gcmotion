import gcmotion as gcm
import numpy as np

# To be passed to objects as parameters for ease
R = 1.65
a = 0.5
qfactor = gcm.qfactor.Hypergeometric(R, a, q0=1.1, q_wall=3.5, n=2)
N = 15  # Number of particles

# fmt: off

params = {

    "R"           :   R,
    "a"           :   a,
    "qfactor"     :   qfactor,
    "Bfield"      :   gcm.bfield.LAR(i = 0, g = 1, B0 = 1),
    "Efield"      :   gcm.efield.Radial(R, a, qfactor, Ea = 73500, peak = 0.98, r_w=1/50),
    "species"     :   "he3",
    "mu"          :   0.0000019173,
    "theta0"      :   np.concat((np.repeat(0, 6), np.repeat(np.pi, 2))),
    "psi0"        :   np.concat((np.linspace(0.81, 1.1, 6), np.linspace(0.9, 1, 2))), #times psi_wall
    "zeta0"       :   0,
    "Pzeta0"      :   -0.0272,
    "t_eval"      :   np.linspace(0, 100000, 1000000) # t0, tf, steps

    
}

# fmt: on
