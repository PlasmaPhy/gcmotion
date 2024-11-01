import gcmotion as gcm
import numpy as np

# To be passed to objects as parameters for ease
R = 1.65
a = 0.5
qfactor = gcm.qfactor.Hypergeometric(R, a, q0=1.1, q_wall=3.5, n=2)
N = 16  # Number of particles

# fmt: off

params = {

    "R"           :   R,
    "a"           :   a,
    "qfactor"     :   qfactor,
    "Bfield"      :   gcm.bfield.LAR(i = 0, g = 1, B0 = 1),
    "Efield"      :   gcm.efield.Radial(R, a, qfactor, Ea = 73500, minimum = 0.9, r_w=1/50),    
    "species"     :   "He3",
    "mu"          :   1e-5,
    "theta0"      :   np.repeat(0, N),
    "psi0"        :   np.linspace(0.5, 1, N), #times psi_wall
    "zeta0"       :   0,
    "Pzeta0"      :   -0.025,
    "t_eval"      :   np.linspace(0, 100000, 100000) # t0, tf, steps
    
}

# fmt: on
