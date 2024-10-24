import gcmotion as gcm
import numpy as np

# To be passed to objects as parameters for ease
R = 6
a = 2
q = gcm.qfactor.Hypergeometric(R, a)
N = 20  # Number of particles

# fmt: off

params = {

    "R"           :   R,
    "a"           :   a,
    "q"           :   gcm.qfactor.Hypergeometric(R, a),
    "Bfield"      :   gcm.bfield.LAR(i = 0, g = 1, B0 = 5),
    "Efield"      :   None,#gcm.efield.Radial(R, a, q, Ea = 75000, minimum = 0.94, waist_width=20),
    "species"     :   "p",
    "mu"          :   1e-4,
    "theta0"      :   np.repeat(0, N),
    "psi0"        :   np.linspace(0.01, 0.2, N), #times psi_wall
    "z0"          :   0,
    "Pz0"         :   -0.01,
    "t_eval"      :   np.linspace(0, 100000, 100000) # t0, tf, steps

    
}

# fmt: on
