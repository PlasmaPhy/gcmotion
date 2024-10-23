import gcmotion as gcm
import numpy as np

# To be passed to objects as parameters for ease
R = 12
a = 2
q = gcm.qfactor.Hypergeometric(R, a)
N = 10  # Number of particles

# fmt: off

params = {

    "R"           :   R,
    "a"           :   a,
    "q"           :   gcm.qfactor.Hypergeometric(R, a),
    "Bfield"      :   gcm.bfield.LAR(i = 0, g = 1, B0 = 5),
    "Efield"      :   gcm.efield.Radial(R, a, q, Ea = 75000, minimum = 0.95, waist_width=50),
    "species"     :   "p",
    "mu"          :   1e-5,
    "theta0"      :   np.concat((np.zeros(6), np.pi*np.ones(3))),
    "psi0"        :   np.concat((np.linspace(0.5, 1.3, 6), np.linspace(0.55, 0.95, 3))), #times psi_wall
    "z0"          :   0,
    "Pz0"         :   -0.025,
    "t_eval"      :   np.linspace(0, 100000, 10000) # t0, tf, steps

    
}

# fmt: on
