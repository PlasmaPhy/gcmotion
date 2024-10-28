import gcmotion as gcm
import numpy as np

# To be passed to objects as parameters for ease
R = 12
a = 2
qfactor = gcm.qfactor.Hypergeometric(R, a)
N = 10  # Number of particles

# fmt: off

params = {

    "R"           :   R,
    "a"           :   a,
    "qfactor"           :   gcm.qfactor.Hypergeometric(R, a),
    "Bfield"      :   gcm.bfield.LAR(i = 0, g = 1, B0 = 5),
    "Efield"      :   gcm.efield.Radial(R, a, qfactor, Ea = 75000, minimum = 0.98, waist_width=20),
    "species"     :   "p",
    "mu"          :   1e-5,
    "theta0"      :   np.concat((np.pi*np.ones(3), np.zeros(5))),
    "psi0"        :   np.concat((np.linspace(0.9, 1, 3), np.linspace(0.52 ,0.87, 5))), #times psi_wall
    "z0"          :   0,
    "Pz0"         :   -0.025,
    "t_eval"      :   np.linspace(0, 100000, 100000) # t0, tf, steps

    
}

# fmt: on
