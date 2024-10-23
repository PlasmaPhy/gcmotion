import gcmotion as gcm
import numpy as np

R, a = 12, 2  # Major/Minor Radius in [m]
q = gcm.qfactor.Hypergeometric(R, a)
Bfield = gcm.bfield.LAR(i=0, g=1, B0=5)
Efield = gcm.efield.Radial(R, a, q, Ea=75000, minimum=0.9, waist_width=50)
