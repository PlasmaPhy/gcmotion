from collections import namedtuple


# Create "Profile" template
Profile = namedtuple(
    "Profile", ["R", "a", "B0", "qfactor", "bfield", "efield"]
)

InitialCondtions = namedtuple(
    "InitialConditions",
    ["species", "mu/muB", "theta0", "zeta0", "Psi0", "Pzeta0", "t_eval"],
)
