import gcmotion as gcm
import gcmotion.plot as gplt
import numpy as np

# Quantity Constructor
species = "p"
smart_init = gcm.SmartPositiveInit(species)
Q = smart_init.get_QuantityConstructor()

# Intermediate Quantities
R = smart_init.R
a = smart_init.a
B0 = smart_init.B0
i = Q(0, "NUPlasma_current")
g = Q(1, "NUPlasma_current")
Ea = Q(735000, "Volts/meter")

# Construct a Tokamak
tokamak = gcm.Tokamak(
    R=R,
    a=a,
    qfactor=gcm.qfactor.SmartPositive(),
    bfield=gcm.bfield.SmartPositive(),
    efield=gcm.efield.Nofield(),
)

# Setup Initial Conditions
init = gcm.InitialConditions(
    species="p",
    muB=Q(1e-3, "NUMagnetic_moment"),
    Pzeta0=Q(-0.045, "NUMagnetic_flux"),
    theta0=0,
    zeta0=0,
    psi0=Q(0.3, "psi_wall"),
    t_eval=Q(np.linspace(0, 3e-4, 10000), "seconds"),
)

# Create the particle and calculate its obrit
particle = gcm.Particle(tokamak=tokamak, init=init)
particle.run()
print(particle)
