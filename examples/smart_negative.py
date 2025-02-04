import gcmotion as gcm
import gcmotion.plot as gplt
import numpy as np

# Quantity Constructor
species = "p"
smart_init = gcm.SmartNegativeInit(species)
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
    # qfactor=gcm.qfactor.Parabolic(a=a, B0=B0, q0=1, q_wall=5),
    qfactor=gcm.qfactor.SmartNegative(),
    # bfield=gcm.bfield.LAR(B0, i, g),
    bfield=gcm.bfield.SmartNegative(),
    efield=gcm.efield.Nofield(),
)

# Setup Initial Conditions
init = gcm.InitialConditions(
    species="p",
    muB=Q(1e-5, "NUMagnetic_moment"),
    Pzeta0=Q(-0.05, "NUMagnetic_flux"),
    theta0=0,
    zeta0=0,
    # psi0=Q(0.35, "psi_wall"),
    psi0=Q(0.005123, "Magnetic_flux"),
    # psi0=Q(0.004889, "Magnetic_flux"),
    t_eval=Q(np.linspace(0, 1e-3, 100000), "seconds"),
)

# Create the particle and calculate its obrit
particle = gcm.Particle(tokamak=tokamak, init=init)
event = gcm.events.when_theta(0, 10)
particle.run(events=[event])
print(particle)

# gplt.particle_poloidal_drift(particle, projection="polar")
gplt.particle_poloidal_drift(particle, psilim=[0.2, 0.25])
