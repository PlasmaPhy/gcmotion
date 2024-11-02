from pint import UnitRegistry

ureg = UnitRegistry(case_sensitive=False)
ureg.setup_matplotlib()

# Capitalize the "T" in "Tesla" because
ureg.define("Tesla = tesla = T")


# fmt: off
def setup_pint(R, B0):

    mp = 1.672621923e-27
    qp = 1.602176634e-19

    ureg.define(f"Proton_mass = {mp} kilogram")  # Proton mass [kg]
    ureg.define(f"Proton_charge = {qp} coulomb")  # Proton charge [C]

    w0 = qp * B0 / mp
    E0 = mp * w0**2 * R**2

    ureg.define(f"NUmeter = {R} meter") # Tokamak major radius 
    ureg.define(f"NUhz = {w0} hz") # Cyclotron frequency 
    ureg.define(f"NUsecond = {1/w0} second") # Inverse cyclotron frequency 
    ureg.define(f"NUTesla = {1} Proton_mass * NUhz / Proton_charge")  # Magnetic field strength 
    ureg.define(f"NUJoule = {E0} Joule") # Energy 

    # Needed SI quantites
    ureg.define("Magnetic_flux = Tesla * m^2")
    ureg.define("Magnetic_moment = Ampere * m^2")
    ureg.define("Plasma_current = Tesla * m")
    ureg.define("Canonical_momentum = Magnetic_flux")

    # Needed NU quantities
    ureg.define("NUMagnetic_flux = NUTesla * NUmeter^2")
    ureg.define("NUPlasma_current = Proton_mass * NUmeter^2 / Proton_charge")
    ureg.define("NUMagnetic_moment = Proton_charge / NUsecond * NUmeter^2")
    ureg.define("NUCanonical_momentum = NUMagnetic_flux")
    ureg.define("NUVolts = NUJoule / Proton_charge")

    return ureg, ureg.Quantity
