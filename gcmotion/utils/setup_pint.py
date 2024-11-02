from pint import UnitRegistry

ureg = UnitRegistry(case_sensitive=False)
ureg.setup_matplotlib()

# Capitalize the "T" in "Tesla" because
ureg.define("Tesla = tesla = T")


# fmt: off
def setup_pint(R, B0):

    mp = 1.672621923e-27
    qp = 1.602176634e-19

    ureg.define(f"NUmass = {mp} kilogram = mp")  # Proton mass [kg]
    ureg.define(f"NUcharge = {qp} Coulombs = qp")  # Proton charge [C]

    w0 = qp * B0 / mp
    E0 = mp * w0**2 * R**2

    ureg.define(f"NUmeter = {R} meters = R") # Tokamak major radius [m]
    ureg.define(f"NUsecond = {1/w0} seconds") # Inverse cyclotron frequency [s]
    ureg.define(f"NUTesla = {mp*w0/qp} Tesla = B0")  # Magnetic field strength [T]
    ureg.define(f"NUEnergy = {E0} Joule") # Energy [J]

    # Needed SI quantites
    ureg.define("Magnetic_flux = Tesla * m^2")
    ureg.define("Magnetic_moment = Ampere * m^2")
    ureg.define("Plasma_current = Tesla * m")
    ureg.define("Canonical_momentum = Magnetic_flux")

    # Needed NU quantities
    ureg.define("NUMagnetic_flux = NUTesla * NUmeter^2")
    ureg.define("NUCurrent = mp * NUmeter**2 / qp")
    ureg.define("NUMagnetic_moment = qp / NUsecond * NUmeter^2")
    ureg.define("NUVolts = NUEnergy / qp")

    return ureg, ureg.Quantity
