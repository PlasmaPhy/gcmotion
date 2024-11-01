from pint import UnitRegistry

ureg = UnitRegistry(case_sensitive=False)

# Capitalize the "T" in "Tesla" because
ureg.define("Tesla = tesla = T")


def pint_setup(R, B0):

    mp = 1.672621923e-27  # Proton mass [kg]
    qp = 1.602176634e-19  # Proton charge [C]

    w0 = qp * B0 / mp
    E0 = mp * w0**2 * R**2

    ureg.define(f"NUmeter = {R} meters")
    ureg.define(f"NUsecond = {1/w0} seconds")
    ureg.define(f"NUTesla = {mp*w0/qp} Tesla")  # = 1 Tesla
    ureg.define(f"NUEnergy = {E0} Joule")

    return ureg
