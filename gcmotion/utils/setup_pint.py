"""
Pint setup
----------

Sets up the `pint <https://pint.readthedocs.io/en/stable/>`_ configuration.
"""

from pint import UnitRegistry
from gcmotion.configuration.particle_attributes import particle_attributes

ureg = UnitRegistry(case_sensitive=False)
ureg.setup_matplotlib()

# fmt: off
def setup_pint(R, a, B0, species):
    r"""Creates all needed [NU] units, as well as some extra [SI] unit 
    aliases and stores them in the UnitRegistry.

    Parameters
    ----------
    R : float
        The tokamak's major radius **in [m]**.
    a : float
        The tokamak's minor radius **in [m]**. Only used to create "psi_wall"
        and "NUpsi_wall" as units of magnetic flux, so it is possible to setup 
        initial :math:`\psi_0` conditions with respect to the wall, instead of
        guessing.
    B0 : float
        The magnetic field strength on the magnetic axis **in [T]**.
    species : str
        The particle species.

    Returns
    -------
    2-tuple
        The updated UnitRegistry and the UnitRegistry.Quantity object, which is 
        used to create all other Quantities.
    """

    # Additional SI quantites (= aliases, for display only)
    ureg.define("Magnetic_flux    = Tesla * m^2   = Tm^2")
    ureg.define("Magnetic_moment  = Ampere * m^2  = keV/T")
    ureg.define("Plasma_current   = Tesla * m     = Tm")

    # Base NU units
    mp = 1.672621923e-27
    qp = 1.602176634e-19
 
    ureg.define(f"Proton_mass   = {mp} kilogram")  # Proton mass [kg]
    ureg.define(f"Proton_charge = {qp} coulomb")  # Proton charge [C]
    
    M = particle_attributes[species.lower() + "_M"]
    Z = particle_attributes[species.lower() + "_Z"]

    w0 = (Z / M) * qp / mp * B0 # s^-1
    E0 = mp * w0**2 * R**2 # Joule

    ureg.define(f"NUsecond = {1/w0} second")  # Time [NU]
    ureg.define(f"NUw0     = {w0} hz")  # Cyclotron frequency
    ureg.define(f"NUmeter  = {R} meter")  # Tokamak major radius
    ureg.define(f"NUJoule  = {E0} Joule")  # Energy
    ureg.define(f"NUkeV    = {qp} NUJoule") # Energy
    ureg.define(f"NUTesla  = {M/Z} Proton_mass * NUw0 / Proton_charge ")  # Magnetic field strength

    # Additional NU quantities
    ureg.define("NUvelocity           = NUmeter * NUw0")
    ureg.define("NUMagnetic_flux      = NUTesla * NUmeter^2                   = NUmf")
    ureg.define("NUPlasma_current     = NUTesla * NUmeter                     = NUpc")
    ureg.define("NUMagnetic_moment    = Proton_charge / NUsecond * NUmeter^2  = NUmu")
    ureg.define("NUVolts              = NUJoule / Proton_charge               = NUV")
    ureg.define("NUVolts_per_NUmeter  = NUVolts / NUmeter                     = NUV/NUm")

    # Also define psi_wall as a unit of Magnetic_flux, to assing psi initial
    # values with respect to it
    ureg.define(f"psi_wall = {B0 * a**2 / 2} Magnetic_flux")
    ureg.define(f"NUpsi_wall = {(a / R)**2 / 2} NUMagnetic_flux") # not really need but sure

    # Assign custom values to Q for easier access.
    ureg.Quantity.w0 = w0
    ureg.Quantity.E0 = E0

    return ureg, ureg.Quantity
