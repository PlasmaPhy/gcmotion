from collections import namedtuple


# System of equations to be solved
def system(theta: float, psi: float, parameters: namedtuple, profile: namedtuple, method: str):

    # Tokamak profile
    qfactor = profile.qfactor
    bfield = profile.bfield
    efield = profile.efield

    # Define quantites for the solver for clarity
    solverqNU = qfactor.solverqNU
    psipNU = qfactor.psipNU
    solverbNU = bfield.solverbNU
    solverPhiderNU = efield.solverPhiderNU

    # Parameters ()
    mu = parameters.mu
    Pzeta0 = parameters.Pzeta0
    Pzeta = Pzeta0

    psi = max(psi, 1e-3)

    # Object methods calls
    q = solverqNU(psi)
    psi_p = psipNU(psi)
    b, b_der, currents, currents_der = solverbNU(psi, theta)
    phi_der_psi, phi_der_theta = solverPhiderNU(psi, theta)

    # Unpack
    b_der_psi, b_der_theta = b_der
    i, g = currents
    i_der, g_der = currents_der
    # Multiply current derivatives by q to get their derivatives with
    # respect to psi instead of psip
    i_der, g_der = q * i_der, q * g_der

    # Intermediate values
    rho = (Pzeta + psi_p) / g
    par = mu + rho**2 * b
    bracket1 = -par * b_der_psi + phi_der_psi
    bracket2 = par * b_der_theta + phi_der_theta
    D = g * q + i + rho * (g * i_der - i * g_der)

    # Canonical Equations
    theta_dot = (1 - rho * g_der) / D * rho * b**2 + q * g / D * bracket1
    psi_dot = -q * g / D * bracket2

    if method == "differential evolution":
        return theta_dot**2 + psi_dot**2

    elif method == "fsolve":
        return [theta_dot, psi_dot]
