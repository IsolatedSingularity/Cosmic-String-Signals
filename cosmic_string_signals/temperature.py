"""Brightness temperature models for cosmic string wake signals."""

import numpy as np


def photon_cmb_temperature(z: float) -> float:
    """CMB photon temperature at redshift z.

    Parameters
    ----------
    z : float
        Cosmological redshift.

    Returns
    -------
    float
        Temperature in Kelvin.
    """
    return 2.72548 * (1 + z)


def kinetic_temperature(z: float, g_mu: float, gamma_v: float) -> float:
    """Kinetic temperature inside the wake at redshift z.

    Parameters
    ----------
    z : float
        Cosmological redshift.
    g_mu : float
        Dimensionless string tension.
    gamma_v : float
        Product of Lorentz factor and string speed (gamma * v_s).

    Returns
    -------
    float
        Kinetic temperature in Kelvin.
    """
    return 20 * (g_mu * 1e6) ** 2 * gamma_v ** 2 * (1100 + 1) / (z + 1)


def gas_temperature(z: float) -> float:
    """Gas temperature scaling with redshift.

    Parameters
    ----------
    z : float
        Cosmological redshift.

    Returns
    -------
    float
        Gas temperature in Kelvin.
    """
    return 0.02 * (z + 1) ** 2


def brightness_temperature(z: float, g_mu: float, gamma_v: float) -> float:
    """21cm brightness temperature fluctuation at redshift z.

    Computes the differential brightness temperature induced by a cosmic
    string wake, following the model in Brandenberger et al.

    Parameters
    ----------
    z : float
        Cosmological redshift.
    g_mu : float
        Dimensionless string tension.
    gamma_v : float
        Product of Lorentz factor and string speed.

    Returns
    -------
    float
        Brightness temperature fluctuation in Kelvin.
    """
    deexcitation_cross_section = 0.16
    t_kin = kinetic_temperature(z, g_mu, gamma_v)
    t_gas = gas_temperature(z)

    u = t_kin if t_kin > 3 * t_gas else 3 * t_gas
    c = 4 if t_kin > 3 * t_gas else 1 + t_kin / t_gas

    t_gamma = photon_cmb_temperature(z)
    result = (
        0.017
        * deexcitation_cross_section
        / (1 + deexcitation_cross_section)
        * (1 - t_gamma / u)
        * np.sqrt(1 + z)
        * c
    )
    return result
