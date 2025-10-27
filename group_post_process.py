"""
Module for post-processing group ids
"""

from typing import Union
from astropy.cosmology import Planck18
from astropy.units import Quantity
import astropy.units as u
import numpy as np


def calc_rvir_from_mvir(
    mvir: Union[Quantity, np.ndarray[Quantity]],
    redshift: Union[float, np.ndarray[float]],
) -> Union[float, np.ndarray[float]]:
    """
    Calculates the virial radius from a given virial mass.

    The mass must be given in whatever mass units and then will be returned in Mpc
    without the astropy units.

    Using the Plank18 cosmology like shark does.
    """
    numerator = 3 * mvir
    denominator = 4 * np.pi * 200 * Planck18.critical_density(redshift)
    rvir_units = (numerator/denominator) ** (1 / 3)
    return rvir_units.to(u.Mpc).value
