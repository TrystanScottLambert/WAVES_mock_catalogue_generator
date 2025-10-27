"""
Module for post-processing group ids
"""

from typing import Union
from dataclasses import dataclass

from astropy.cosmology import Planck18
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
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
    rvir_units = (numerator / denominator) ** (1 / 3)
    return rvir_units.to(u.Mpc).value


@dataclass
class Group:
    """
    Data class representing a group.
    """
    mass: float
    ra: float
    dec: float
    zcos: float

    def __post_init__(self):
        self.rvir = calc_rvir_from_mvir(self.mass * u.solMass, self.zcos)
        distance = Planck18.comoving_distance(self.zcos)
        c = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg, distance=distance)
        self.center = np.array(
            [c.cartesian.x.value, c.cartesian.y.value, c.cartesian.z.value]
        )

    def overlap(self, other_group: "Group") -> bool:
        """
        Determins group is overlapping with other_group
        """
        difference = self.center - other_group.center
        distance_between_centers = np.sqrt((difference[0])**2 + (difference[1])**2 + (difference[2])**2)
        total_radii = self.rvir + other_group.rvir
        overlapped = distance_between_centers < total_radii
        return overlapped