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
import networkx as nx


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
    id: str

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
        distance_between_centers = np.sqrt(
            (difference[0]) ** 2 + (difference[1]) ** 2 + (difference[2]) ** 2
        )
        total_radii = self.rvir + other_group.rvir
        overlapped = distance_between_centers < total_radii
        return overlapped


def join_groups(groups: list[Group]) -> dict[str, int]:
    """
    Iteratively joins groups and returns the mapping to their new ids.
    """
    mapping = {}
    edges = []
    all_ids = [group.id for group in groups]
    group_graph = nx.Graph()
    group_graph.add_nodes_from(all_ids)
    for i in range(len(groups) - 1):
        for j in range(i + 1, len(groups)):
            if groups[i].overlap(groups[j]):
                edges.append((groups[i].id, groups[j].id))

    group_graph.add_edges_from(edges)
    real_groups = list(nx.connected_components(group_graph))
    counter = 1
    for group in real_groups:
        if len(group) == 1:
            mapping[group.pop()] = -1
        else:
            for g in group:
                mapping[g] = counter
            counter+=1
    return mapping
