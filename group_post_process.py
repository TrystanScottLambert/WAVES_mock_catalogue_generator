"""
Module for post-processing group ids - fully vectorized version
"""

from typing import Union

from astropy.cosmology import Planck18
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import networkx as nx
import polars as pl
from scipy.spatial import cKDTree


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


def skycoord_to_cartesian_vectorized(ra, dec, zcos):
    """
    Vectorized conversion of sky coordinates to cartesian.
    """
    distance = Planck18.comoving_distance(zcos)
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=distance)
    centers = np.column_stack(
        [c.cartesian.x.value, c.cartesian.y.value, c.cartesian.z.value]
    )
    return centers


def join_groups(
    centers: np.ndarray, rvirs: np.ndarray, ids: np.ndarray
) -> dict[str, int]:
    """
    Fully vectorized version using KDTree for spatial queries.
    No Group objects needed!

    Parameters:
    -----------
    centers : np.ndarray
        Nx3 array of cartesian coordinates
    rvirs : np.ndarray
        N array of virial radii
    ids : np.ndarray
        N array of group IDs (as strings)
    """
    # Build KDTree for fast spatial queries
    tree = cKDTree(centers)

    # Find all pairs within possible overlap distance
    max_search_radius = 2 * np.max(rvirs)

    print(
        f"Finding overlapping groups (max search radius: {max_search_radius:.3f} Mpc)..."
    )

    # Query pairs within max_search_radius
    pairs = tree.query_pairs(r=max_search_radius, output_type="ndarray")

    print(f"Found {len(pairs)} candidate pairs, checking actual overlaps...")

    # Vectorized overlap check
    i_indices = pairs[:, 0]
    j_indices = pairs[:, 1]

    # Calculate distances for all pairs at once
    distances = np.linalg.norm(centers[i_indices] - centers[j_indices], axis=1)
    total_radii = rvirs[i_indices] + rvirs[j_indices]

    # Find actually overlapping pairs
    overlapping = distances < total_radii
    overlapping_pairs = pairs[overlapping]

    print(f"Found {len(overlapping_pairs)} actual overlaps")

    # Convert indices to IDs for graph
    edges = [(ids[i], ids[j]) for i, j in overlapping_pairs]

    # Build graph and find connected components
    group_graph = nx.Graph()
    group_graph.add_nodes_from(ids)
    group_graph.add_edges_from(edges)

    real_groups = list(nx.connected_components(group_graph))

    mapping = {}
    counter = 1
    for group in real_groups:
        if len(group) == 1:
            mapping[group.pop()] = -1
        else:
            for g in group:
                mapping[g] = counter
            counter += 1

    return mapping


def add_fof_ids(galaxies_file: str, groups_file: str):
    """
    Builds a list of joined groups which represent the ids that can be used to tune group-finders.
    Reads in the already generated mock parquet files.
    """
    print("Reading parquet files...")
    # Read in
    df_galaxies = pl.read_parquet(galaxies_file)
    df_groups = pl.read_parquet(groups_file)

    # Extract group properties as numpy arrays
    print("Extracting group properties...")
    ras = df_groups["ra"].to_numpy()
    decs = df_groups["dec"].to_numpy()
    zcos = df_groups["zcos"].to_numpy()
    mvir = df_groups["mvir"].to_numpy()
    ids = df_groups["id_group_sky"].cast(str).to_numpy()

    print(f"Processing {len(ids)} groups...")

    # Vectorized calculations - all at once!
    print("Computing virial radii (vectorized)...")
    rvirs = calc_rvir_from_mvir(mvir * u.solMass, zcos)

    print("Converting to cartesian coordinates (vectorized)...")
    centers = skycoord_to_cartesian_vectorized(ras, decs, zcos)

    # Generate mapping: id_group_sky (str) → id_fof (int)
    # No Group objects created at all!
    new_group_mapping = join_groups(centers, rvirs, ids)

    # Apply mapping to both dataframes
    print("Applying FoF mapping...")

    # Groups file: all IDs should be in mapping
    df_groups = df_groups.with_columns(
        pl.col("id_group_sky").cast(str).replace(new_group_mapping).alias("id_fof")
    )

    df_galaxies = df_galaxies.with_columns(
        pl.when(
            (pl.col("id_fof").count().over("id_fof") == 1) & (pl.col("id_fof") != -1)
        )
        .then(-1)
        .otherwise(pl.col("id_fof"))
        .alias("id_fof")
    )

    print("Writing results...")
    df_galaxies.write_parquet(galaxies_file)
    df_groups.write_parquet(groups_file)
    print("Done!")
