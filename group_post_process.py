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


def assign_ungrouped_galaxies(
    df_galaxies: pl.DataFrame,
    group_centers: np.ndarray,
    group_rvirs: np.ndarray,
    group_ids: np.ndarray,
) -> pl.DataFrame:
    """
    Assign galaxies with id_group_sky == -1 to the nearest group if within that group's radius.

    Parameters:
    -----------
    df_galaxies : pl.DataFrame
        DataFrame containing galaxy data
    group_centers : np.ndarray
        Nx3 array of group cartesian coordinates
    group_rvirs : np.ndarray
        N array of group virial radii
    group_ids : np.ndarray
        N array of group IDs (as strings)
    """
    # Get ungrouped galaxies
    ungrouped_mask = df_galaxies["id_group_sky"] == -1

    if ungrouped_mask.sum() == 0:
        print("No ungrouped galaxies to assign")
        return df_galaxies

    print(f"Found {ungrouped_mask.sum()} ungrouped galaxies, attempting assignment...")

    # Extract ungrouped galaxy positions
    ungrouped_df = df_galaxies.filter(ungrouped_mask)
    ungrouped_ras = ungrouped_df["ra"].to_numpy()
    ungrouped_decs = ungrouped_df["dec"].to_numpy()
    ungrouped_zcos = ungrouped_df["zcos"].to_numpy()

    # Convert to cartesian
    ungrouped_positions = skycoord_to_cartesian_vectorized(
        ungrouped_ras, ungrouped_decs, ungrouped_zcos
    )

    # Build KDTree for groups
    tree = cKDTree(group_centers)

    # For each ungrouped galaxy, find all groups within max possible radius
    max_search_radius = np.max(group_rvirs)

    # Query for nearby groups (returns list of arrays for each point)
    nearby_indices = tree.query_ball_point(ungrouped_positions, r=max_search_radius)

    # Process each ungrouped galaxy
    new_assignments = []
    assigned_count = 0

    for gal_idx, nearby_group_indices in enumerate(nearby_indices):
        if len(nearby_group_indices) == 0:
            # No groups nearby
            new_assignments.append(-1)
            continue

        # Calculate distances to all nearby groups
        gal_pos = ungrouped_positions[gal_idx]
        distances = np.linalg.norm(
            group_centers[nearby_group_indices] - gal_pos, axis=1
        )

        # Check which groups contain this galaxy
        within_radius = distances < group_rvirs[nearby_group_indices]

        if not np.any(within_radius):
            # Not within any group radius
            new_assignments.append(-1)
        else:
            # Assign to nearest group that contains it
            valid_distances = distances[within_radius]
            valid_indices = np.array(nearby_group_indices)[within_radius]
            nearest_idx = valid_indices[np.argmin(valid_distances)]
            new_assignments.append(group_ids[nearest_idx])
            assigned_count += 1

    print(f"Successfully assigned {assigned_count} ungrouped galaxies to groups")

    # Update the dataframe
    # Create a mapping for ungrouped galaxies
    ungrouped_galaxy_indices = (
        df_galaxies.with_row_index().filter(ungrouped_mask)["index"].to_numpy()
    )

    # Update id_group_sky for those that were assigned
    id_group_sky_array = df_galaxies["id_group_sky"].to_numpy().copy()
    for i, new_id in enumerate(new_assignments):
        if new_id != -1:
            id_group_sky_array[ungrouped_galaxy_indices[i]] = int(new_id)

    df_galaxies = df_galaxies.with_columns(
        pl.Series("id_group_sky", id_group_sky_array)
    )

    return df_galaxies


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
            mapping[group.pop()] = counter
            counter += 1
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

    # PRELIMINARY STEP: Assign ungrouped galaxies to nearby groups
    print("\n=== Assigning ungrouped galaxies ===")
    df_galaxies = assign_ungrouped_galaxies(df_galaxies, centers, rvirs, ids)
    print("=== Ungrouped galaxy assignment complete ===\n")

    # Generate mapping: id_group_sky (str) → id_fof (int)
    # No Group objects created at all!
    new_group_mapping = join_groups(centers, rvirs, ids)
    new_group_mapping["-1"] = -1

    # Apply mapping to both dataframes
    print("Applying FoF mapping...")

    # Groups file: all IDs should be in mapping
    df_groups = df_groups.with_columns(
        pl.col("id_group_sky").cast(str).replace(new_group_mapping).alias("id_fof")
    )

    df_galaxies = df_galaxies.with_columns(
        pl.col("id_group_sky").cast(str).replace(new_group_mapping).alias("id_fof")
    )

   

    df_galaxies = df_galaxies.with_columns(pl.col("id_fof").cast(int))
    df_groups = df_groups.with_columns(pl.col("id_fof").cast(int))

    # All galaxies that only show up once need to be assigned isolated.
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
