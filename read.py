"""
General functions used for building the WAVES mock catalogue.
"""

import os
from collections import defaultdict

import h5py
import numpy as np


def read_lightcone(
    model_dir: str, sub_dir: str, fields: dict, sub_volumes: np.ndarray, file_name: str
) -> list[np.ndarray]:
    """Read the mock file for the given model/subvolume as efficiently as possible."""
    data = defaultdict(list)
    for sub_volume in sub_volumes:
        full_name = os.path.join(
            model_dir, sub_dir, f"{file_name}_{sub_volume:02d}.hdf5"
        )
        print(f"Reading galaxies data from {full_name}")
        with h5py.File(full_name, "r") as f:
            for group_name, data_names in fields.items():
                group = f[group_name]
                for data_name in data_names:
                    full_data_path = f"{group_name}/{data_name}"
                    data[full_data_path].append(group[data_name][()])

    for key in data:
        data[key] = np.concatenate(data[key], axis=0)

    return list(data.values())


def read_photometry_data_hdf5(
    model_dir: str, sub_dir: str, fields: dict, sub_volumes: np.ndarray, sed_file: str
) -> tuple[np.ndarray, list]:
    """Read the Sting-SED*.hdf5 file for the given model/subvolume"""

    data = defaultdict(list)
    ids = []

    for subv in sub_volumes:
        full_name = os.path.join(model_dir, sub_dir, f"{sed_file}_{subv:02d}.hdf5")
        print(f"Reading galaxies data from {full_name}")

        with h5py.File(full_name, "r") as f:
            ids.append(f["id_galaxy_sky"][()])

            for group_name, data_names in fields.items():
                group = f[group_name]
                for data_name in data_names:
                    full_data_path = f"{group_name}/{data_name}"
                    data[full_data_path].append(group[data_name][()])

    ids = np.concatenate(ids, axis=0)
    for key in data:
        data[key] = np.concatenate(data[key], axis=1)

    return ids, list(data.values())


def read_filter_names(model_dir: str, sub_dir: str, sed_file: str) -> list[str]:
    """
    Returns a list of the filter names for the given sed file.
    """
    full_name = os.path.join(model_dir, sub_dir, f"{sed_file}_00.hdf5")
    with h5py.File(full_name, "r") as f:
        filter_names = f["filters"][:]
        filter_names = [
            filter_name.decode() for filter_name in filter_names
        ]  # removing byte encoding
    return filter_names
