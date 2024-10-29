"""
General functions used for building the WAVES mock catalogue.
"""

from collections import defaultdict

import h5py
import numpy as np

from load import Config


def read_lightcone(config: Config, source_type: str) -> list[np.ndarray]:
    """Read the mock file for the given model/subvolume as efficiently as possible."""

    if source_type == "gal":
        fields = config.gal_props_read
    elif source_type == "group":
        fields = config.group_props_read
    else:
        raise ValueError(f'type must be either "group" or "gal", not {source_type}')

    data = defaultdict(list)
    for sub_volume in config.dirs.sub_volumes:
        full_name = config.print_full_file_name("mock", sub_volume)
        print(f"Reading galaxies data from: {full_name}")
        with h5py.File(full_name, "r") as f:
            for group_name, data_names in fields.items():
                group = f[group_name]
                for data_name in data_names:
                    data[data_name].append(group[data_name][()])

    for key in data:
        data[key] = np.concatenate(data[key], axis=0)

    return data


def read_filter_names(config: Config) -> list[str]:
    """
    Returns a list of the filter names for the given sed file.
    """
    full_name = config.print_full_file_name("sed", 0)
    with h5py.File(full_name, "r") as f:
        filter_names = f["filters"][:]
        filter_names = [filter_name.decode() for filter_name in filter_names]
    return filter_names


def combine_filters_and_data(filter_names: list[str], filter_data: dict) -> dict:
    """
    Takes the data dictionary from the 'read_photometry_data_hdf5' function and combines it with
    the filter names from 'read_filter_names' and creates a new dictionary where each key represents
    each combination of filters (total_ap_dust_FUV_GALEX, or bulge_t_ab_dust_Band8_ALMA).
    """
    # Restructuring the dictionary
    new_data = {}

    for old_key, arrays in filter_data.items():
        components = old_key.split("/")

        for i, filter_name in enumerate(filter_names):
            new_key = f"{components[2]}_{components[1]}_{filter_name}"
            new_data[new_key] = arrays[i]
    return new_data


def read_photometry_data_hdf5(config: Config) -> tuple[np.ndarray, list]:
    """Read the Sting-SED*.hdf5 file for the given model/subvolume"""

    data = defaultdict(list)
    ids = []

    for sub_volume in config.dirs.sub_volumes:
        full_name = config.print_full_file_name("sed", sub_volume)
        print(f"Reading galaxies data from {full_name}")

        with h5py.File(full_name, "r") as file:
            ids.append(file["id_galaxy_sky"][()])

            for group_name, data_names in config.sed_fields.items():
                group = file[group_name]
                for data_name in data_names:
                    full_data_path = f"{group_name}/{data_name}"
                    data[full_data_path].append(group[data_name][()])

    ids = np.concatenate(ids, axis=0)
    for key in data:
        data[key] = np.concatenate(data[key], axis=1)

    filters = read_filter_names(config)
    data = combine_filters_and_data(filters, data)
    return ids, data
