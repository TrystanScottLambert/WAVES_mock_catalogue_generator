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

def combine_filters_and_data(filter_names: list[str], filter_data: dict) -> dict:
    """
    Takes the data dictionary from the 'read_photometry_data_hdf5' function and combines it with
    the filter names from 'read_filter_names' and creates a new dictionary where each key represents
    each combination of filters (total_ap_dust_FUV_GALEX, or bulge_t_ab_dust_Band8_ALMA).
    """
    # Restructuring the dictionary
    new_data = {}

    for old_key, arrays in filter_data.items():
        components = old_key.split('/')

        for i, filter_name in enumerate(filter_names):
            new_key = f'{components[2]}_{components[1]}_{filter_name}'
            new_data[new_key] = arrays[i]
    return new_data

def read_photometry_data_hdf5(
    model_dir: str, sub_dir: str, fields: dict, sub_volumes: np.ndarray, sed_file: str
) -> tuple[np.ndarray, list]:
    """Read the Sting-SED*.hdf5 file for the given model/subvolume"""

    data = defaultdict(list)
    ids = []

    for subv in sub_volumes:
        full_name = os.path.join(model_dir, sub_dir, f"{sed_file}_{subv:02d}.hdf5")
        print(f"Reading galaxies data from {full_name}")

        with h5py.File(full_name, "r") as file:
            ids.append(file["id_galaxy_sky"][()])

            for group_name, data_names in fields.items():
                group = file[group_name]
                for data_name in data_names:
                    full_data_path = f"{group_name}/{data_name}"
                    data[full_data_path].append(group[data_name][()])

    ids = np.concatenate(ids, axis=0)
    for key in data:
        data[key] = np.concatenate(data[key], axis=1)

    filters = read_filter_names(model_dir, sub_dir, sed_file)
    data = combine_filters_and_data(filters, data)
    return ids, data



if __name__ == '__main__':
    lightcone_dir = '/scratch/pawsey0119/clagos/Stingray/output/medi-SURFS/Shark-TreeFixed-ReincPSO-kappa0p002/deep-optical-final/'
    subvols = np.arange(10)
    subdir = 'split/'
    fields = {'galaxies': ('dec', 'ra', 'zobs',
                           'id_galaxy_sky','sfr_burst','sfr_disk','mstars_bulge','mstars_disk','rstar_bulge_apparent',
                           'rstar_disk_apparent','id_group_sky','dc', 'mvir_hosthalo', 'type')}
    
    sed_fields =  {'SED/ap_dust': ('total', 'bulge_t'), 'SED/ab_dust': ('total', 'bulge_t')}
    sed_file = "Sting-SED-VST-eagle-rr14"

    data = read_lightcone(lightcone_dir, subdir, fields, subvols, 'mock')
    sed_data = read_photometry_data_hdf5(lightcone_dir, subdir, sed_fields, np.arange(3), sed_file)
    filters = read_filter_names(lightcone_dir, subdir, sed_file)