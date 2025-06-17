"""
General functions used for building the WAVES mock catalogue.
"""

from collections import defaultdict
import warnings
import glob

import h5py
import numpy as np

from load import Config


def read_spectra(
    file_name: str, match: bool = False, match_ids: np.ndarray[float] = None
) -> np.ndarray[np.ndarray]:
    """
    Reads a spectra file and returns a 2d array where every row is the galaxy the columns are
    the wavelength.

    |sky_id| val1 val2 val3 ....

    The wavelengths can/should be read separately.
    """
    print(f'Getting spectra from: {file_name}')
    spectra = h5py.File(file_name)
    sky_ids = spectra["id_galaxy_sky"][:]

    if match is True:
        sky_id_matches = np.intersect1d(sky_ids, match_ids)
        if len(sky_id_matches) == 0:
            return None
        idx_matches = np.array(
            [np.where(sky_ids == sky_id)[0][0] for sky_id in sky_id_matches]
        )
        all_spectra = spectra["spectra"][:, idx_matches]
        sky_ids = sky_ids[idx_matches]

    else:
        warnings.warn("Not matching the indicies will take a long time.")
        all_spectra = spectra["spectra"][:]

    spectra_table = all_spectra.T  # have every row be a galaxy.
    full_table = np.hstack((sky_ids.reshape(len(sky_ids), 1), spectra_table))
    return full_table


def read_all_spectra(
    directory: str, file_stub: str, matching_ids: np.ndarray[int] = None
) -> np.ndarray[np.ndarray]:
    """
    Reads the spectra files and concatenates the entire thing into a single parquet file.
    """
    spectra_files = np.sort(glob.glob(f"{directory}/{file_stub}*.hdf5"))
    if matching_ids is not None:
        spectra_results = [ read_spectra(file, match=True, match_ids=matching_ids)
                for file in spectra_files]
        # Filter values that had no overlap.
        spectra_results = [result for result in spectra_results if result is not None]
        spectra_table = np.vstack(spectra_results)
    else:
        spectra_table = np.vstack([read_spectra(file) for file in spectra_files])

    return spectra_table


def read_lightcone(config: Config, source_type: str) -> dict[np.ndarray]:
    """
    Reads in the mock data using the group/gal to read values in the config file.
    loops over all subvolumes stored in the config file and adds the columns of the selected
    data together into numpy arrays. This is then returned as a dictionary where every key is 
    the column that needs to be read in and value is the concatenation of all the data across
    all the subvolumes.
    """

    if source_type == "gal":
        fields = config.gal_props_read
    elif source_type == "group":
        fields = config.group_props_read
        print('fields: ', fields)
    else:
        raise ValueError(f'Type must be either "group" or "gal", not "{source_type}"')

    data = defaultdict(list)
    for sub_volume in config.dirs.sub_volumes:
        full_name = config.print_full_file_name("mock", sub_volume, mock_or_sed='mock')
        print(f"Reading data from: {full_name}")
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


def read_photometry_data_hdf5(config: Config) -> tuple[np.ndarray, dict]:
    """Read the Sting-SED*.hdf5 file for the given model/subvolume"""

    data = defaultdict(list)
    ids = []

    for sub_volume in config.dirs.sub_volumes:
        full_name = config.print_full_file_name("sed", sub_volume, mock_or_sed='sed')
        print(f"Reading data from: {full_name}")

        with h5py.File(full_name, "r") as file:
            ids.append(file["id_galaxy_sky"][()])

            for group_name, data_names in config.sed_fields.items():
                group = file[group_name]
                for data_name in data_names:
                    full_data_path = f"{group_name}/{data_name}"
                    data[full_data_path].append(group[data_name][()])

    ids = np.concatenate(ids)
    for key in data:
        data[key] = np.concatenate(data[key], axis=1)

    filters = read_filter_names(config)
    data = combine_filters_and_data(filters, data)
    return ids, data
