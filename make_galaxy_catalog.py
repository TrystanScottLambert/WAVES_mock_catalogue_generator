"""
Creating a mock catalogue from the SHARK runs on pawsey.
"""

import h5py
import numpy as np

from group_post_process import add_fof_ids
from load import Config, load_all

from property_dictionaries import GALAXY_PROPERTIES, GROUP_PROPERTIES
from read import read_lightcone, read_photometry_data_hdf5
from table_formats import GalaxyTable, GroupTable
from write import write_to_parquet


def filter_based_on_mag(
    config_object: Config, sed_data: dict, galaxy_data: dict
) -> dict:
    """
    Based on the magnitude limit stored in the config_object we filter the galaxy and sed data.
    """
    mag_filter = config_object.cat_details.mag_filter
    mag_limit = config_object.cat_details.mag_cut
    cut_sed_idxs = np.where(sed_data[mag_filter] < mag_limit)[0]
    cut_gal_idxs = np.where(
        galaxy_data["zobs"] < config_object.cat_details.redshift_cut
    )[0]
    cut_idxs = np.intersect1d(cut_sed_idxs, cut_gal_idxs)

    cut_sed_data = {key: value[cut_idxs] for key, value in sed_data.items()}
    cut_galaxy_data = {key: value[cut_idxs] for key, value in galaxy_data.items()}
    return cut_sed_data, cut_galaxy_data


def main():
    """
    Main function to to manage scoping.
    """
    config = load_all()  # this will also perform the input validation.

    # Reading the data from the hdf5 files
    galaxy_data = read_lightcone(config, "gal")
    group_data = read_lightcone(config, "group")
    _, sed_data = read_photometry_data_hdf5(config)

    sed_data, galaxy_data = filter_based_on_mag(config, sed_data, galaxy_data)

    # Working out calculated properties
    galaxy_data = GalaxyTable(galaxy_data, GALAXY_PROPERTIES, config.cosmo)
    group_data = GroupTable(group_data, GROUP_PROPERTIES, config.cosmo)

    # Writing
    galaxy_write_fields = config.gal_props_write
    group_write_fields = config.group_props_write

    test_sub_volume = config.dirs.sub_volumes[0]
    full_name = config.print_full_file_name("mock", test_sub_volume, mock_or_sed="mock")

    if len(group_write_fields) == 0:
        print("writing all read properties for groups. No selection found in config")
        with h5py.File(full_name, "r") as f:
            group_write_fields = list(f["groups"].keys())

    if len(galaxy_write_fields) == 0:
        print("writing all read properties for galaxies. No selection found in config")
        with h5py.File(full_name, "r") as f:
            galaxy_write_fields = list(f["galaxies"].keys())

    galaxy_header, galaxy_data_to_write = galaxy_data.sample(
        list_of_columns=galaxy_write_fields
    )
    group_header, group_data_to_write = group_data.sample(
        list_of_columns=group_write_fields
    )
    write_to_parquet(
        [galaxy_data_to_write, sed_data],
        galaxy_header,
        config.cat_details,
        config.galaxy_outfile_name,
    )
    write_to_parquet(
        [group_data_to_write],
        group_header,
        config.cat_details,
        config.group_outfile_name,
    )

    # Post-processing for groups
    add_fof_ids(config.galaxy_outfile_name, config.group_outfile_name)


if __name__ == "__main__":
    main()
