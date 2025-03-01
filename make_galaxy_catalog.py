"""
Creating a mock catalogue from the SHARK runs on pawsey.
"""

import numpy as np

from load import load_all, Config
from read import read_lightcone, read_photometry_data_hdf5
from write import write_to_parquet
from property_dictionaries import GALAXY_PROPERTIES, GROUP_PROPERTIES
from table_formats import GalaxyTable, GroupTable

def filter_based_on_mag(config_object: Config, sed_data: dict, galaxy_data: dict) -> dict:
    """
    Based on the magnitude limit stored in the config_object we filter the galaxy and sed data.
    """
    mag_filter = config_object.cat_details.mag_filter
    mag_limit = config_object.cat_details.mag_cut
    cut_idxs = np.where(sed_data[mag_filter] < mag_limit)[0]
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
    galaxy_header, galaxy_data_to_write = galaxy_data.sample(
        list_of_columns=config.gal_props_write
    )
    group_header, group_data_to_write = group_data.sample(
        list_of_columns=config.group_props_write
    )
    write_to_parquet(
        [galaxy_data_to_write, sed_data],
        galaxy_header,
        config.cat_details,
        config.galaxy_outfile_name,
    )
    write_to_parquet(
        [group_data_to_write], group_header, config.cat_details, config.group_outfile_name
    )


if __name__ == "__main__":
    main()
