"""
Creating a mock catalogue from the SHARK runs on pawsey.
"""

import numpy as np

from load import load_all
from read import read_lightcone, read_photometry_data_hdf5
from write import write_to_parquet
from property_dictionaries import GALAXY_PROPERTIES, GROUP_PROPERTIES
from table_formats import GalaxyTable, GroupTable



def main():
    """
    Main function to to manage scoping.
    """
    config = load_all()  # this will also perform the input validation.

    # Reading the data from the hdf5 files
    galaxy_data = read_lightcone(config, "gal")
    group_data = read_lightcone(config, "group")
    _, sed_data = read_photometry_data_hdf5(config)

    # filtering the data
    cat_details = config.cat_details
    indicies = np.where(
        (sed_data[cat_details.mag_filter] < cat_details.mag_cut)
        & (sed_data[cat_details.mag_filter] > 0)
    )[0]
    sed_data = {key: value[indicies] for key, value, in sed_data.items()}
    galaxy_data = {key: value[indicies] for key, value in galaxy_data.items()}

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
        cat_details,
        config.galaxy_outfile_name,
    )
    write_to_parquet(
        [group_data_to_write], group_header, cat_details, config.group_outfile_name
    )


if __name__ == "__main__":
    main()
