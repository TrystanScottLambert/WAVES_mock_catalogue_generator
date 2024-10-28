"""
Creating a mock catalogue from the SHARK runs on pawsey.
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM

from read import read_lightcone, read_photometry_data_hdf5
from write import CatalogueDetails, write_catagloue
from property_dictionaries import GALAXY_PROPERTIES, GROUP_PROPERTIES
from table_formats import GalaxyTable, GroupTable


# Inputs
cosmo = FlatLambdaCDM(Om0=0.3, H0=67.51)
cat_details = CatalogueDetails(
    area=107.889,
    mag_filter="total_ap_dust_Z_VISTA",
    mag_cut=21.2,
    redshift_cut=3,
    version="b0.0.1",
)

GROUP_FIELDS = {
    "groups": (
        "id_group_sky",
        "ra",
        "dec",
        "zobs",
        'subvolume',
        'tile', 
        'id_halo_sam',
        'flag',
        'snapshot',
    )
}

GALAXY_FIELDS = {
    "galaxies": (
        "dec",
        "ra",
        "zobs",
        "id_galaxy_sky",
        "sfr_burst",
        "sfr_disk",
        "mstars_bulge",
        "mstars_disk",
        "rstar_bulge_apparent",
        "rstar_disk_apparent",
        "id_group_sky",
        "dc",
        "mvir_hosthalo",
        "type",
        'tile',
        'subvolume',
        'id_halo_sam',
        'snapshot',
    )
}

GALAXY_PROPS_TO_WRITE = [
    "unique_group_id",
    "ra",
    "dec",
    "zobs",
    'flag',
]

PROPS_TO_WRITE = [
        "id_galaxy_sky",
        "ra",
        "dec",
        "zobs",
        "log_mstar_total",
        "log_sfr_total",
        "bt",
        "re",
        'unique_group_id',
    ]

SED_FIELDS = {"SED/ap_dust": ["total"], "SED/ab_dust": ["total"]}
LIGHTCONE_DIR = "/scratch/pawsey0119/clagos/Stingray/output/medi-SURFS/Shark-TreeFixed-ReincPSO-kappa0p002/deep-optical-final/"
SUB_DIR = "split/"
SED_FILE = "Sting-SED-VST-eagle-rr14"
SUB_VOLUMES = np.arange(2)

def main():
    """
    Main function to to manage scoping.
    """
    # Reading the data from the hdf5 files
    galaxy_data = read_lightcone(LIGHTCONE_DIR, SUB_DIR, GALAXY_FIELDS, SUB_VOLUMES, "mock")
    group_data = read_lightcone(LIGHTCONE_DIR, SUB_DIR, GROUP_FIELDS, SUB_VOLUMES, "mock")
    _, sed_data = read_photometry_data_hdf5(
        LIGHTCONE_DIR, SUB_DIR, SED_FIELDS, SUB_VOLUMES, SED_FILE)

    # filtering the data
    indicies = np.where((sed_data[cat_details.mag_filter] < cat_details.mag_cut)
        & (sed_data[cat_details.mag_filter] > 0))[0]
    sed_data = {key: value[indicies] for key, value, in sed_data.items()}
    galaxy_data = {key: value[indicies] for key, value in galaxy_data.items()}

    # Working out calculated properties
    galaxy_data = GalaxyTable(galaxy_data, GALAXY_PROPERTIES, cosmo)
    group_data = GroupTable(group_data, GROUP_PROPERTIES, cosmo)

    # Writing
    galaxy_header, galaxy_data_to_write = galaxy_data.sample(list_of_columns=PROPS_TO_WRITE)
    group_header, group_data_to_write = group_data.sample(list_of_columns=GALAXY_PROPS_TO_WRITE)
    write_catagloue([galaxy_data_to_write, sed_data], galaxy_header, cat_details, "delete.txt")
    write_catagloue([group_data_to_write], group_header, cat_details, 'test_delete.txt')


if __name__ == "__main__":
    #main()
    galaxy_data = read_lightcone(LIGHTCONE_DIR, SUB_DIR, GALAXY_FIELDS, SUB_VOLUMES, "mock")
    group_data = read_lightcone(LIGHTCONE_DIR, SUB_DIR, GROUP_FIELDS, SUB_VOLUMES, "mock")
    _, sed_data = read_photometry_data_hdf5(
        LIGHTCONE_DIR, SUB_DIR, SED_FIELDS, SUB_VOLUMES, SED_FILE)

    # filtering the data
    indicies = np.where((sed_data[cat_details.mag_filter] < cat_details.mag_cut)
        & (sed_data[cat_details.mag_filter] > 0))[0]
    sed_data = {key: value[indicies] for key, value, in sed_data.items()}
    galaxy_data = {key: value[indicies] for key, value in galaxy_data.items()}

    # Working out calculated properties
    galaxy_data = GalaxyTable(galaxy_data, GALAXY_PROPERTIES, cosmo)
    group_data = GroupTable(group_data, GROUP_PROPERTIES, cosmo)

    # Writing
    galaxy_header, galaxy_data_to_write = galaxy_data.sample(list_of_columns=PROPS_TO_WRITE)
    group_header, group_data_to_write = group_data.sample(list_of_columns=GALAXY_PROPS_TO_WRITE)
    write_catagloue([galaxy_data_to_write, sed_data], galaxy_header, cat_details, "delete.txt")
    write_catagloue([group_data_to_write], group_header, cat_details, 'test_delete.txt')
