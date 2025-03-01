"""
Handles loading in the configuration file and input checking.
"""

import os
import errno
import warnings
from dataclasses import dataclass
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import yaml

from write import CatalogueDetails
from property_dictionaries import GROUP_PROPERTIES, GALAXY_PROPERTIES


@dataclass
class FileStrings:
    """
    Stores all the string information that we need to locate the files.
    """

    lightcone_directory: str
    sub_directory: str
    lightcone_file: str
    sed_file: str
    sub_volumes: np.ndarray


@dataclass
class Config:
    """
    Stores the final input settings.
    """

    cosmo: FlatLambdaCDM
    cat_details: CatalogueDetails
    gal_props_read: dict
    group_props_read: dict
    gal_props_write: list
    group_props_write: list
    sed_fields: dict
    dirs: FileStrings
    outfile_prefix: str

    def __post_init__(self):
        """
        Create the outfile names post init
        """
        self.galaxy_outfile_name = f"{self.outfile_prefix}_gals.parquet"
        self.group_outfile_name = f"{self.outfile_prefix}_groups.parquet"

    def dump_directory(self, file_type: str) -> tuple[str, str, str]:
        """
        Takes either 'sed' or 'mock' and then returns the parent directory and the file name.
        This is so that we can more easily handle the file paths in functions in the read module.
        """
        if file_type == "sed":
            file_name = self.dirs.sed_file
        elif file_type == "mock":
            file_name = self.dirs.lightcone_file
        else:
            raise ValueError(f'file_type must be "sed" or "mock" not "{file_type}"')
        return self.dirs.lightcone_directory, self.dirs.sub_directory, file_name

    def print_full_file_name(self, file_type: str, sub_volume: int) -> str:
        """
        Returns the full path name for the given file type and sub_volume
        """
        parent_dir, sub_dir, file_name = self.dump_directory(file_type)
        return os.path.join(parent_dir, sub_dir, f"{file_name}_{sub_volume:02d}.hdf5")


def load_cosmo(input_params: dict) -> FlatLambdaCDM:
    """
    Reads in the Cosmology parameters and creates an object.
    """
    return FlatLambdaCDM(
        H0=input_params["Cosmology"]["H0"], Om0=input_params["Cosmology"]["Om0"]
    )


def load_cat_details(input_params: dict) -> CatalogueDetails:
    """
    Reads in the catagloue details and creates an object
    """
    cat_details = CatalogueDetails(
        area=input_params["Catalogue_Details"]["area"],
        mag_filter=input_params["Catalogue_Details"]["mag_filter"],
        mag_cut=input_params["Catalogue_Details"]["mag_cut"],
        redshift_cut=input_params["Catalogue_Details"]["redshift_cut"],
        version=input_params["Catalogue_Details"]["version"],
    )
    return cat_details


def remove_duplicates_in_list(possibly_duplicated_list: list) -> list:
    """
    Removes duplicated values in a list but preserves the order. 
    """
    return list(dict.fromkeys(possibly_duplicated_list))


def load_read_properties(input_parameters: dict) -> tuple[dict]:
    """
    Reads in the properties that need to read in galaxy and groups.
    """

    group_fields = {
        "groups": tuple(remove_duplicates_in_list(input_parameters["Properties_To_Read_In"]["groups"]))
    }
    galaxy_fields = {
        "galaxies": tuple(
            remove_duplicates_in_list(input_parameters["Properties_To_Read_In"]["galaxies"]))
    }

    # Checking that these fields actually even exist in the first place.
    bad_group_fields = [
        key for key in group_fields["groups"] if key not in GROUP_PROPERTIES.keys()
    ]
    bad_galaxy_fields = [
        key for key in galaxy_fields["galaxies"] if key not in GALAXY_PROPERTIES.keys()
    ]
    if len(bad_group_fields) > 0 or len(bad_galaxy_fields) > 0:
        raise AttributeError(
            f"{len(bad_group_fields)} bad group field(s) and {len(bad_galaxy_fields)} bad galaxy field(s) found. \nBad group fields: {bad_group_fields} \nBad galaxy fields: {bad_galaxy_fields}"
        )

    # warn user if multiple values are found in the config
    if len(set(input_parameters["Properties_To_Read_In"]["groups"])) != len(
        input_parameters["Properties_To_Read_In"]["groups"]
    ):
        warnings.warn(
            "Repeated values in the config file for group properties to read in. Value will be read in only once. Consider editing config."
        )
    if len(set(input_parameters["Properties_To_Read_In"]["galaxies"])) != len(
        input_parameters["Properties_To_Read_In"]["galaxies"]
    ):
        warnings.warn(
            "Repeated values in the config file for galaxy properties to read in. Value will be read in only once. Consider editing config."
        )

    return group_fields, galaxy_fields


def load_write_properties(input_parameters: dict) -> tuple[list]:
    """
    Reads in the properties that need to be written for the galaxies and groups.
    """
    # We need to check that we can even make these properties. This is so far down the line
    # that its important to do. Recommend creating Write_Properties_YAML straight from the
    group_props = remove_duplicates_in_list(input_parameters["Properties_To_Write"]["groups"])
    gal_props = remove_duplicates_in_list(input_parameters["Properties_To_Write"]["galaxies"])

    # warn if multiple values are found
    if len(set(input_parameters["Properties_To_Write"]["groups"])) != len(
        input_parameters["Properties_To_Read_In"]["groups"]
    ):
        warnings.warn(
            "Repeated values in the config file for group properties to write. Value will be read in only once. Consider editing config."
        )
    if len(set(input_parameters["Properties_To_Write"]["groups"])) != len(
        input_parameters["Properties_To_Read_In"]["groups"]
    ):
        warnings.warn(
            "Repeated values in the config file for group properties to write. Value will be read in only once. Consider editing config."
        )
    return group_props, gal_props


def load_subvolumes(input_parameters: dict) -> np.ndarray:
    """
    Reads in the subvolume parameters and either generates the array from 0 to the given number,
    or converts the passed list as an array
    """
    if isinstance(input_parameters["Sub_Volumes"], list):
        val = np.array(input_parameters["Sub_Volumes"])
    elif isinstance(input_parameters["Sub_Volumes"], (int, float)):
        val = np.arange(input_parameters["Sub_Volumes"])
    else:
        raise ValueError("Sub_Volumes must be either or list or a single integer/float")
    return val


def load_directory_string(input_parameters: dict) -> FileStrings:
    """
    Reads in the values required to find the files and returns a FileStrings Object.
    """
    sub_volumes = load_subvolumes(input_parameters)
    file_strings = FileStrings(
        lightcone_directory=input_parameters["Lightcone_Directory"],
        sub_directory=input_parameters["Sub_Directory"],
        sed_file=input_parameters["SED_file"],
        lightcone_file=input_parameters["Lightcone_file"],
        sub_volumes=sub_volumes,
    )

    # Checking that the files exist.
    sed_files = [
        os.path.join(
            file_strings.lightcone_directory,
            file_strings.sub_directory,
            f"{file_strings.sed_file}_{sub_volume:02d}.hdf5",
        )
        for sub_volume in sub_volumes
    ]

    mock_files = [
        os.path.join(
            file_strings.lightcone_directory,
            file_strings.sub_directory,
            f"{file_strings.lightcone_file}_{sub_volume:02d}.hdf5",
        )
        for sub_volume in sub_volumes
    ]
    for sed_file in sed_files:
        if not os.path.isfile(sed_file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), sed_file)

    for mock_file in mock_files:
        if not os.path.isfile(mock_file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), mock_file)

    return file_strings


def validate_input_file(input_parameters: dict) -> None:
    """
    Checks that the input file has correct values before running.
    """
    required_settings = [
        "Cosmology",
        "Catalogue_Details",
        "Properties_To_Read_In",
        "Properties_To_Write",
        "SED_fields",
        "Lightcone_Directory",
        "Sub_Directory",
        "SED_file",
        "Lightcone_file",
        "Sub_Volumes",
        "Outfile_Prefix",
    ]
    additional = [
        key for key in input_parameters.keys() if key not in required_settings
    ]
    missing = [key for key in required_settings if key not in input_parameters.keys()]
    if len(missing) > 0:
        raise ValueError(
            f"The following keywords are missing from the config file: {missing}"
        )

    if len(additional) > 0:
        warnings.warn(
            f"There are unknown settings in the YAML file and will be ignored: {additional}"
        )


def load_all() -> Config:
    """
    Main function which loads all the settings and returns the 'Config' class.
    """
    with open("config.yml", encoding="utf-8") as file:
        settings = yaml.safe_load(file)

    validate_input_file(
        settings
    )  # checks that all the settings are there in the first place.
    cosmo = load_cosmo(settings)  # checks cosmology is correct.
    file_strings = load_directory_string(settings)  # checks files actually exist.
    group_read_props, gal_read_props = load_read_properties(settings)
    group_write_props, gal_write_props = load_write_properties(settings)
    cat_details = load_cat_details(settings)

    return Config(
        cosmo=cosmo,
        cat_details=cat_details,
        gal_props_read=gal_read_props,
        gal_props_write=gal_write_props,
        group_props_read=group_read_props,
        group_props_write=group_write_props,
        sed_fields=settings["SED_fields"],
        dirs=file_strings,
        outfile_prefix=settings["Outfile_Prefix"],
    )
