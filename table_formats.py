"""
Formats the data into tables including any calculated properties.
"""

from dataclasses import dataclass

import numpy as np
from astropy.cosmology import FlatLambdaCDM


@dataclass
class DataDescription:
    """
    Dataclass to store the column name, description, and data associated with each column.
    """

    column_name: str
    description: str
    data: np.ndarray


class CalculatedTable:
    """
    Base class for tables that take in inputs and return caclulated properties.
    Used specifically for the groups and galaxies.
    """

    def __init__(
        self, scrapped_dict: dict, description_dict: dict, cosmology: FlatLambdaCDM
    ) -> None:
        """
        Scrapped dict is the scrapped data from the various different Hdf5 files. Keys represent
        the properties and the data arrays are the values.
        """
        for key, value in scrapped_dict.items():
            setattr(self, key, DataDescription(key, description_dict[key], value))
        self.cosmo = cosmology

    def list_columns(self) -> list[str]:
        """
        Lists all available columns, including attributes and calculated properties.
        """
        return [
            name
            for name in dir(self)
            if isinstance(getattr(self, name), DataDescription)
        ]

    def get_column(self, column_name: str) -> DataDescription:
        """
        Retrieves the DataDescription for the given column name.
        """
        return getattr(self, column_name)

    def sample(self, list_of_columns: list[str] = None) -> tuple[dict, dict]:
        """
        Takes a subsample of column names that want to be written and returns the header to be
        written to the file as well as a dictionary that can be combined with sed data or other
        galaxy data read in from elsewhere.

        If list_of_columns is None then the entire data set is generated (default)
        """
        if list_of_columns is None:
            list_of_columns = self.list_columns()

        columns = [self.get_column(col_name) for col_name in list_of_columns]
        header_dict = {column.column_name : column.description for column in columns}
        writeable_dict = {column.column_name : column.data for column in columns}
        return header_dict, writeable_dict


class GalaxyTable(CalculatedTable):
    """
    Manages the galaxy data that we read in from the mock hdf5 files.
    """

    @property
    def log_mstar_total(self) -> DataDescription:
        """
        log of the total stellar mass
        """
        column_name = "log_mstar_total"
        description = "The total stellar mass of the system."
        value = np.log10(
            (self.mstars_bulge.data + self.mstars_disk.data) / self.cosmo.h
        )
        return DataDescription(column_name, description, value)

    @property
    def log_sfr_total(self) -> DataDescription:
        """
        log of the total star formation rate in [...]
        """
        column_name = "log_sfr_total"
        description = "The total star formation rate []"
        value = np.log10(
            (self.sfr_burst.data + self.sfr_disk.data) / self.cosmo.h / 1e9
        )
        return DataDescription(column_name, description, value)

    @property
    def re(self) -> DataDescription:
        """
        What is this? This needs to be written out explicitly and with units.
        """
        column_name = "re"
        description = "Effective radius? "

        total_mass = self.mstars_bulge.data + self.mstars_disk.data
        weighted_re_sums = (
            self.rstar_bulge_apparent.data * self.mstars_bulge.data
            + self.rstar_disk_apparent.data * self.mstars_disk.data
        )

        value = weighted_re_sums / total_mass
        return DataDescription(column_name, description, value)

    @property
    def bt(self) -> DataDescription:
        """
        Bulge-to-Total mass ratio
        """
        column_name = "bt"
        description = "Bulge-to-total mass ratio"
        value = self.mstars_bulge.data / (
            self.mstars_bulge.data + self.mstars_disk.data
        )
        return DataDescription(column_name, description, value)

class GroupTable(CalculatedTable):
    """
    Managing the data that was read from the hdf5 files.
    """
