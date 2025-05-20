"""
Formats the data into tables including any calculated properties.
"""

from dataclasses import dataclass

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy import interpolate
from scipy.stats import truncnorm


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

    @property
    def bulge_axis_ratio(self) -> DataDescription:
        """
        In Shark, bulges are perfectly spherical, so we sample their axis ratios from 
        https://ui.adsabs.harvard.edu/link_gateway/2023A&A...671A.102E/doi:10.1051/0004-6361/202245042
        """
        column_name = "bulge_axis_ratio"
        description = "Bulge axis ratio"
        loc_bulge, scale_bulge, clip_a_bulge, clip_b_bulge = 0.7, 0.3, 0., 1.
        a_bulge, b_bulge = (clip_a_bulge - loc_bulge) / scale_bulge, (clip_b_bulge - loc_bulge) / scale_bulge
        value = truncnorm.rvs(a_bulge, b_bulge, loc=loc_bulge, scale=scale_bulge, size=len(self.mstars_bulge.data))
        return DataDescription(column_name, description, value)

    @property
    def disk_axis_ratio(self) -> DataDescription:
        """
        Disk axis ratios are approximated by considering the half mass radii as the projected major axes, while the projected minor axes are
        computed with the equation in section 2.2 of Lagos+19 (depends on disk radius and inclination).
        https://ui.adsabs.harvard.edu/link_gateway/2019MNRAS.489.4196L/doi:10.1093/mnras/stz2427
        The value 7.3 comes from the scale height-to-scale length observed relation in local galaxy discs (Kregel, van der Kruit & de Grijs 2002).
        The formula assumes sizes in comoving kpc, while the Shark quantities are in comoving Mpc (hence the 10^3 factor).
        """
        column_name = "disk_axis_ratio"
        description = "Disk axis ratio"
        value = np.sin(self.inclination.data*np.pi/180) * ((self.rstar_disk_intrinsic.data*10**3/self.cosmo.h) - (self.rstar_disk_intrinsic.data*10**3/self.cosmo.h)/7.3) + (self.rstar_disk_intrinsic.data*10**3/self.cosmo.h)/7.3
        return DataDescription(column_name, description, value)

    def bulge_half_light_radius(self, key='half-mass') -> DataDescription:
        """
        This function computes the half-light radius of the bulge from its intrinsic half-mass radius.
        If key = 'half-mass', then we assume the two to be the same.
        If key = 'Suess+19', then we assume the relation in Suess+19 between half-light and half-mass radius.
        Suess+19 compares the two quantities in kpc. The ratio evolves with redshift, being one beyond z=2.2.
        """
        column_name = "bulge_r50"
        description = "Bulge half-light radius"
        if key == 'half-mass':
            value = self.rstar_bulge_intrinsic.data
        elif key == 'Suess+19':
            z_values=np.array([0.25, 0.75, 1.25, 1.75, 2.25])
            ratio_half_light_half_mass_suess19 = np.array([4/3, 2.8/1.9, 1.8/1.4, 1.6/1.3, 1.])
            f_interp = interpolate.interp1d(z_values, ratio_half_light_half_mass_suess19, kind='linear', fill_value=(ratio_half_light_half_mass_suess19[0],ratio_half_light_half_mass_suess19[-1]), bounds_error=False)
            ratio_half_light_half_mass = f_interp(self.zobs.data)
            value = self.rstar_bulge_intrinsic.data * ratio_half_light_half_mass
        else:
            raise KeyError
        return DataDescription(column_name, description, value)

    def disk_half_light_radius(self, key='half-mass') -> DataDescription:
        """
        This function computes the half-light radius of the disk from its intrinsic half-mass radius.
        If key = 'half-mass', then we assume the two to be the same.
        If key = 'Suess+19', then we assume the relation in Suess+19 between half-light and half-mass radius.
        Suess+19 compares the two quantities in kpc. The ratio evolves with redshift, being roughly one beyond z=2.2.
        """
        column_name = "disk_r50"
        description = "Disk half-light radius"
        if key == 'half-mass':
            value = self.rstar_disk_intrinsic.data
        elif key == 'Suess+19':
            z_values=np.array([0.25, 0.75, 1.25, 1.75, 2.25])
            ratio_half_light_half_mass_suess19 = np.array([6/4, 5/3.3, 4.3/2.9, 3.8/3, 3.1/2.9])
            f_interp = interpolate.interp1d(z_values, ratio_half_light_half_mass_suess19, kind='linear', fill_value=(ratio_half_light_half_mass_suess19[0],ratio_half_light_half_mass_suess19[-1]), bounds_error=False)
            ratio_half_light_half_mass = f_interp(self.zobs.data)
            value = self.rstar_disk_intrinsic.data * ratio_half_light_half_mass
        else:
            raise KeyError
        return DataDescription(column_name, description, value)

class GroupTable(CalculatedTable):
    """
    Managing the data that was read from the hdf5 files.
    """
