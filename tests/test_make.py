"""
Tesing the functions in the make_galaxy_catalog.py module.
"""

import unittest
import numpy as np
import numpy.testing as npt
from astropy.cosmology import FlatLambdaCDM

from make_galaxy_catalog import filter_based_on_mag
from load import FileStrings, Config
from write import CatalogueDetails

cat_details = CatalogueDetails(1, "r_ap_dust", 9, 0.6, "1.0")
dirs = FileStrings(
    lightcone_directory="./",
    sub_directory="",
    lightcone_file="test_example",
    sed_file="test_sed_example",
    sub_volumes=np.array([0, 1]),
)

test_config = Config(
    FlatLambdaCDM(H0=70, Om0=0.3),
    cat_details,
    {"galaxies": ["ra", "dec", "id_galaxy_sky"]},
    {"groups": ["ra", "id_group_sky"]},
    ["ra", "dec", "angsep"],
    ["ra", "dec", "flag"],
    {"SED/ab_dust": ["total"], "SED/ap_dust": ["total"]},
    dirs=dirs,
    outfile_prefix="test_out",
)


class TestFilteringMagnitudes(unittest.TestCase):
    """
    Testing the filter_based_on_mag function.
    """

    test_sed_data = {
        "b_ap_dust": np.array([2, 3, 4]),
        "b_ab_dust": np.array([-2, -3, -4]),
        "r_ap_dust": np.array([8, 9, 7]),
        "r_ab_dust": np.array([-7, -8, -9]),
    }

    test_galax_data = {
        "id": np.arange(3),
        "dec": np.arange(3) - 20,
        "ra": np.arange(3) + 100,
        "zobs": np.repeat(0.2, 3),
    }

    def test_simple_case(self):
        """
        Testing that when we pass config with a cut < r_ap_dust < 9 all the things
        are filtered correctly and all igore the last two numbres.
        """
        correct_sed = {
            "b_ap_dust": np.array([2, 4]),
            "b_ab_dust": np.array([-2, -4]),
            "r_ap_dust": np.array([8, 7]),
            "r_ab_dust": np.array([-7, -9]),
        }

        correct_galaxy = {
            "id": np.array([0, 2]),
            "dec": np.array([-20, -18]),
            "ra": np.array([100, 102]),
        }

        val_sed, val_gal = filter_based_on_mag(
            test_config, self.test_sed_data, self.test_galax_data
        )
        for key, item in correct_sed.items():
            npt.assert_array_equal(val_sed[key], item)

        for key, item in correct_galaxy.items():
            npt.assert_array_equal(val_gal[key], item)


if __name__ == "__main__":
    unittest.main()
