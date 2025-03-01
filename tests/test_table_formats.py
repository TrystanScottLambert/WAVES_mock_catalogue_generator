"""
Testing the table_formats.py module.
"""

import unittest
import numpy as np
import numpy.testing as npt
from astropy.cosmology import FlatLambdaCDM

from table_formats import CalculatedTable, DataDescription, GalaxyTable


class TestCaclulatedTable(unittest.TestCase):
    """
    Testing the CalculatedTable base class.
    """

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    scrapped_dict = {
        "ra": np.arange(5),
        "dec": np.arange(5) - 10,
    }

    description_dict = {"ra": "This is the ra.", "dec": "This is the dec."}
    test_table = CalculatedTable(scrapped_dict, description_dict, cosmo)

    def test_attributes_created_correctly(self):
        """
        Testing that the correct attributes are being passed and set correctly.
        The calculated table should create ra and dec values as attributes of the table.
        These should also be DataDescription dataclasses which include the descriptions.
        """
        self.assertEqual(type(self.test_table.ra), DataDescription)
        self.assertEqual(type(self.test_table.dec), DataDescription)
        npt.assert_equal(self.test_table.ra.data, self.scrapped_dict["ra"])
        npt.assert_equal(self.test_table.dec.data, self.scrapped_dict["dec"])
        self.assertEqual(self.test_table.ra.column_name, "ra")
        self.assertEqual(self.test_table.dec.column_name, "dec")
        self.assertEqual(self.test_table.ra.description, "This is the ra.")
        self.assertEqual(self.test_table.dec.description, "This is the dec.")

    def test_cosmo_correct(self):
        """
        Testing that the cosmology that has been put in is correct.
        """
        read_cosmo = self.test_table.cosmo.luminosity_distance(np.arange(1000))
        actual = self.cosmo.luminosity_distance(np.arange(1000))
        npt.assert_array_equal(read_cosmo, actual)


class TestGalaxyTable(unittest.TestCase):
    """
    Testing the Galaxy Table which inherrits the CalculatedTable.
    """

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    scrapped_dict = {
        "mstars_bulge": np.arange(1, 6),
        "mstars_disk": np.arange(1, 6) + 10,
    }

    description_dict = {
        "mstars_bulge": "This is bulge.",
        "mstars_disk": "This is disk.",
    }
    test_table = GalaxyTable(scrapped_dict, description_dict, cosmo)

    def test_props_readin(self):
        """
        Testing that the calcualted properties are generated but also the read in properties.
        """
        self.assertEqual(type(self.test_table.mstars_bulge), DataDescription)
        self.assertEqual(type(self.test_table.mstars_disk), DataDescription)
        npt.assert_equal(
            self.test_table.mstars_bulge.data, self.scrapped_dict["mstars_bulge"]
        )
        npt.assert_equal(
            self.test_table.mstars_disk.data, self.scrapped_dict["mstars_disk"]
        )
        self.assertEqual(self.test_table.mstars_bulge.column_name, "mstars_bulge")
        self.assertEqual(self.test_table.mstars_disk.column_name, "mstars_disk")
        self.assertEqual(self.test_table.mstars_bulge.description, "This is bulge.")
        self.assertEqual(self.test_table.mstars_disk.description, "This is disk.")

    def test_calc_props_readin(self):
        """
        Testing that the calculated properties are being done correctly.
        """
        correct_data = np.log10(
            (self.scrapped_dict["mstars_bulge"] + self.scrapped_dict["mstars_disk"])
            / self.cosmo.h
        )
        self.assertEqual(type(self.test_table.log_mstar_total), DataDescription)
        npt.assert_equal(self.test_table.log_mstar_total.data, correct_data)
        self.assertEqual(self.test_table.log_mstar_total.column_name, "log_mstar_total")
        self.assertEqual(
            self.test_table.log_mstar_total.description,
            "The total stellar mass of the system.",
        )


if __name__ == "__main__":
    unittest.main()
