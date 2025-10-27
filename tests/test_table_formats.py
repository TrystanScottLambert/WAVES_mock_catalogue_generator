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

    def test_list_columns(self):
        """
        Testing the list columns method.
        """
        correct = ["dec", "ra"]  # alphabetical order
        self.assertListEqual(correct, self.test_table.list_columns())

    def test_get_column(self):
        """
        Testing the get column method.
        """
        correct_data = self.scrapped_dict["ra"]
        self.assertEqual(type(self.test_table.get_column("ra")), DataDescription)
        npt.assert_array_equal(correct_data, self.test_table.get_column("ra").data)
        self.assertEqual(
            "This is the ra.", self.test_table.get_column("ra").description
        )
        self.assertEqual("ra", self.test_table.get_column("ra").column_name)

    def sample(self):
        """
        Testing the sample method.
        """
        # if we pass the scrapped dict that's what we'll get as the writeable dict.
        correct_header_dict = {"ra": "This is the ra.", "dec": "This is the dec."}
        head_dict, write_dict = self.test_table.sample(
            ["dec", "ra"]
        )  # doesn't matter the order.
        for key, item in head_dict.items():
            self.assertEqual(correct_header_dict[key], item)

        for key, item in write_dict.items():
            self.assertEqual(self.scrapped_dict[key], item)


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
        Testing that the properties that are read in.
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

    def test_sample(self):
        """
        Testing that adding the calculated properties doesn't mess up the sampling.
        """
        correct_header = {
            "mstars_bulge": "This is bulge.",
            "log_mstar_total": "The total stellar mass of the system.",
        }
        correct_data = np.log10(
            (self.scrapped_dict["mstars_bulge"] + self.scrapped_dict["mstars_disk"])
            / self.cosmo.h
        )
        correct_write = {
            "mstars_bulge": self.scrapped_dict["mstars_bulge"],
            "log_mstar_total": correct_data,
        }

        head_dict, write_dict = self.test_table.sample(
            ["log_mstar_total", "mstars_bulge"]
        )
        for key, item in head_dict.items():
            self.assertEqual(correct_header[key], item)

        for key, item in write_dict.items():
            npt.assert_array_equal(correct_write[key], item)


if __name__ == "__main__":
    unittest.main()
