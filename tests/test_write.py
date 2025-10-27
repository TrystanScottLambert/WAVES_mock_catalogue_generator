"""
Tests for the 'write' module.
"""

import tempfile
import unittest
import sys
import os

import numpy as np
import numpy.testing as npt
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from write import (
    stamp_preamble,
    CatalogueDetails,
    write_spectra_table_to_parquet,
    write_to_parquet,
)


HEADER_PREAMBLE_DIR = "header_preamble/header_preamble_v1.head"


class TestStampPreamble(unittest.TestCase):
    """
    Testing the stamp_preamble function which is meant to add the preamble to file objects.
    """

    temp_file_name = "test.txt"

    def setUp(self):
        """Setup: Create a temporary file before each test"""
        open(self.temp_file_name, "a", encoding="utf-8").close()

    def tearDown(self):
        """Teardown: Close and delete the temporary file after each test"""
        os.remove(self.temp_file_name)

    def test_stamp_preamble(self):
        """Call the function to stamp the preamble on the temp file"""
        with open(self.temp_file_name, "a", encoding="utf-8") as temp_file:
            stamp_preamble(temp_file)

        # Read the contents of both the stamped temp file and the original preamble file
        with open(self.temp_file_name, encoding="utf-8") as stamped_file, open(
            HEADER_PREAMBLE_DIR, encoding="utf-8"
        ) as preamble_file:
            stamped_content = stamped_file.read()
            preamble_content = preamble_file.read()

        self.assertEqual(stamped_content, preamble_content)


class TestCatalogueDetails(unittest.TestCase):
    """
    Testing the CatalogueDetails class and methods.
    """

    def setUp(self):
        # Create an instance of CatalogueDetails with example values
        self.catalogue = CatalogueDetails(
            area=500.0, mag_filter="g", mag_cut=22.5, redshift_cut=0.5, version=1.0
        )
        # Expected string output
        self.expected_output = (
            "# Lightcone Details: \n"
            "# area: 500.0 deg2 \n"
            "# g <= 22.5 \n"
            "# redshift < 0.5 \n"
            "# VERSION: v1.0 \n"
        )

    def test_stamp_details(self):
        """Testing the stamp_details method in the class."""
        with tempfile.NamedTemporaryFile(delete=False, mode="w+") as temp_file:
            self.catalogue.stamp_details(temp_file)
            temp_file.flush()
            temp_file.seek(0)
            stamped_content = temp_file.read()
        self.assertEqual(stamped_content, self.expected_output)


class TestWritingSpectra(unittest.TestCase):
    """
    Testing the write functionailty of the spectra files.
    """

    def test_write(self):
        """
        creating a fake 2d array with values and ensureing it's written.
        """
        wavelength = np.array([11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
        ids = np.array([0, 1, 2, 3]).reshape(4, 1)
        data = np.random.rand(4, 10)
        table = np.hstack((ids, data))
        correct_column_headers = np.append(
            np.array(["id_galaxy_sky"]), np.arange(11, 21).astype(int).astype(str)
        )
        write_spectra_table_to_parquet(table, wavelength, "test.parquet")
        test_df = pd.read_parquet("test.parquet", engine="pyarrow")
        df_ids = np.array(test_df["id_galaxy_sky"])
        npt.assert_array_equal(df_ids, np.array([0, 1, 2, 3]))
        self.assertEqual(len(test_df), 4)
        npt.assert_array_equal(np.array(test_df.columns), correct_column_headers)
        npt.assert_array_equal(np.array(test_df["id_galaxy_sky"]), ids.reshape(4))
        os.remove("test.parquet")


class TestWriteToParquet(unittest.TestCase):
    """
    Testing the write_to_parquet function.
    """

    gal_dict = {
        "ra": np.random.random(10),
        "dec": np.random.random(10),
        "z": np.random.random(10),
    }
    sed_dict = {"b_ab_dust": np.random.random(10), "b_ap_dust": np.random.random(10)}

    write_dicts = [gal_dict, sed_dict]
    unit_header = {"ra": "This is RA", "dec": "This is Dec", "z": "This is Redshift"}
    cat_details = CatalogueDetails(1, "g", 20, 0.6, "1.0")

    outfile = "test_write.parquet"

    def test_write(self):
        """
        Testing that the writing process does not scramble things.
        """
        # write the parquet file and we will read it back in and see that it is correct.
        write_to_parquet(
            self.write_dicts, self.unit_header, self.cat_details, self.outfile
        )

        test_df = pd.read_parquet(self.outfile)
        npt.assert_array_equal(np.array(test_df["ra"]), self.gal_dict["ra"])
        npt.assert_array_equal(np.array(test_df["dec"]), self.gal_dict["dec"])
        npt.assert_array_equal(np.array(test_df["z"]), self.gal_dict["z"])
        npt.assert_array_equal(
            np.array(test_df["b_ab_dust"]), self.sed_dict["b_ab_dust"]
        )
        npt.assert_array_equal(
            np.array(test_df["b_ap_dust"]), self.sed_dict["b_ap_dust"]
        )

    def tearDown(self):
        os.remove(self.outfile)


if __name__ == "__main__":
    unittest.main()
