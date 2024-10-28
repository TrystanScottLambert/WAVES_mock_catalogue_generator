"""
Tests for the 'write' module.
"""

import tempfile
import unittest

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from write import stamp_preamble, CatalogueDetails, write_catagloue


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


class TestWriteCatalogue(unittest.TestCase):
    """
    Testing the write_catalogue function.
    """

    def setUp(self):

        self.temp_file_name = "temp.txt"

        self.writable_dicts = [
            {"col1": [1, 2, 3], "col2": [4, 5, 6]},
            {"col3": [7, 8, 9]},
        ]
        self.unit_header = "# Units: \n"
        self.cat_details = CatalogueDetails(
            area=500.0, mag_filter="g", mag_cut=22.5, redshift_cut=0.5, version=1.0
        )
        self.delimeter = " "
        with open(HEADER_PREAMBLE_DIR, encoding="utf8") as preamble_file:
            preamble = preamble_file.read()

        self.expected_output = preamble + (
            "# Lightcone Details: \n"
            "# area: 500.0 deg2 \n# g <= 22.5 \n"
            "# redshift < 0.5 \n# VERSION: v1.0 \n"
            "# Units: \n"
            "col1 col2 col3 \n"
            "1 4 7 \n"
            "2 5 8 \n"
            "3 6 9 \n"
        )

    def test_write_catagloue(self):
        """
        Testing writing functionality.
        """
        with open(self.temp_file_name, "w", encoding="utf-8") as temp_file:
            outfile = temp_file.name

            write_catagloue(
                writable_dicts=self.writable_dicts,
                unit_header=self.unit_header,
                cat_details=self.cat_details,
                outfile=outfile,
                delimeter=self.delimeter,
            )

            with open(outfile, "r", encoding="utf-8") as f:
                output = f.read()
            self.assertEqual(output, self.expected_output)

    def tearDown(self):
        os.remove(self.temp_file_name)


if __name__ == "__main__":
    unittest.main()
