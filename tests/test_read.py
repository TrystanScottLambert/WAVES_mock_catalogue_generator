"""
Testing the functions in the 'read' module.
"""

import os
import sys
import unittest
import tempfile
import numpy as np
import h5py

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from read import read_lightcone, read_filter_names, combine_filters_and_data


class TestReadLightcone(unittest.TestCase):
    """
    Testing the 'read_lightcone' function.
    """

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.model_dir = self.test_dir.name
        self.sub_dir = "subdir"
        os.makedirs(os.path.join(self.model_dir, self.sub_dir), exist_ok=True)

        # Define the sample fields structure and sub_volumes array
        self.fields = {"galaxies": ["dec", "ra", "zobs"]}
        self.sub_volumes = np.array([0, 1])  # Example subvolume identifiers
        self.file_name = "test_file"

        # Create small HDF5 files in the temporary directory
        for sub_volume in self.sub_volumes:
            file_path = os.path.join(
                self.model_dir, self.sub_dir, f"{self.file_name}_{sub_volume:02d}.hdf5"
            )
            with h5py.File(file_path, "w") as f:
                grp = f.create_group("galaxies")
                grp.create_dataset("dec", data=np.random.rand(10))
                grp.create_dataset("ra", data=np.random.rand(10))
                grp.create_dataset("zobs", data=np.random.rand(10))

    def tearDown(self):
        self.test_dir.cleanup()

    def test_read_lightcone(self):
        """Call the read_lightcone function"""
        data = read_lightcone(
            model_dir=self.model_dir,
            sub_dir=self.sub_dir,
            fields=self.fields,
            sub_volumes=self.sub_volumes,
            file_name=self.file_name,
        )

        # Assert that all specified fields are in the output
        for field in self.fields["galaxies"]:
            self.assertIn(field, data)
            # Check that data arrays are non-empty and have the expected length
            self.assertEqual(
                len(data[field]), 20
            )  # Should be 10 entries per subvolume * 2 sub_volumes
            self.assertTrue(isinstance(data[field], np.ndarray))


class TestReadFilterNames(unittest.TestCase):
    """
    Testing the read_filters function.
    """

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.model_dir = self.test_dir.name
        self.sub_dir = "subdir"
        os.makedirs(os.path.join(self.model_dir, self.sub_dir), exist_ok=True)

        # Create a mock HDF5 file
        self.sed_file = "test_sed_file"
        self.file_path = os.path.join(
            self.model_dir, self.sub_dir, f"{self.sed_file}_00.hdf5"
        )

        # Create an HDF5 file with mock filter names
        with h5py.File(self.file_path, "w") as f:
            filter_names = np.array([b"filter_a", b"filter_b", b"filter_c"], dtype="S")
            f.create_dataset("filters", data=filter_names)

    def tearDown(self):
        # Clean up the temporary directory
        self.test_dir.cleanup()

    def test_read_filter_names(self):
        """
        Testing the read_filter function
        """
        filter_names = read_filter_names(
            model_dir=self.model_dir, sub_dir=self.sub_dir, sed_file=self.sed_file
        )

        expected_filter_names = ["filter_a", "filter_b", "filter_c"]
        self.assertEqual(filter_names, expected_filter_names)


class TestCombineFiltersAndData(unittest.TestCase):
    """
    Testing the combine_filters_and_data function.
    """
    def test_combine_filters_and_data(self):
        """
        Simple test.
        """
        filter_names = ["FUV_GALEX", "NUV_GALEX", "Band8_ALMA"]

        filter_data = {
            "SED/ab_dust/total": [np.array([1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9])],
            "SED/ab_dust/bulge": [np.array([10, 11, 12]), np.array([13, 14, 15]), np.array([16, 17, 18])]
        }

        expected_output = {
            "total_ab_dust_FUV_GALEX": np.array([1, 2, 3]),
            "total_ab_dust_NUV_GALEX": np.array([4, 5, 6]),
            "total_ab_dust_Band8_ALMA": np.array([7, 8, 9]),
            "bulge_ab_dust_FUV_GALEX": np.array([10, 11, 12]),
            "bulge_ab_dust_NUV_GALEX": np.array([13, 14, 15]),
            "bulge_ab_dust_Band8_ALMA": np.array([16, 17, 18])
        }

        result = combine_filters_and_data(filter_names, filter_data)

        for key in expected_output:
            with self.subTest(key=key):
                np.testing.assert_array_equal(result[key], expected_output[key])



if __name__ == "__main__":
    unittest.main()
