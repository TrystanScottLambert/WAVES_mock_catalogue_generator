"""
Testing the read module.
"""

import os
import unittest
from unittest.mock import patch, MagicMock
from astropy.cosmology import FlatLambdaCDM
import h5py
import numpy as np
import numpy.testing as npt

from load import Config, FileStrings
from write import CatalogueDetails
from read import (
    read_lightcone,
    combine_filters_and_data,
    read_filter_names,
    read_spectra,
    read_all_spectra
)


cat_details = CatalogueDetails(1, "g", 20, 0.6, "1.0")
dirs = FileStrings(
    lightcone_directory="light_dir",
    sub_directory="sub_dir",
    lightcone_file="light_file",
    sed_file="sed_file",
    sub_volumes=np.array([0, 1]),
)

test_config = Config(
    FlatLambdaCDM(H0=70, Om0=0.3),
    cat_details,
    {"galaxies": ["ra", "dec"]},
    {"groups": ["ra", "dec"]},
    ["ra", "dec", "angsep"],
    ["ra", "dec", "flag"],
    {"SED/ab_dust": ["total"]},
    dirs=dirs,
    outfile_prefix='test_out',
)

class TestReadLightCone(unittest.TestCase):
    """
    Testing the read_lightcone function.
    """

    def test_input_validation(self):
        """Tests that if the wrong source type is entered an error is raised."""
        with self.assertRaises(ValueError) as context:
            read_lightcone('this doesn"t matter', "not_correct_source_type")
        self.assertEqual(
            str(context.exception),
            'Type must be either "group" or "gal", not "not_correct_source_type"',
        )

    @patch("h5py.File", autospec=True)
    def test_read_lightcone(self, mock_h5py_file):
        """Test reading data with correct input by focusing on mocking h5py.File."""

        # Mock data that should be returned by the HDF5 file structure
        mock_data = {
            "galaxies": {"ra": np.array([1.1, 1.2]), "dec": np.array([2.1, 2.2])}
        }

        # Create a mock for the HDF5 file structure with dataset access
        mock_file = {}
        for group_name, datasets in mock_data.items():
            mock_group = {}
            for dataset_name, data in datasets.items():
                dataset_mock = MagicMock()
                dataset_mock.__getitem__.return_value = (
                    data  # Return the array directly
                )
                mock_group[dataset_name] = dataset_mock
            mock_file[group_name] = mock_group

        # Set the return value of `h5py.File` context manager to this mock file structure
        mock_h5py_file.return_value.__enter__.return_value = mock_file

        # Update test_config to only include a single subvolume
        test_config.dirs.sub_volumes = np.array([0])

        result = read_lightcone(test_config, "gal")
        expected_result = {"ra": np.array([1.1, 1.2]), "dec": np.array([2.1, 2.2])}
        npt.assert_equal(result, expected_result)


class TestCombineFiltersAndData(unittest.TestCase):
    """
    Testing the combine_filters_and_data function
    """

    test_names = ["r", "g", "b"]
    test_data = {
        "SED/ab_dust/total": [
            np.array([20, 21]),
            np.array([19, 18]),
            np.array([22, 23]),
        ]
    }

    def test_run(self):
        "testing that the function returns the correct output"
        correct = {
            "total_ab_dust_r": np.array([20, 21]),
            "total_ab_dust_g": np.array([19, 18]),
            "total_ab_dust_b": np.array([22, 23]),
        }
        val = combine_filters_and_data(self.test_names, self.test_data)
        npt.assert_equal(correct, val)


class TestReadFilterNames(unittest.TestCase):
    """
    Testing the read_filter function
    """

    @patch("h5py.File", autospec=True)
    def test_read_filter_names(self, mock_h5py_file):
        """Test that filter names are read and decoded correctly."""

        # Define mock filter names as byte strings like the structure of the HDF5 file dataset
        mock_filter_names = [b"u", b"g", b"r", b"i", b"z"]

        # Create a mock for the HDF5 file structure with `filters` dataset
        mock_file = MagicMock()
        mock_file.__getitem__.return_value = (
            mock_filter_names  # `f["filters"][:]` returns mock_filter_names
        )
        mock_h5py_file.return_value.__enter__.return_value = {"filters": mock_file}

        result = read_filter_names(test_config)
        expected_result = ["u", "g", "r", "i", "z"]
        self.assertEqual(result, expected_result)


class TestReadSpectra(unittest.TestCase):
    """
    Testing the reading of the spectra files (the spectra that get passed into prospect.)
    """
    def setUp(self):
        """
        Setting up a random hdf5 file which mimics the one on pawsey
        """
        id_galaxy_sky = np.array([1, 2, 3, 4, 5])
        spectra = np.random.rand(10, 5)  # 5 galaxies, each with 10 spectral values
        wavegrid = np.linspace(4000, 7000, 10)  # Example wavelength grid
        with h5py.File("test_spectra_1.hdf5", "w") as f:
            f.create_dataset("id_galaxy_sky", data=id_galaxy_sky)
            f.create_dataset("spectra", data=spectra)
            f.create_dataset("wavegrid", data=wavegrid)

        id_galaxy_sky = np.array([11, 21, 31, 41])
        spectra = np.random.rand(10, 4)  # 5 galaxies, each with 10 spectral values
        with h5py.File("test_spectra_2.hdf5", "w") as f:
            f.create_dataset("id_galaxy_sky", data=id_galaxy_sky)
            f.create_dataset("spectra", data=spectra)
            f.create_dataset("wavegrid", data=wavegrid)

    def test_no_matching(self):
        """
        Testing when no matches are passed.
        """
        table = read_spectra("test_spectra_1.hdf5")
        self.assertEqual(table.shape[0], 5)
        self.assertEqual(table.shape[1], 11) # 10 + 1 including the id

    def test_no_matches(self):
        """
        Testing the case where the matching ids have zero overlap.
        """
        good_list = np.arange(20, 40)
        table = read_spectra("test_spectra_1.hdf5", match=True, match_ids=good_list)
        self.assertIsNone(table)

    def test_matching(self):
        """
        Testing that we can match the catalog to only read galaxies that 
        we want.
        """
        good_list = np.arange(3, 7) # This isn't a subset of the ids.
        good_list = np.append(good_list, np.array([1])) # just to have a non sequential value.
        table = read_spectra("test_spectra_1.hdf5", match=True, match_ids=good_list)

        npt.assert_array_equal(table[:,0], np.array([1, 3, 4, 5]))
        self.assertEqual(table.shape[0], 4)
        self.assertEqual(table.shape[1], 11)

    def test_read_all_no_matching(self):
        """
        Testing reading all the spectra over multiple hdf5 files.
        """
        directory = './'
        file_stub = 'test_spectra_'
        table = read_all_spectra(directory, file_stub)
        self.assertEqual(table.shape[0], 9)
        self.assertEqual(table.shape[1], 11)
        npt.assert_array_equal(table[:, 0], np.array([1, 2, 3, 4, 5, 11, 21, 31, 41]))

    def test_read_all_matches(self):
        """
        Testing the case where we have only some galaxies in hdf5 that we want to match to.
        """
        directory = './'
        file_stub = 'test_spectra_'
        good_list = np.array([2, 4, 11, 31, 102])
        table = read_all_spectra(directory, file_stub, matching_ids=good_list)
        self.assertEqual(table.shape[0], 4)
        self.assertEqual(table.shape[1], 11)
        npt.assert_array_equal(table[:, 0], good_list[:-1])

    def test_read_all_blank(self):
        """
        Testing the case where one of the files has no overlap with matching ids.
        """
        directory = './'
        file_stub = 'test_spectra_'
        good_list = np.array([2, 4, 102]) # no ids from test_spectra_2
        table = read_all_spectra(directory, file_stub, good_list)
        self.assertEqual(table.shape[0], 2)
        self.assertEqual(table.shape[1], 11)
        npt.assert_array_equal(table[:, 0], np.array([2, 4]))

    def tearDown(self):
        os.remove('test_spectra_1.hdf5')
        os.remove('test_spectra_2.hdf5')

if __name__ == "__main__":
    unittest.main()
