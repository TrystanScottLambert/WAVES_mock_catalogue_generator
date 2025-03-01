"""
Testing the load module.
"""

import os
import warnings
from unittest.mock import patch
import unittest
import numpy as np
import numpy.testing as npt
from astropy.cosmology import FlatLambdaCDM
import yaml

from write import CatalogueDetails
from load import (
    load_cosmo,
    load_cat_details,
    load_read_properties,
    load_write_properties,
    load_subvolumes,
    load_directory_string,
    validate_input_file,
    load_all,
    remove_duplicates_in_list,
    Config,
    FileStrings
)

class TestFileString(unittest.TestCase):
    """
    Testing the FileStrings class. Simple data class that only requires testing reading.
    """
    def test_read(self):
        """
        Testing that everything reads in correctly.

        Very basic tests that will only break if the names of the fields change later.
        """
        test_class = FileStrings('light_dir', 'sub_dir', 'light_file', 'sed_file', np.arange(2))
        self.assertEqual(test_class.lightcone_directory, 'light_dir')
        self.assertEqual(test_class.sub_directory, 'sub_dir')
        self.assertEqual(test_class.lightcone_file, 'light_file')
        self.assertEqual(test_class.sed_file, 'sed_file')
        npt.assert_equal(test_class.sub_volumes, np.arange(2))

class TestConfig(unittest.TestCase):
    """
    Testing the Config  data class
    """

    test_config_class = Config(
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3),
            cat_details=CatalogueDetails(1, 'g', 20, 0.1, '1.0'),
            group_props_read={'groups':['ra', 'dec']},
            gal_props_read= {'galaxies': ['zobs']},
            gal_props_write=['zobs', 'zcos'],
            group_props_write=['angsep'],
            sed_fields={'SED/dust': ['total']},
            dirs = FileStrings('light_dir', 'sub_dir', 'light_file', 'sed_file', np.arange(2)),
            outfile_prefix='test'
        )

    def test_read(self):
        """
        Testing that initilizing the class works. This test will only fail if fields are renamed.
        """
        # Testing Cosmology
        self.assertIsInstance(self.test_config_class.cosmo, FlatLambdaCDM)
        self.assertEqual(self.test_config_class.cosmo.H0.value, 70)
        self.assertEqual(self.test_config_class.cosmo.Om0, 0.3)

        # Testing Cat-Details
        self.assertIsInstance(self.test_config_class.cat_details, CatalogueDetails)
        self.assertEqual(self.test_config_class.cat_details.area, 1)
        self.assertEqual(self.test_config_class.cat_details.mag_filter, 'g')
        self.assertEqual(self.test_config_class.cat_details.mag_cut, 20)
        self.assertEqual(self.test_config_class.cat_details.redshift_cut, 0.1)
        self.assertEqual(self.test_config_class.cat_details.version, '1.0')

        # Testing groups_props_read and gal props read
        self.assertIsInstance(self.test_config_class.group_props_read, dict)
        self.assertIsInstance(self.test_config_class.gal_props_read, dict)
        self.assertEqual(self.test_config_class.group_props_read['groups'], ['ra', 'dec'])
        self.assertEqual(self.test_config_class.gal_props_read['galaxies'], ['zobs'])

        # Testing group_prop_write and gal_prop_write
        self.assertIsInstance(self.test_config_class.group_props_write, list)
        self.assertIsInstance(self.test_config_class.gal_props_write, list)
        self.assertEqual(self.test_config_class.gal_props_write, ['zobs', 'zcos'])
        self.assertEqual(self.test_config_class.group_props_write, ['angsep'])

        # Testing sed fields
        self.assertIsInstance(self.test_config_class.sed_fields, dict)
        self.assertEqual(self.test_config_class.sed_fields['SED/dust'], ['total'])

        # Testing directory strings. This is already tested above so we just need to test instance.
        self.assertIsInstance(self.test_config_class.dirs, FileStrings)

    def test_dump_directory_names(self):
        """
        Testing the dump directory method of the Config class
        Specifically if the file returns the correct sed or mock file when promted.
        """
        # Testing SED
        parent_dir, sub_dir, file = self.test_config_class.dump_directory('sed')
        self.assertEqual(parent_dir, 'light_dir')
        self.assertEqual(sub_dir, 'sub_dir')
        self.assertEqual(file, 'sed_file')

        # Testing lightcone
        _, _, file = self.test_config_class.dump_directory('mock')
        self.assertEqual(file, 'light_file')

    def test_dump_directory_validation(self):
        """
        Testing that a ValueError is rasied when we don't use 'mock' or 'sed'
        """
        with self.assertRaises(ValueError) as context:
            self.test_config_class.dump_directory('garbage_input')

        correct_error_string = 'file_type must be "sed" or "mock" not "garbage_input"'
        self.assertEqual(str(context.exception), correct_error_string)

    def test_print_full_file_name(self):
        """
        Testing the print full file name method in the Config class
        """
        val = self.test_config_class.print_full_file_name('mock', 0)
        required = os.path.join('light_dir', 'sub_dir', 'light_file_00.hdf5')
        self.assertEqual(val, required)


class TestLoadingFunctions(unittest.TestCase):
    """
    Testing all the loading functions in the load module. 
    """

    with open('tests/config_testcase.yml', encoding='utf-8') as file:
        settings = yaml.safe_load(file)

    def test_load_cosmology(self):
        """
        Testing the load cosmology function. 
        """
        cosmo = load_cosmo(self.settings)
        self.assertIsInstance(cosmo, FlatLambdaCDM)
        self.assertEqual(cosmo.H0.value, 67.51)
        self.assertEqual(cosmo.Om0, 0.3)

    def test_load_cat_details(self):
        """
        Testing that load_cat_details function.
        """
        cat_details = CatalogueDetails(107.889, 'total_ap_dust_Z_VISTA', 21.2, 3, 'a0.0.1')
        val = load_cat_details(self.settings)
        self.assertIsInstance(val, CatalogueDetails)
        self.assertEqual(val.area, cat_details.area)
        self.assertEqual(val.mag_cut, cat_details.mag_cut)
        self.assertEqual(val.mag_filter, cat_details.mag_filter)
        self.assertEqual(val.redshift_cut, cat_details.redshift_cut)
        self.assertEqual(val.version, cat_details.version)

    def test_load_read_properties(self):
        """
        Testing the load_read_properties function
        """
        # reading in correctly.
        real_gal_dict = {'galaxies': ('ra', 'dec')} # alphabetical
        real_group_dict = {'groups': ('id_group_sky', 'zobs')}
        val_groups, val_gals = load_read_properties(self.settings)
        self.assertIsInstance(val_groups, dict)
        self.assertIsInstance(val_gals, dict)
        self.assertDictEqual(real_gal_dict, val_gals)
        self.assertDictEqual(real_group_dict, val_groups)

        # test assertion is working correctly
        broken_settings = self.settings
        broken_settings['Properties_To_Read_In']['galaxies'].append('test_gals')
        broken_settings['Properties_To_Read_In']['groups'].append('test_groups')
        with self.assertRaises(AttributeError) as context:
            load_read_properties(broken_settings)
        correct_error_string = "1 bad group field(s) and 1 bad galaxy field(s) found. \nBad group fields: ['test_groups'] \nBad galaxy fields: ['test_gals']"
        self.assertEqual(str(context.exception), correct_error_string)

    def test_load_write_properties(self):
        """
        Testing the load_write_properties
        """
        group_val, gal_val = load_write_properties(self.settings)
        gal_correct = ['id_galaxy_sky', 'ra','dec']
        group_correct = ['unique_group_id', 'flag']
        self.assertEqual(group_val, group_correct)
        self.assertEqual(gal_val, gal_correct)

    def test_load_sub_volumes(self):
        """
        Testing the load sub_volumes function.
        """
        # Testing working for list of sub volumes
        test_settings = self.settings
        test_settings['Sub_Volumes']  = [1, 3, 4, 5]
        correct_answer = np.array([1, 3, 4, 5])
        npt.assert_equal(load_subvolumes(test_settings), correct_answer)

        # Testing working for a single value of sub volumes
        test_settings['Sub_Volumes'] = 64
        correct_answer = np.arange(64)
        npt.assert_equal(correct_answer, load_subvolumes(test_settings))

        # Test Validation
        test_settings['Sub_Volumes'] = 'garbage_input'
        with self.assertRaises(ValueError) as context:
            load_subvolumes(test_settings)
        correct_error_message = "Sub_Volumes must be either or list or a single integer/float"
        self.assertEqual(correct_error_message, str(context.exception))


class TestLoadDirectoryString(unittest.TestCase):
    """
    Testing the load_directory_string function
    """
    def setUp(self):
        self.input_params = {
            "Lightcone_Directory": "/fake/path/lightcone",
            "Sub_Directory": "subdir",
            "SED_file": "sed_file",
            "Lightcone_file": "mock_file",
            "Sub_Volumes": [0, 1, 2]  # Test with a few sub-volumes
        }

    @patch("os.path.isfile")
    def test_load_directory_string_files_exist(self, mock_isfile):
        """
        Test working when files are there.
        """
        # Mock `os.path.isfile` to return True for all file checks
        mock_isfile.return_value = True

        # Call the function and check it does not raise an exception
        file_strings = load_directory_string(self.input_params)

        # Validate that FileStrings object is created with expected values
        self.assertIsInstance(file_strings, FileStrings)
        self.assertEqual(file_strings.lightcone_directory, "/fake/path/lightcone")
        self.assertEqual(file_strings.sub_directory, "subdir")
        self.assertEqual(file_strings.sed_file, "sed_file")
        self.assertEqual(file_strings.lightcone_file, "mock_file")
        np.testing.assert_array_equal(file_strings.sub_volumes, np.array([0, 1, 2]))

        # Check that `os.path.isfile` was called for each expected file
        expected_sed_files = [
            os.path.join("/fake/path/lightcone", "subdir", f"sed_file_{i:02d}.hdf5")
            for i in [0, 1, 2]
        ]
        expected_mock_files = [
            os.path.join("/fake/path/lightcone", "subdir", f"mock_file_{i:02d}.hdf5")
            for i in [0, 1, 2]
        ]

        # Check if mock_isfile was called with expected file paths
        expected_calls = [((file,),) for file in (expected_sed_files + expected_mock_files)]
        mock_isfile.assert_has_calls(expected_calls, any_order=True)

    @patch("os.path.isfile")
    def test_load_directory_string_files_missing(self, mock_isfile):
        """
        Test working when files are not there.
        """
        # Mock `os.path.isfile` to return False for all files, simulating missing files
        mock_isfile.return_value = False

        # Verify that FileNotFoundError is raised
        with self.assertRaises(FileNotFoundError):
            load_directory_string(self.input_params) 


class TestValidateInputFile(unittest.TestCase):
    """
    Testing the validate_input_file function
    """
    def setUp(self):
        # Define a base set of input parameters with all required settings
        self.valid_input_parameters = {
            "Cosmology": "Planck18",
            "Catalogue_Details": {},
            "Properties_To_Read_In": [],
            "Properties_To_Write": [],
            "SED_fields": [],
            "Lightcone_Directory": "/path/to/lightcone",
            "Sub_Directory": "subdir",
            "SED_file": "sed_file",
            "Lightcone_file": "lightcone_file",
            "Sub_Volumes": [0, 1, 2],
            "Outfile_Prefix": "test",
        }

    def test_validate_input_file_all_settings_present(self):
        """
        Test running correctly when everything is correct.
        """
        # Test with all required settings; no error or warning should be raised
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")  # Catch all warnings
            validate_input_file(self.valid_input_parameters)
            # Check that no warnings are raised
            self.assertEqual(len(w), 0, "Unexpected warnings raised")

    def test_validate_input_file_missing_required_settings(self):
        """
        Test crashes if a required value is missing from the config file.
        """
        # Remove a required setting and test for ValueError
        missing_setting_input = self.valid_input_parameters.copy()
        missing_setting_input.pop("SED_file")  # Remove a required setting

        with self.assertRaises(ValueError) as cm:
            validate_input_file(missing_setting_input)

        # Check if the error message contains the missing key
        self.assertIn("SED_file", str(cm.exception))
        self.assertIn("missing from the config file", str(cm.exception))

    def test_validate_input_file_with_additional_keys(self):
        """
        Test warning if some random value appears that isn't required.
        """
        # Add an extra, unrecognized key to the input parameters
        extra_key_input = self.valid_input_parameters.copy()
        extra_key_input["Unknown_Setting"] = "extra_value"

        # Check for warning about unknown settings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")  # Catch all warnings
            validate_input_file(extra_key_input)

            # Verify that a warning was raised
            self.assertEqual(len(w), 1, "Expected a warning but none was raised")
            self.assertTrue(
                issubclass(w[-1].category, UserWarning),
                "Expected a UserWarning for additional settings",
            )
            # Check if the warning message mentions the extra key
            self.assertIn("Unknown_Setting", str(w[-1].message))
            self.assertIn("will be ignored", str(w[-1].message))

    def test_validate_input_file_missing_and_additional_keys(self):
        """
        Test with both missing and additional keys
        """
        complex_input = self.valid_input_parameters.copy()
        complex_input.pop("Sub_Directory")  # Remove a required setting
        complex_input["Extra_Key"] = "extra_value"  # Add an unrecognized key

        # Check that a ValueError is raised for missing keys, warning is ignored
        with self.assertRaises(ValueError) as cm:
            validate_input_file(complex_input)

        self.assertIn("Sub_Directory", str(cm.exception))
        self.assertIn("missing from the config file", str(cm.exception))


class TestLoadAllFunction(unittest.TestCase):
    """
    Testing the load_all function.
    """
    @patch("load.validate_input_file")
    @patch("load.load_cosmo")
    @patch("load.load_directory_string")
    @patch("load.load_read_properties")
    @patch("load.load_write_properties")
    @patch("load.load_cat_details")
    @patch("builtins.open")
    @patch("yaml.safe_load")
    def test_load_all(
        self, mock_yaml_load, mock_open, mock_cat_details, 
        mock_load_write_props, mock_load_read_props, 
        mock_load_dir, mock_load_cosmo, mock_validate
    ):
        # Mock YAML content
        mock_yaml_load.return_value = {
            "Cosmology": {"H0": 70, "Om0": 0.3},
            "Catalogue_Details": {"area": 100, "mag_filter": "r", "mag_cut": 24.5, "redshift_cut": 1.0, "version": "1.0"},
            "Properties_To_Read_In": {"groups": ["group_id"], "galaxies": ["galaxy_id"]},
            "Properties_To_Write": {"groups": ["group_mass"], "galaxies": ["galaxy_luminosity"]},
            "SED_fields": {"field1": "data1"},
            "Lightcone_Directory": "/path/to/lightcone/",
            "Sub_Directory": "/path/to/sub/",
            "SED_file": "sed_file",
            "Lightcone_file": "lightcone_file",
            "Sub_Volumes": [0, 1, 2],
            "Outfile_Prefix": "test",
        }

        # Mock returns for dependent functions
        mock_load_cosmo.return_value = FlatLambdaCDM(H0=70, Om0=0.3)
        mock_load_dir.return_value = FileStrings(
            lightcone_directory="/path/to/lightcone/",
            sub_directory="/path/to/sub/",
            sed_file="sed_file",
            lightcone_file="lightcone_file",
            sub_volumes=[0, 1, 2]
        )
        mock_load_read_props.return_value = ({"groups": ["group_id"]}, {"galaxies": ["galaxy_id"]})
        mock_load_write_props.return_value = (["group_mass"], ["galaxy_luminosity"])
        mock_cat_details.return_value = CatalogueDetails(
            area=100, mag_filter="r", mag_cut=24.5, redshift_cut=1.0, version="1.0"
        )

        # Run the function under test
        config = load_all()

        # Assertions
        self.assertIsInstance(config, Config)
        self.assertEqual(config.cosmo.H0.value, 70)
        self.assertEqual(config.cosmo.Om0, 0.3)
        self.assertEqual(config.cat_details.area, 100)
        self.assertEqual(config.gal_props_read, {"galaxies": ["galaxy_id"]})
        self.assertEqual(config.group_props_read, {"groups": ["group_id"]})
        self.assertEqual(config.gal_props_write, ["galaxy_luminosity"])
        self.assertEqual(config.group_props_write, ["group_mass"])
        self.assertEqual(config.sed_fields, {"field1": "data1"})
        self.assertEqual(config.dirs.lightcone_directory, "/path/to/lightcone/")
        self.assertEqual(config.dirs.sub_directory, "/path/to/sub/")
        self.assertEqual(config.dirs.sed_file, "sed_file")
        self.assertEqual(config.dirs.lightcone_file, "lightcone_file")
        self.assertListEqual(config.dirs.sub_volumes, [0, 1, 2])


class TestRemoveDuplicatedList(unittest.TestCase):
    """
    Testing the _remove_duplicates_in_list function.
    """

    def test_list_no_dupes(self):
        """
        Testing when the list doesn't have any duplicate values.
        """
        test_list = [5, 6, 7, 8, 1, 3]
        val_list = remove_duplicates_in_list(test_list)
        self.assertListEqual(val_list, test_list)

    def test_list_with_dupes(self):
        """
        Testing that list works with duplicated values.
        """
        test_list = [5, 6, 6, 7, 1, 2, 3, 2]
        val_list = remove_duplicates_in_list(test_list)
        self.assertListEqual(val_list, [5, 6, 7, 1, 2, 3])

if __name__ == "__main__":
    unittest.main()
