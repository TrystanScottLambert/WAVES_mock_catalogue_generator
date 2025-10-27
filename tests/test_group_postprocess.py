"""
Unit tests for the group post-processing module.
"""

import unittest
import astropy.units as u
import numpy as np

from group_post_process import calc_rvir_from_mvir


class TestMvir(unittest.TestCase):
    """
    Testing the calc_rvir_from_mvir function
    """

    def test_floats(self):
        """
        Testing that function works for floats.
        """
        mvir = 1e12 * u.solMass
        redshift = 0.3
        result = calc_rvir_from_mvir(mvir, redshift)
        answer = 0.18987214
        self.assertAlmostEqual(result, answer, places=4)

    def test_array(self):
        """
        Testing that the function works on arrays of quantities
        """
        mvirs = np.array([1e12, 1e15]) * u.solMass
        redshifts = np.array([0.3, 0.2])
        results = calc_rvir_from_mvir(mvirs, redshifts)
        answers = [0.18987214, 1.97123138]
        for res, ans in zip(results, answers):
            self.assertAlmostEqual(res, ans, places=4)


if __name__ == "__main__":
    unittest.TestCase()
