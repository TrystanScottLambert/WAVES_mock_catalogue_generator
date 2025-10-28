"""
Unit tests for the group post-processing module.
"""

import unittest
import astropy.units as u
import numpy as np


from group_post_process import (
    calc_rvir_from_mvir,
    skycoord_to_cartesian_vectorized,
    join_groups,
)


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


class TestGrouping(unittest.TestCase):
    """
    Testing the join_groups function.
    """

    def test_simple(self):
        """
        Simple test case of two groups far away and one isoalted group
        """
        masses = np.repeat(1e15, 7) * u.solMass

        ra = np.array([0, 0, 0, 180, 180, 180, 180])
        dec = np.array([0, 0, 0, 0, 0, 0, -45])
        z = np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1])
        ids = ["a", "b", "c", "d", "e", "f", "g"]
        centers = skycoord_to_cartesian_vectorized(ra, dec, z)
        rvirs = calc_rvir_from_mvir(masses, z)
        mapping = join_groups(centers, rvirs, ids)
        self.assertEqual(mapping["a"], 1)
        self.assertEqual(mapping["b"], 1)
        self.assertEqual(mapping["c"], 1)
        self.assertEqual(mapping["d"], 2)
        self.assertEqual(mapping["e"], 2)
        self.assertEqual(mapping["f"], 2)
        self.assertEqual(mapping["g"], -1)


if __name__ == "__main__":
    unittest.TestCase()
