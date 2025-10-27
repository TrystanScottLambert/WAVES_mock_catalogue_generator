"""
Unit tests for the group post-processing module.
"""

import unittest
import astropy.units as u
import numpy as np

from group_post_process import calc_rvir_from_mvir, Group, join_groups


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


class TestOverlap(unittest.TestCase):
    """
    Testing that the overlap method in a group works.
    """

    group_1 = Group(1e15, 0.0, 0.0, 0.2, "a")
    group_2 = Group(1e15, 0.0, 0.0, 0.2, "b")
    group_3 = Group(1e15, 180, 0, 0.2, "c")
    group_4 = Group(1e15, 0.0, 0.0, 1, "d")

    def test_trivial_case(self):
        """test true for same group"""
        self.assertTrue(self.group_1.overlap(self.group_2))
        self.assertTrue(self.group_2.overlap(self.group_1))

    def test_on_sky_sep(self):
        """
        Testing fails when we separate group on sky.
        """
        self.assertFalse(self.group_1.overlap(self.group_3))

    def test_behind(self):
        """
        Testing fail when groups are behind one another.
        """
        self.assertFalse(self.group_1.overlap(self.group_4))


class TestGrouping(unittest.TestCase):
    """
    Testing the join_groups function.
    """

    def test_simple(self):
        """
        Simple test case of two groups far away and one isoalted group
        """
        group_1 = Group(1e15, 0, 0, 0.2, "a")
        group_2 = Group(1e15, 0, 0, 0.2, "b")
        group_3 = Group(1e15, 0, 0, 0.2, "c")

        group_4 = Group(1e15, 180, 0, 0.2, "d")
        group_5 = Group(1e15, 180, 0, 0.2, "e")
        group_6 = Group(1e15, 180, 0, 0.2, "f")

        group_7 = Group(1e15, 180, -45, 1, "g")  # isolated

        groups = [group_1, group_4, group_5, group_7, group_3, group_2, group_6]
        mapping = join_groups(groups)
        self.assertEqual(mapping["a"], 1)
        self.assertEqual(mapping['b'], 1)
        self.assertEqual(mapping['c'], 1)
        self.assertEqual(mapping['d'], 2)
        self.assertEqual(mapping['e'], 2)
        self.assertEqual(mapping['f'], 2)
        self.assertEqual(mapping['g'], -1)


if __name__ == "__main__":
    unittest.TestCase()
