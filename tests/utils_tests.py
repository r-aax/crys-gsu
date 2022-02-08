import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../crysgsu/')
import utils


class TestUtils(unittest.TestCase):

    def test_flatten(self):
        self.assertEqual(utils.flatten([]), [])
        self.assertEqual(utils.flatten([1]), [1])
        self.assertEqual(utils.flatten([1, 2]), [1, 2])
        self.assertEqual(utils.flatten([1, [2, 3]]), [1, 2, 3])
        self.assertEqual(utils.flatten([1, [2, [3], 4], 5]), [1, 2, 3, 4, 5])

    def test_has_filename_timestamp(self):
        self.assertFalse(utils.has_filename_timestamp(''))
        self.assertFalse(utils.has_filename_timestamp('short'))
        self.assertFalse(utils.has_filename_timestamp('name_000000aaa000'))
        self.assertFalse(utils.has_filename_timestamp('nameX000111222333'))
        self.assertTrue(utils.has_filename_timestamp('name_000000000000'))
        self.assertTrue(utils.has_filename_timestamp('name_111222333444'))
        self.assertFalse(utils.has_filename_timestamp('name_111222333444_r'))

    def test_get_filename_and_timestamp_pair(self):
        self.assertEqual(utils.get_filename_and_timestamp_pair('name'),
                         ('name', '000000000000'))
        self.assertEqual(utils.get_filename_and_timestamp_pair('name_000111222333'),
                         ('name', '000111222333'))

    def test_dummy_split_filename(self):
        self.assertEqual(utils.dummy_split_filename('bunny_00001_111222333444.cry'),
                         ('bunny', 1, 111222333444, 'cry'))

    def test_is_filename_correct_crys_crytxt_file(self):
        self.assertTrue(utils.is_filename_correct_crys_cry_file('bunny_00004_444333222111.cry', 'bunny'))
        self.assertTrue(utils.is_filename_correct_crys_txt_file('bunny_00004_444333222111.txt', 'bunny'))
