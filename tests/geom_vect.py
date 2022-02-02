import unittest
from src.geom.vect import Vect

class TestVectArithmetics(unittest.TestCase):

    def setUp(self):
        self.zero = Vect(0, 0, 0)
        self.one = Vect(1, 1, 1)
        self.x_two = Vect(2, 0, 0)

    def test_add(self):
        self.assertEqual(self.zero + self.one, self.one)
        self.assertEqual(self.one + self.x_two, Vect(3, 1, 1))

    def test_sub(self):
        self.assertEqual(self.zero - self.x_two, Vect(-2, 0, 0))
        self.assertEqual(self.one - self.zero, self.one)

    def test_mul(self):
        self.assertEqual(self.zero * 5, self.zero)
        self.assertEqual(self.one * 5, Vect(5, 5, 5))
        self.assertEqual(self.one * self.x_two, 2)
