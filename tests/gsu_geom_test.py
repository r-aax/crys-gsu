import unittest
from src.geom.vect import Vect
from src.gsu_geom import Triangle


class TestUtils(unittest.TestCase):

    def test_no_intersection_by_boxes(self):
        t1 = Triangle(Vect(0.0, 0.0, 0.0), Vect(1.0, 0.0, 0.0), Vect(0.0, 1.0, 0.0))
        t2 = Triangle(Vect(6.0, 0.0, 0.0), Vect(5.0, 1.0, 0.0), Vect(5.0, 0.0, 0.0))
        self.assertTrue(t1.is_no_intersection_with_triangle_by_boxes(t2))

    def test_area(self):
        t = Triangle(Vect(0.0, 0.0, 0.0), Vect(1.0, 0.0, 0.0), Vect(0.0, 1.0, 0.0))
        self.assertTrue(abs(t.area() - 0.5) < 1.0e-10)


t = TestUtils()
t.test_area()
t.test_no_intersection_by_boxes()
