import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../crysgsu/')
from geom import Vect, Triangle


class TestUtils(unittest.TestCase):

    def test_point_in_triangle(self):
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        p = Vect(-5, 0, 0)
        self.assertEqual(tri1.point_in_triangle(p), 0)

        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        p = Vect(0, 0, 0)
        self.assertEqual(tri1.point_in_triangle(p), 1)

        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        p = Vect(5, 0, 0)
        self.assertEqual(tri1.point_in_triangle(p), 2)

        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        p = Vect(5, 1, 0)
        self.assertEqual(tri1.point_in_triangle(p), 3)

    def test_is_good_triangle(self):
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        self.assertTrue(tri1.all_corners_smaller(90))

    def test_intersection_with_triangle(self):
        # тест 1 - одинаковые треугольники в одной плоскости, разные объекты, пересечения по 3 вершинам
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 3)

        # тест 2 - одинаковые треугольники в одной плоскости, один объект, нет пересечений
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        res = tri1.intersection_with_triangle(tri1)
        self.assertTrue(len(res) == 0)

        # тест 3 - треугольники в одной плоскости, один вложен в другой, пересечения по 3 вершинам
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, 9, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 3)

        # тест 4 - треугольники в одной плоскости, одина вержина, пересечения по 1 вершине и 2 точки на ребре
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 5, 0), Vect(10, 5, 0), Vect(5, 15, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 3)

        # тест 5 - треугольники в одной плоскости, две вержины, пересечения по 2 вершинам и 2 точки на ребре
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, -10, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 4)

        # тест 6 - треугольники в одной плоскости, не пересекаются, нет пересечений
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(15, 0, 0), Vect(20, 0, 0), Vect(15, 10, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 0)

        # тест 7 - треугольники в разных плоскостях, одной вершиной
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(5, 5, 0), Vect(10, 0, 1), Vect(5, 10, 1))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 1)

        # тест 8 - треугольники в разных плоскостях, двумя вершинами - одним ребром, пересечения по 2 вершинам
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, 10, 1))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 2)

        # тест 9 - треугольники в разных плоскостях, пересекаются, у каждого по одной точке пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(5, 0.1, -1), Vect(15, 0.1, -1), Vect(5, 10, 1))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 2)

        # тест 10 - треугольники в разных плоскостях, один протыкает вершиной второго, две точки пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(1, 0.1, -1), Vect(9, 0.1, -1), Vect(5, 10, 1))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 2)

        # тест 11 - треугольники в разных плоскостях, одно общее ребро, две точки пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(1, 0, 0), Vect(9, 0, 0), Vect(5, 10, 1))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 2)

        # тест 12 - треугольники в одной плоскости, одно общее ребро, две точки пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(1, 0, 0), Vect(9, 0, 0), Vect(5, -10, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 2)

        # тест 13 - треугольники в разных плоскостях, одно общее ребро, одной длины, две точки пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 1))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 2)

        # тест 14 - треугольники в одной плоскости, одно общее ребро, одной длины, две точки пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, -10, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 2)

        # тест 15 - треугольники в разных плоскостях, одна общая вершина, 1 точка пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 10, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 1)

        # тест 16 - треугольники в одной плоскости, одна общая вершина, 1 точка пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 20, 0), Vect(10, 20, 0), Vect(5, 10, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 1)

        # тест 17 - треугольники в разных плоскостях, пересечение ребер, две точки пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 0, -1))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 2)

        # тест 18 - треугольники в разных плоскостях, параллельны, нет пересечений
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 10, 1))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 0)

        # тест 19 - треугольники в одной плоскости, одина вершина у обоих внутри другого, 4 точки пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(5, 1, 0), Vect(0, 11, 0), Vect(10, 11, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 4)

        # тест 20 - треугольники в одной плоскости, звезда Давида, 6 точек пересечения
        tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
        tri2 = Triangle(Vect(0, 7, 0), Vect(10, 7, 0), Vect(5, -3, 0))
        res = tri1.intersection_with_triangle(tri2)
        self.assertTrue(len(res) == 6)


t = TestUtils()
t.test_point_in_triangle()
t.test_is_good_triangle()
t.test_intersection_with_triangle()
