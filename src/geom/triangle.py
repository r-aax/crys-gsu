"""
Triangle realization.
"""

import math
import numpy as np
from vect import Vect
from segment import Segment
# from geom.vect import Vect
# from geom.segment import Segment

# ==================================================================================================


class Triangle:
    """
    Triangle realization.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, a, b, c):
        """
        Constructor.
        :param a: First point.
        :param b: Second point.
        :param c: Third point.
        """

        self.Points = [a, b, c]
        self.BackRef = None

    # ----------------------------------------------------------------------------------------------

    def __str__(self):
        """
        String representation.
        :return: String.
        """

        return 'Triangle {0}, {1}, {2}'.format(self.Points[0], self.Points[1], self.Points[2])

    # ----------------------------------------------------------------------------------------------

    def a(self):
        """
        Get first point.
        :return: First point.
        """

        return self.Points[0]

    # ----------------------------------------------------------------------------------------------

    def b(self):
        """
        Get second point.
        :return: Second point.
        """

        return self.Points[1]

    # ----------------------------------------------------------------------------------------------

    def c(self):
        """
        Get third point.
        :return: Third point.
        """

        return self.Points[2]

    # ----------------------------------------------------------------------------------------------

    def area(self):
        """
        Area.
        :return: Area.
        """

        return 0.5 * Vect.cross_product(self.b() - self.a(), self.c() - self.b()).mod()

    # ----------------------------------------------------------------------------------------------

    def normal_orth(self):
        """
        Get normal with length 1.0.
        :return: Normal orth.
        """

        n = Vect.cross_product(self.b() - self.a(), self.c() - self.b())

        return n.orth()

    # ----------------------------------------------------------------------------------------------

    def centroid(self):
        """
        Get centroid.
        :return: Centroid.
        """

        return (self.a() + self.b() + self.c()) / 3.0

    # ----------------------------------------------------------------------------------------------

    def intersection_with_segment(self, s):
        """
        Calculate intersection with segment.
        :param s: Segment.
        :return:  Intersection points.
        """

        # Get all points.
        a, b, c = self.a(), self.b(), self.c()
        p, q = s.a(), s.b()

        #
        # Point (x, y, z) inside triangle can be represented as
        # x = x_a + (x_b - x_a) * alfa + (x_c - x_a) * beta
        # y = y_a + (y_b - y_a) * alfa + (y_c - y_a) * beta
        # z = z_a + (z_b - z_a) * alfa + (z_c - z_a) * beta
        #    where (x_a, y_a, z_a) - coordinates of point a,
        #          (x_b, y_b, z_b) - coordinates of point b,
        #          (x_c, y_c, z_c) - coordinates of point c,
        #          alfa >= 0,
        #          beta >= 0,
        #          alfa + beta <= 1.
        # ...
        # x = x_a + x_ba * alfa + x_ca * beta
        # y = y_a + y_ba * alfa + y_ca * beta
        # z = z_a + z_ba * alfa + z_ca * beta
        #
        # Point (x, y, z) on segment can be represented as
        # x = x_p + (x_q - x_p) * phi
        # y = y_p + (y_q - y_p) * phi
        # x = z_p + (z_q - z_p) * phi
        #   where (x_p, y_p, z_p) - coordinates of point p,
        #         (x_q, y_q, z_q) - coordinates of point q,
        #         0 <= phi <= 1.
        # ...
        # x = x_p + x_qp * phi
        # y = y_p + y_qp * phi
        # x = z_p + z_qp * phi
        #
        # So to find intersection we have to solve system
        # x_a + x_ba * alfa + x_ca * beta = x_p + x_qp * phi
        # y_a + y_ba * alfa + y_ca * beta = y_p + y_qp * phi
        # z_a + z_ba * alfa + z_ca * beta = z_p + z_qp * phi
        # ...
        # x_ba * alfa + x_ca * beta + (-x_qp) * phi = x_p - x_a
        # y_ba * alfa + y_ca * beta + (-y_qp) * phi = y_p - y_a
        # z_ba * alfa + z_ca * beta + (-z_qp) * phi = z_p - z_a
        # ...
        # x_ba * alfa + x_ca * beta + x_pq * phi = x_pa
        # y_ba * alfa + y_ca * beta + y_pq * phi = y_pa
        # z_ba * alfa + z_ca * beta + z_pq * phi = z_pa
        #
        # Matrix view of this system can be written in the following view
        # [x_ba x_ca x_pq]     [alfa]     [x_pa]
        # [y_ba y_ca y_pq]  X  [beta]  =  [y_pa]
        # [z_ba z_ca z_pq]     [phi ]     [z_pa]
        #

        # Vectors differences.
        ba, ca, pq, pa = b - a, c - a, p - q, p - a

        m = np.array([ba.coords_list(), ca.coords_list(), pq.coords_list()])
        m = np.transpose(m)
        d = np.linalg.det(m)

        if abs(d) < 1e-10:
            # TODO:
            # To analyze case of many intersection points.
            return []

        im = np.linalg.inv(m)
        [alfa, beta, phi] = im.dot(np.array(pa.coords_list()))

        # Check solution.
        is_inside_tri = (alfa >= 0.0) and (beta >= 0.0) and (alfa + beta <= 1.0)
        is_inside_seg = (phi >= 0.0) and (phi <= 1.0)

        if is_inside_tri and is_inside_seg:
            return [p - pq * phi]
        else:
            return []

    # ----------------------------------------------------------------------------------------------

    def point_in_triangle(self, point):
        """
        Point in triangle or not.
        :param point: Vect object
        :return: [Vect object]
        """

        # Get all points.
        a, b, c = self.a(), self.b(), self.c()

        # А point in a triangle if the sum of the formed areas is equal to the area of the original triangle.
        full_area = round(self.area(), 10)
        ar = round(Triangle(a, b, point).area() + Triangle(a, point, c).area() + Triangle(point, b, c).area(), 10)
        if full_area == ar:
            return [point]
        return []

    # ----------------------------------------------------------------------------------------------

    def in_one_plane(self, tri):
        """
        If all points of triangle 2 with triangle 1 are in the same plane, then the triangles lie in the same plane.
        :param tri: triangle 2
        :return: True or False
        """

        # Get all points.
        a, b, c = self.a(), self.b(), self.c()
        ta, tb, tc = tri.a(), tri.b(), tri.c()

        # если треугольники в одной плоскости
        A = np.array([[ta.X - a.X, ta.Y - a.Y, ta.Z - a.Z],
                      [b.X - a.X, b.Y - a.Y, b.Z - a.Z],
                      [c.X - a.X, c.Y - a.Y, c.Z - a.Z]])
        B = np.array([[tb.X - a.X, tb.Y - a.Y, tb.Z - a.Z],
                      [b.X - a.X, b.Y - a.Y, b.Z - a.Z],
                      [c.X - a.X, c.Y - a.Y, c.Z - a.Z]])
        C = np.array([[tc.X - a.X, tc.Y - a.Y, tc.Z - a.Z],
                      [b.X - a.X, b.Y - a.Y, b.Z - a.Z],
                      [c.X - a.X, c.Y - a.Y, c.Z - a.Z]])

        A = np.linalg.det(A)
        B = np.linalg.det(B)
        C = np.linalg.det(C)

        return A + B + C == 0

    # ----------------------------------------------------------------------------------------------

    def intersection_with_triangle(self, tri):
        """
        Triangles intersect if one triangle line intersects another triangle.
        :param tri: Triangle.
        :return: Triangle intersection.
        """

        # если это один и тот же объект
        if self == tri:
            return []

        # Get all points.
        a, b, c = self.a(), self.b(), self.c()
        ta, tb, tc = tri.a(), tri.b(), tri.c()

        # если треугольники в одной плоскости
        if self.in_one_plane(tri):

            pit = self.point_in_triangle
            pit2 = tri.point_in_triangle

            # треугольники совпадают
            if (a == ta) and (b == tb) and (c == tc):
                return [a, b, c]

            # прочие случаи
            res = pit(ta) + pit(tb) + pit(tc) + pit2(a) + pit2(b) + pit2(c)
            if len(res) == 2:
                if res[0] == res[1]:
                    return [res[0]]
            elif len(res) == 4:
                return res[:2]
            return res

        # если треугольники в разных плоскостях (могут быть параллельными)
        int11 = self.intersection_with_segment(Segment(ta, tb))
        int12 = self.intersection_with_segment(Segment(tb, tc))
        int13 = self.intersection_with_segment(Segment(tc, ta))
        if int11 or int12 or int13:
            res = int11 + int12 + int13
            if len(res) == 2:
                if res[0] == res[1]:
                    return [res[0]]
            return res
        else:
            int21 = tri.intersection_with_segment(Segment(a, b))
            int22 = tri.intersection_with_segment(Segment(b, c))
            int23 = tri.intersection_with_segment(Segment(c, a))
            if int21 or int22 or int23:
                res = int21 + int22 + int23
                if len(res) == 2:
                    if res[0] == res[1]:
                        return [res[0]]
                return res
            else:
                return []

# ==================================================================================================


if __name__ == '__main__':

    # тест 1 - треугольники в одной плоскости, разные объекты
    print()
    print('тест 1 - треугольники в одной плоскости, разные объекты')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 2 - треугольники в одной плоскости, один объект
    print()
    print('тест 2 - треугольники в одной плоскости, один объект')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri1)
    if res:
        print(res)

    # тест 3 - треугольники в одной плоскости, один вложен в другой
    print()
    print('тест 3 - треугольники в одной плоскости, один вложен в другой')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, 9, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 4 - треугольники в одной плоскости, одина вержина
    print()
    print('тест 4 - треугольники в одной плоскости, одина вержина')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 5, 0), Vect(10, 5, 0), Vect(5, 15, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 5 - треугольники в одной плоскости, две вержины
    print()
    print('тест 5 - треугольники в одной плоскости, две вержины')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, -10, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 6 - треугольники в одной плоскости, не пересекаются
    print()
    print('тест 6 - треугольники в одной плоскости, не пересекаются')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(15, 0, 0), Vect(20, 0, 0), Vect(15, 10, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 7 - треугольники в разных плоскостях, один объект
    print()
    print('тест 7 - треугольники в разных плоскостях, один объект')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri1)
    if res:
        print(res)

    # тест 8 - треугольники в разных плоскостях, одной вершиной
    print()
    print('тест 8 - треугольники в разных плоскостях, одной вершиной')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(5, 5, 0), Vect(10, 0, 1), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 9 - треугольники в разных плоскостях, двумя вершинами - одним ребром
    print()
    print('тест 9 - треугольники в разных плоскостях, двумя вершинами - одним ребром')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 10 - треугольники в разных плоскостях, одним ребром - одна точка пересечения
    print()
    print('тест 10 - треугольники в разных плоскостях, одним ребром - одна точка пересечения')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(5, 0.1, -1), Vect(15, 0.1, -1), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 11 - треугольники в разных плоскостях, двумя ребрами - две точки пересечения
    print()
    print('тест 11 - треугольники в разных плоскостях, двумя ребрами - две точки пересечения')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, -1), Vect(9, 0.1, -1), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 12 - треугольники в разных плоскостях, одно общее ребро
    print()
    print('тест 12 - треугольники в разных плоскостях, одно общее ребро')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0, 0), Vect(9, 0, 0), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 13 - треугольники в одной плоскости, одно общее ребро
    print()
    print('тест 13 - треугольники в одной плоскости, одно общее ребро')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0, 0), Vect(9, 0, 0), Vect(5, -10, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 14 - треугольники в разных плоскостях, одно общее ребро, одной длины
    print()
    print('тест 14 - треугольники в разных плоскостях, одно общее ребро, одной длины')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 15 - треугольники в одной плоскости, одно общее ребро, одной длины
    print()
    print('тест 15 - треугольники в одной плоскости, одно общее ребро, одной длины')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, -10, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res[0], res[1])

    # тест 16 - треугольники в разных плоскостях, одна общая вершина
    print()
    print('тест 16 - треугольники в разных плоскостях, одна общая вершина')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 17 - треугольники в одной плоскости, одна общая вершина
    print()
    print('тест 17 - треугольники в одной плоскости, одна общая вершина')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 20, 0), Vect(10, 20, 0), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 18 - треугольники в разных плоскостях, пересечение ребер
    print()
    print('тест 18 - треугольники в разных плоскостях, пересечение ребер')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 0, -1))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 19 - треугольники в разных плоскостях, параллельны
    print()
    print('тест 19 - треугольники в разных плоскостях, параллельны')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

    # тест 20 - треугольники в одной плоскости, одина вержина у обоих
    print()
    print('тест 20 - треугольники в одной плоскости, одина вержина у обоих')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(5, 1, 0), Vect(0, 11, 0), Vect(10, 11, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res[0], res[1])

    # тест 21 - треугольники в одной плоскости, звезда Давида
    print()
    print('тест 21 - треугольники в одной плоскости, звезда Давида')
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 7, 0), Vect(10, 7, 0), Vect(5, -3, 0))
    res = tri1.intersection_with_triangle(tri2)
    if res:
        print(res)

# ==================================================================================================
