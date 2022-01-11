"""
Triangle realization.
"""

import math
import numpy as np
from sympy import Float
import itertools
from geom.vect import Vect
from geom.segment import Segment

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
        ba, ca, pq, pa, qp = b - a, c - a, p - q, p - a, q - p

        m = np.array([ba.coords_list(), ca.coords_list(), pq.coords_list()])
        m = np.transpose(m)
        d = np.linalg.det(m)

        def determinant(i, j):
            """

            :param i: Axis i.
            :param j: Axis j.
            :return: Determinant.
            """

            det = np.array([[ba.Coords[i], ca.Coords[i]],
                            [ba.Coords[j], ca.Coords[j]]])
            return round(np.linalg.det(det), 10)

        def solve_the_linear_equation(a, b):
            """
            Solution of the equation 0 <= a*x + b.

            :param a: Param a.
            :param b: Param b.
            :return: Range of solutions for x.
            """

            if a > 0:
                return [-b / a, float('Inf')]
            elif a < 0:
                return [float('-Inf'), -b / a]
            else:
                if b >= 0:
                    return [float('-Inf'), float('Inf')]
                else:
                    return []

        def intersection_of_intervals(data, aps):
            """

            :param aps: Rounding accuracy.
            :param data: List of intervals.
            :return: A list of one interval.
            """

            res = []
            if len(data) == 3:
                res = [np.max(data, axis=0)[0], np.min(data, axis=0)[1]]
                if res[0] + aps > res[1] or res[0] > res[1] + aps or \
                        (res[0] - aps > 1 and res[1] > 1) or (res[1] + aps < 0 and res[0] < 0):
                    return []
                elif math.isclose(res[0], 1) and res[1] > 1:
                    return [res[0]]
                elif math.isclose(res[1], 0) and res[0] < 0:
                    return [res[1]]
                elif res[0] - aps <= 0 and res[1] + aps >= 1:
                    res = [0, 1]
                elif res[0] - aps <= 0 and res[1] < 1:
                    res = [0, res[1]]
                elif res[0] > 0 and res[1] + aps >= 1:
                    res = [res[0], 1]

                try1 = data[0][0] - aps <= res[0] and data[0][1] + aps >= res[1]
                try2 = data[1][0] - aps <= res[0] and data[1][1] + aps >= res[1]
                try3 = data[2][0] - aps <= res[0] and data[2][1] + aps >= res[1]
                if not (try1 and try2 and try3):
                    res = []
            return res

        if abs(d) < 1e-10:

            # one plane
            if self.point_in_one_plane(p) and self.point_in_one_plane(q):

                """
                a - alf*ba + bet*ca = p + fi*qp; alf >= 0, bet >= 0, alf + bet <= 1, 0 <= fi <= 1
    
                alf*ba + bet*ca + fi*pq = pa
                alf*ba + bet*ca = fi*qp + pa
    
                alf*ba.x + bet*ca.x = fi*qp.x + pa.x
                alf*ba.y + bet*ca.y = fi*qp.y + pa.y
    
                det = [[ba.x, ca.x],
                       [ba.y, ca.y]]
    
                det_alf = [[fi*qp.x + pa.x, ca.x],
                           [fi*qp.y + pa.y, ca.y]]
                det_bet = [[ba.x, fi*qp.x + pa.x],
                           [ba.y, fi*qp.y + pa.y]]
    
                alf = det_alf / det
                bet = det_bet / det
    
                alf >= 0; bet >= 0; alf + bet <= 1; 0 <= fi <= 1
                """

                for i in range(2):
                    for j in range(1, 3):
                        if i != j:

                            det = determinant(i, j)
                            if det == 0:
                                continue

                            res = []

                            # alf
                            a_alf = np.array([[qp[i], ca[i]],
                                              [qp[j], ca[j]]])
                            a_alf = np.linalg.det(a_alf) / det
                            b_alf = np.array([[pa[i], ca[i]],
                                              [pa[j], ca[j]]])
                            b_alf = np.linalg.det(b_alf) / det
                            res_a = solve_the_linear_equation(a_alf, b_alf)
                            if res_a:
                                res = res + [res_a]

                            # bet
                            a_bet = np.array([[ba[i], qp[i]],
                                              [ba[j], qp[j]]])
                            a_bet = np.linalg.det(a_bet) / det
                            b_bet = np.array([[ba[i], pa[i]],
                                              [ba[j], pa[j]]])
                            b_bet = np.linalg.det(b_bet) / det
                            res_b = solve_the_linear_equation(a_bet, b_bet)
                            if res_b:
                                res = res + [res_b]

                            # alf + bet
                            a_alf_bet = -a_alf - a_bet
                            b_alf_bet = -b_alf - b_bet + 1
                            res_ab = solve_the_linear_equation(a_alf_bet, b_alf_bet)
                            if res_ab:
                                res = res + [res_ab]

                            fi = intersection_of_intervals(res, 1.0e-10) # TODO it is magic constant
                            if not fi:
                                return []
                            else:
                                return [p + qp * i for i in fi]
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
        The method answers the question in which part of the triangle is this point.

        Parameters
        ----------
        point - Vect object

        Returns
        -------
        0 (Int) - the point does not belong to the triangle
        1 (Int) - the point belongs to the vertex of the triangle
        2 (Int) - the point belongs to the edge of the triangle
        3 (Int) - the point belongs to the inner space of the triangle

        """

        # Get all points.
        a, b, c = self.a(), self.b(), self.c()

        # А point in a triangle if the sum of the formed areas is equal to the area of the original triangle.
        ar1 = Triangle(a, b, point).area()
        ar2 = Triangle(a, point, c).area()
        ar3 = Triangle(point, b, c).area()
        if math.isclose(self.area(), ar1 + ar2 + ar3):
            return bool(ar1) + bool(ar2) + bool(ar3)
        return 0

    # ----------------------------------------------------------------------------------------------

    def point_in_one_plane(self, p):
        """
        a point in one plane if::

                  [p.X - a.X, p.Y - a.Y, p.Z - a.Z]
          det =   [b.X - a.X, b.Y - a.Y, b.Z - a.Z] = 0
                  [c.X - a.X, c.Y - a.Y, c.Z - a.Z]

        :param p: Vect
        :return: True or False
        """

        # Get all points.
        a, b, c = self.a(), self.b(), self.c()

        det = np.array([[p.X - a.X, p.Y - a.Y, p.Z - a.Z],
                        [b.X - a.X, b.Y - a.Y, b.Z - a.Z],
                        [c.X - a.X, c.Y - a.Y, c.Z - a.Z]])

        det = np.linalg.det(det)

        if det == 0:
            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------------

    def in_one_plane(self, tri):
        """
        If all points of triangle 2 with triangle 1 are in the same plane, then the triangles lie in the same plane.

        :param tri: triangle 2
        :return: True or False
        """

        # Get all points.
        ta, tb, tc = tri.a(), tri.b(), tri.c()

        ta = self.point_in_one_plane(ta)
        tb = self.point_in_one_plane(tb)
        tc = self.point_in_one_plane(tc)
        return ta and tb and tc

    # ----------------------------------------------------------------------------------------------

    def intersection_with_triangle(self, tri):
        """
        Triangles intersect if one triangle line intersects another triangle.

        :param tri: Triangle.
        :return: Triangle intersection.
        """

        # if it is the same object
        if id(self) == id(tri):
            return []

        # Get all points.
        a, b, c = self.a(), self.b(), self.c()
        ta, tb, tc = tri.a(), tri.b(), tri.c()

        int11 = self.intersection_with_segment(Segment(ta, tb))
        int12 = self.intersection_with_segment(Segment(tb, tc))
        int13 = self.intersection_with_segment(Segment(tc, ta))
        int21 = tri.intersection_with_segment(Segment(a, b))
        int22 = tri.intersection_with_segment(Segment(b, c))
        int23 = tri.intersection_with_segment(Segment(c, a))
        res = int11 + int12 + int13 + int21 + int22 + int23
        if res:
            res.sort()
            return [res[0]] + [res[i] for i in range(1, len(res)) if res[i] != res[i - 1]]
        return res

    # ----------------------------------------------------------------------------------------------

    def all_corners_smaller(self, corner):
        """

        Parameters
        ----------
        corner - an angle whose value is greater than is not allowed. takes values from 1 to 180 (Float)

        Returns
        -------
        triangle is good or no

        """

        a, b, c = self.a(), self.b(), self.c()
        ba = b - a
        ca = c - a
        corn1 = Vect.get_angle_in_deg(ba, ca) < corner
        ab = a - b
        cb = c - b
        corn2 = Vect.get_angle_in_deg(ab, cb) < corner
        ac = a - c
        bc = b - c
        corn3 = Vect.get_angle_in_deg(ac, bc) < corner
        return corn1 and corn2 and corn3

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def sorting_by_the_selected_axis(triangles_for_sorting, axis):
        """

        Parameters
        ----------
        triangles_for_sorting - list of triangles to sort
        axis - sorting axis

        Returns
        -------
        an array of triangles sorted by the selected axis

        """

        triangles_for_sorting.sort(key=lambda tri: tri.centroid()[axis])
        return triangles_for_sorting

    # ----------------------------------------------------------------------------------------------

    def __eq__(self, other):
        """
        Compares coordinate triangles

        Parameters
        ----------
        other - Triangle object

        Returns
        -------
        Result (True or False)

        """

        if self.centroid() == other.centroid():
            if self.a() == other.a() or self.a() == other.b() or self.a() == other.c():
                if self.b() == other.a() or self.b() == other.b() or self.b() == other.c():
                    return True
        return False

    # ----------------------------------------------------------------------------------------------

    def points_for_split_the_triangle(self):
        """

        Divides the triangle into two on the longest side.

        Returns
        -------
        A list of two points ['the vertex opposite the long side', 'center of the long side'] (Vect object).

        """

        a, b, c = self.a(), self.b(), self.c()
        ab = a.dist_to(b)
        bc = b.dist_to(c)
        ca = c.dist_to(a)
        if ab > bc:
            if ab > ca:
                return [c, (a + b) / 2]
            else:
                return [b, (a + c) / 2]
        else:
            if bc > ca:
                return [a, (c + b) / 2]
            else:
                return [b, (a + c) / 2]

    # ----------------------------------------------------------------------------------------------

    def split_the_triangle_on_the_split_point(self, position, point):
        """
        Splits the triangle into new triangles depending on the split point.

        Parameters
        ----------
        point - the split point
        position - position of the point in the triangle

        Returns
        -------
        List of new triangles.

        """

        if position == 0:
            return self
        elif position == 1:
            return self
        elif position == 2:
            if Vect.point_on_vector(self.a(), self.b(), point):
                return [Triangle(self.a(), point, self.c()), Triangle(point, self.b(), self.c())]
            elif Vect.point_on_vector(self.b(), self.c(), point):
                return [Triangle(self.a(), self.b(), point), Triangle(self.a(), point, self.c())]
            elif Vect.point_on_vector(self.c(), self.a(), point):
                return [Triangle(self.a(), self.b(), point), Triangle(point, self.b(), self.c())]
        elif position == 3:
            # TODO дописать разбиение треугольника, когда точка внутри него
            pass
        return []

    # ----------------------------------------------------------------------------------------------

    def rearranging_intersecting_triangles(self, other):
        """

        Parameters
        ----------
        other - another triangle

        Returns
        -------
        list of objects triangle or []

        """

        res = []

        # получаем список точек пересечений обоих треугольников
        list_of_points = self.intersection_with_triangle(other)

        # перестраиваем каждый их двух треугольников
        for tri in [self, other]:
            list_of_new_triangles = [tri]
            points_in_triangel = []

            # находим точки пересечений для этого треугольника
            for point in list_of_points:
                answer = tri.point_in_triangle(point)
                if answer != 0 and answer != 1:
                    points_in_triangel.append([answer, point])

            if points_in_triangel:
                # сортируем, чтобы точки в ребрах были первыми
                points_in_triangel.sort()

                # для каждой точки проводим разбиение треугольника
                # TODO тип расположения точки в треугольнике (0, 1, 2 или 3) указаны для исходного треугольника,
                #  а для генерируемых нужно относительно них пересчитывать - сейчас не так
                #  скорее всего придется генератор заменить на цикл с опознанием расположения точек
                #  для каждого нового треугольника

                for p in points_in_triangel:
                    list_of_new_triangles = list(itertools.chain(*[t.split_the_triangle_on_the_split_point(*p) \
                                                                   for t in list_of_new_triangles]))

                # записываем список треугольников в общий список для вывода
                res += list_of_new_triangles

        if res:
            # проверяем на наличие одинаковых треугольников и удаляем их из списка
            # в классе TrianglesCloud написать статический, который принимает список треугольников,
            # проверяет список на наличие одинаковых треугольников, удаляет одинаковые,
            # возвращает очищенный от повторений список треугольников
            return TrianglesCloud.removing_repetitions(res)

        return []

# ==================================================================================================


if __name__ == '__main__':

    # from vect import Vect
    # from segment import Segment

    # rearranging_intersecting_triangles
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, 9, 0))
    res = tri1.rearranging_intersecting_triangles(tri2)
    assert(res == [])

    # point_in_triangle
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    p = Vect(-5, 0, 0)
    assert(tri1.point_in_triangle(p) == 0)

    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    p = Vect(0, 0, 0)
    assert(tri1.point_in_triangle(p) == 1)

    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    p = Vect(5, 0, 0)
    assert(tri1.point_in_triangle(p) == 2)

    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    p = Vect(5, 1, 0)
    assert(tri1.point_in_triangle(p) == 3)

    # is_good_triangle
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    assert(tri1.all_corners_smaller(90))

    # intersection_with_triangle
    # тест 1 - одинаковые треугольники в одной плоскости, разные объекты, пересечения по 3 вершинам
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 3)

    # тест 2 - одинаковые треугольники в одной плоскости, один объект, нет пересечений
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri1)
    assert(len(res) == 0)

    # тест 3 - треугольники в одной плоскости, один вложен в другой, пересечения по 3 вершинам
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, 9, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 3)

    # тест 4 - треугольники в одной плоскости, одина вержина, пересечения по 1 вершине и 2 точки на ребре
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 5, 0), Vect(10, 5, 0), Vect(5, 15, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 3)

    # тест 5 - треугольники в одной плоскости, две вержины, пересечения по 2 вершинам и 2 точки на ребре
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, -10, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 4)

    # тест 6 - треугольники в одной плоскости, не пересекаются, нет пересечений
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(15, 0, 0), Vect(20, 0, 0), Vect(15, 10, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 0)

    # тест 7 - треугольники в разных плоскостях, одной вершиной
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(5, 5, 0), Vect(10, 0, 1), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 1)

    # тест 8 - треугольники в разных плоскостях, двумя вершинами - одним ребром, пересечения по 2 вершинам
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, 0), Vect(9, 0.1, 0), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 2)

    # тест 9 - треугольники в разных плоскостях, пересекаются, у каждого по одной точке пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(5, 0.1, -1), Vect(15, 0.1, -1), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 2)

    # тест 10 - треугольники в разных плоскостях, один протыкает вершиной второго, две точки пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0.1, -1), Vect(9, 0.1, -1), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 2)

    # тест 11 - треугольники в разных плоскостях, одно общее ребро, две точки пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0, 0), Vect(9, 0, 0), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 2)

    # тест 12 - треугольники в одной плоскости, одно общее ребро, две точки пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(1, 0, 0), Vect(9, 0, 0), Vect(5, -10, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 2)

    # тест 13 - треугольники в разных плоскостях, одно общее ребро, одной длины, две точки пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 2)

    # тест 14 - треугольники в одной плоскости, одно общее ребро, одной длины, две точки пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, -10, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 2)

    # тест 15 - треугольники в разных плоскостях, одна общая вершина, 1 точка пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 1)

    # тест 16 - треугольники в одной плоскости, одна общая вершина, 1 точка пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 20, 0), Vect(10, 20, 0), Vect(5, 10, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 1)

    # тест 17 - треугольники в разных плоскостях, пересечение ребер, две точки пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 0, -1))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 2)

    # тест 18 - треугольники в разных плоскостях, параллельны, нет пересечений
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 0, 1), Vect(10, 0, 1), Vect(5, 10, 1))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 0)

    # тест 19 - треугольники в одной плоскости, одина вершина у обоих внутри другого, 4 точки пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(5, 1, 0), Vect(0, 11, 0), Vect(10, 11, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 4)

    # тест 20 - треугольники в одной плоскости, звезда Давида, 6 точек пересечения
    tri1 = Triangle(Vect(0, 0, 0), Vect(10, 0, 0), Vect(5, 10, 0))
    tri2 = Triangle(Vect(0, 7, 0), Vect(10, 7, 0), Vect(5, -3, 0))
    res = tri1.intersection_with_triangle(tri2)
    assert(len(res) == 6)

# ==================================================================================================
