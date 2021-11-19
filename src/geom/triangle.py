"""
Triangle realization.
"""

import math
import numpy as np
from geom.vect import Vect

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
        # ...
        # x = x_a + x_ba * alfa + x_ca * beta
        # y = y_a + y_ba * alfa + y_ca * beta
        # z = z_a + z_ba * alfa + z_ca * beta
        #
        # Point (x, y, z) on segment can be represented as
        # x = x_p + (x_q - x_p) * phi
        # y = y_p + (y_q - y_p) * phi
        # x = z_p + (z_q - z_p) * phi
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
        # Matrix of this system can be written in the following view
        # [x_ba x_ca x_pq
        #  y_ba y_ca y_pq
        #  z_ba z_ca z_pq]

        m = np.array([[b.X - a.X, c.X - a.X, p.X - q.X],
                      [b.Y - a.Y, c.Y - a.Y, p.Y - q.Y],
                      [b.Z - a.Z, c.Z - a.Z, p.Z - q.Z]])
        d = np.linalg.det(m)

        if abs(d) < 1e-10:
            # TODO:
            # To analyze case of many intersection points.
            return []

        im = np.linalg.inv(m)
        r = im.dot(np.array([p.X - a.X, p.Y - a.Y, p.Z - a.Z]))
        alfa, beta, phi = r[0], r[1], r[2]

        # Check solution.
        is_inside_tri = (alfa >= 0.0) and (beta >= 0.0) and (alfa + beta <= 1.0)
        is_inside_seg = (phi >= 0.0) and (phi <= 1.0)

        if is_inside_tri and is_inside_seg:
            return [p + (q - p) * phi]
        else:
            return []

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
