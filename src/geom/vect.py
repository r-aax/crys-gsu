"""
Vect realization.
"""

import math
import numpy as np

# ==================================================================================================


class Vect:
    """
    Vector of three coordinates.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, x=0.0, y=0.0, z=0.0):
        """
        Constructor.

        :param x: Coord X.
        :param y: Coord Y.
        :param z: Coord Z.
        """

        self.Coords = [x, y, z]

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def from_iterable(t):
        """
        Constructor from iterable object.

        :param t: Tuple.
        :return:  Vector.
        """

        return Vect(t[0], t[1], t[2])

    # ----------------------------------------------------------------------------------------------

    @property
    def X(self):
        """
        Coord X.

        :return: Coord X.
        """
        return self.Coords[0]

    # ----------------------------------------------------------------------------------------------

    @X.setter
    def X(self, value):
        """
        Set X coord.

        :param value: Value
        """
        self.Coords[0] = value

    # ----------------------------------------------------------------------------------------------

    @property
    def Y(self):
        """
        Coord Y.

        :return: Coord Y.
        """
        return self.Coords[1]

    # ----------------------------------------------------------------------------------------------

    @Y.setter
    def Y(self, value):
        """
        Set Y coord.

        :param value: Value
        """
        self.Coords[1] = value

    # ----------------------------------------------------------------------------------------------

    @property
    def Z(self):
        """
        Coord Z.

        :return: Coord Z.
        """
        return self.Coords[2]

    # ----------------------------------------------------------------------------------------------

    @Z.setter
    def Z(self, value):
        """
        Set Z coord.

        :param value: Value
        """
        self.Coords[2] = value

    # ----------------------------------------------------------------------------------------------

    def __getitem__(self, item):
        """
        Get coordinate.

        :param item: Index.
        :return:     Coordinate.
        """

        return self.Coords[item]

    # ----------------------------------------------------------------------------------------------

    def __setitem__(self, key, value):
        """
        Set coordinate.

        :param key:   Index.
        :param value: Value.
        """

        self.Coords[key] = value

    # ----------------------------------------------------------------------------------------------

    def __str__(self):
        """
        String representation.

        :return: String.
        """

        return 'Vect ({0}, {1}, {2})'.format(self.X, self.Y, self.Z)

    # ----------------------------------------------------------------------------------------------

    def coords_tuple(self):
        """
        Get coordinates as a tuple.

        :return: Coordinates as a tuple.
        """

        return self.X, self.Y, self.Z

    # ----------------------------------------------------------------------------------------------

    def coords_list(self):
        """
        Get coordinates as a list.

        :return: Coordinates as a list.
        """

        return [self.X, self.Y, self.Z]

    # ----------------------------------------------------------------------------------------------

    def rounded_coords_tuple(self, digits):
        """
        Get tuple of rounded coordinates.

        :param digits: Digits number.
        :return:       Tuple of rounded coordinates.
        """

        return round(self.X, digits), round(self.Y, digits), round(self.Z, digits)

    # ----------------------------------------------------------------------------------------------

    def round_vect(self, digits):
        """
        Get Vect of rounded Points.

        :param digits: Digits number.
        :return: Vect with rounded coordinates.
        """

        self.X = round(self.X, digits)
        self.Y = round(self.Y, digits)
        self.Z = round(self.Z, digits)

    # ----------------------------------------------------------------------------------------------

    def __add__(self, v):
        """
        Addition of two vectors.

        :param v: Vector.
        :return:  Result.
        """

        return Vect(self.X + v.X, self.Y + v.Y, self.Z + v.Z)

    # ----------------------------------------------------------------------------------------------

    def __sub__(self, v):
        """
        Subtraction of two vectors.

        :param v: Vector.
        :return:  Result.
        """

        return Vect(self.X - v.X, self.Y - v.Y, self.Z - v.Z)

    # ----------------------------------------------------------------------------------------------

    def __mul__(self, k):
        """
        Multiplication vector on number.

        :param k: Value.
        :return:  Result (vector).
        """

        return Vect(self.X * k, self.Y * k, self.Z * k)

    # ----------------------------------------------------------------------------------------------

    def __truediv__(self, k):
        """
        Division on float value.

        :param k: Value.
        :return:  Result (vector).
        """

        return self * (1.0 / k)

    # ----------------------------------------------------------------------------------------------

    def __lt__(self, other):
        """
        Compares coordinate values.

        :param other: Vector.
        :return: Result (True or False)
        """

        if abs(self.X - other.X) < 1.0e-10 and abs(self.Y - other.Y) < 1.0e-10 and self.Z < other.Z:
            return True
        elif abs(self.X - other.X) < 1.0e-10 and self.Y < other.Y:
            return True
        elif self.X < other.X:
            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------------

    def __eq__(self, other):
        """
        Compares coordinate values.

        :param other: Vector.
        :return: Result (True or False)
        """

        return self.is_near(other, 1.0e-10)

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def dot_product(a, b):
        """
        Dot product of two vectors.

        :param a: First vector.
        :param b: Second vector.
        :return:  Dot product.
        """

        return a.X * b.X + a.Y * b.Y + a.Z * b.Z

    # ----------------------------------------------------------------------------------------------

    def mod2(self):
        """
        Square of mod.

        :return: Square of mod.
        """

        return Vect.dot_product(self, self)

    # ----------------------------------------------------------------------------------------------

    def mod(self):
        """
        Mod.

        :return: Mod.
        """

        return math.sqrt(self.mod2())

    # ----------------------------------------------------------------------------------------------

    def dist_to(self, v):
        """
        Distance to another vector.

        :param v: Vector.
        :return:  Distance.
        """

        return (self - v).mod()

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def point_on_vector(v1, v2, p):
        """
        is the point on the vector or not

        Parameters
        ----------
        v1 - the starting point of the vector
        v2 - the end point of the vector
        p - point

        Returns
        -------
        True or False

        """

        return math.isclose(v1.dist_to(v2), v1.dist_to(p) + p.dist_to(v2))


    # ----------------------------------------------------------------------------------------------

    def is_near(self, v, eps):
        """
        Check if vector near another vector.

        :param v:   Vector.
        :param eps: Epsilon.
        :return:    True - if near, False - otherwise.
        """

        return self.dist_to(v) < eps

    # ----------------------------------------------------------------------------------------------

    def orth(self):
        """
        Get orth of the vector.

        :return: Orth.
        """

        return self / self.mod()

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def angle_cos(a, b):
        """
        Cosine of angle between vectors.

        :param a: First vector.
        :param b: Second vector.
        :return:  Angle cosine.
        """

        return Vect.dot_product(a, b) / (a.mod() * b.mod())

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def get_angle_in_rad(v1, v2):
        """

        Parameters
        ----------
        v1 - vector 1
        v2 - vector 2

        Returns
        -------
        the angle between vectors 1 and 2 in rad.

        """

        return np.arccos(Vect.angle_cos(v1, v2))

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def get_angle_in_deg(v1, v2):
        """

        Parameters
        ----------
        v1 - vector 1
        v2 - vector 2

        Returns
        -------
        the angle between vectors 1 and 2 in deg.

        """

        return np.rad2deg(Vect.get_angle_in_rad(v1, v2))

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def cross_product(a, b):
        """
        Cross product of two vectors.

        :param a: First vector.
        :param b: Second vector.
        :return:  Cross product.
        """

        return Vect(a.Y * b.Z - a.Z * b.Y,
                    a.Z * b.X - a.X * b.Z,
                    a.X * b.Y - a.Y * b.X)

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
