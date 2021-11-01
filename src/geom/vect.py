"""
Vect realization.
"""

import math

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

        self.X = x
        self.Y = y
        self.Z = z

    # ----------------------------------------------------------------------------------------------

    def __str__(self):
        """
        String representation.
        :return: String.
        """

        return 'Vect ({0}, {1}, {2})'.format(self.X, self.Y, self.Z)

    # ----------------------------------------------------------------------------------------------

    def __add__(self, v):
        """
        Addition of two vectors.
        :param v: Vector.
        :return: Result.
        """

        return Vect(self.X + v.X, self.Y + v.Y, self.Z + v.Z)

    # ----------------------------------------------------------------------------------------------

    def __sub__(self, v):
        """
        Subtraction of two vectors.
        :param v: Vector.
        :return: Result.
        """

        return Vect(self.X - v.X, self.Y - v.Y, self.Z - v.Z)

    # ----------------------------------------------------------------------------------------------

    def __mul__(self, k):
        """
        Multiplication vector on number.
        :param k: Value.
        :return: Result (vector).
        """

        return Vect(self.X * k, self.Y * k, self.Z * k)

    # ----------------------------------------------------------------------------------------------

    def __truediv__(self, k):
        """
        Division on float value.
        :param k: Value.
        :return: Result (vector).
        """

        return self * (1.0 / k)

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def dot_product(a, b):
        """
        Dot product of two vectors.
        :param a: First vector.
        :param b: Second vector.
        :return: Dot product.
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
        :return: Distance.
        """

        return (self - v).mod()

    # ----------------------------------------------------------------------------------------------

    def is_near(self, v, eps):
        """
        Check if vector near another vector.
        :param v: Vector.
        :param eps: Epsilon.
        :return: true - if near, false - otherwise.
        """

        return self.dist_to(v) < eps

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def angle_cos(a, b):
        """
        Cosine of angle between vectors.
        :param a: First vector.
        :param b: Second vector.
        :return: Angle cosine.
        """

        return Vect.dot_product(a, b) / (a.mod() * b.mod())

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def cross_product(a, b):
        """
        Cross product of two vectors.
        :param a: First vector.
        :param b: Second vector.
        :return: Cross product.
        """

        return Vect(a.Y * b.Z - a.Z * b.Y,
                    a.Z * b.X - a.X * b.Z,
                    a.X * b.Y - a.Y * b.X)

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
