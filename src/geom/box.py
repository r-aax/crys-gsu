"""
Box realization.
"""

import math
from geom.vect import Vect

# ==================================================================================================


class Box:
    """
    Box realization.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, ps):
        """
        Constructor from points.
        :param ps: Points.
        """

        xs = [p.X for p in ps]
        ys = [p.Y for p in ps]
        zs = [p.Z for p in ps]

        # Init.
        self.MinX = min(xs)
        self.MinY = min(ys)
        self.MinZ = min(zs)
        self.MaxX = max(xs)
        self.MaxY = max(ys)
        self.MaxZ = max(zs)

    # ----------------------------------------------------------------------------------------------

    def __str__(self):
        """
        String representation.
        :return: String.
        """

        return 'Box: X({0} - {1}), Y({2} - {3}), Z({4} - {5})'.format(self.MinX, self.MaxX,
                                                                      self.MinY, self.MaxY,
                                                                      self.MinZ, self.MaxZ)

    # ----------------------------------------------------------------------------------------------

    def is_inside(self, p):
        """
        Check if point inside box.
        :param p: Point.
        :return: True - if point is inside box, False - otherwise.
        """

        return (p.X >= self.MinX) and (p.X <= self.MaxX) \
               and (p.Y >= self.MinY) and (p.Y <= self.MaxY) \
               and (p.Z >= self.MinZ) and (p.Z <= self.MaxZ)

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================

