"""
Box realization.
"""

import math
from geom.vect import Vect
from geom.triangle import Triangle

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

    @staticmethod
    def from_points(ps):
        """
        Constructor from points list.

        :param ps: Points list.
        :return:   Box.
        """

        return Box(ps)

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def from_triangles(ts):
        """
        Constructor from triangles list.

        :param ts: Triangles list.
        :return:   Box.
        """

        # Extract points from all triangles (2-dimensional list).
        pss = [t.Points for t in ts]

        # Merge lists and create box from merged points list.
        return Box([p for ps in pss for p in ps])

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

    def is_point_inside(self, p):
        """
        Check if point inside box.

        :param p: Point.
        :return:  True - if point is inside box, False - otherwise.
        """

        return (p.X >= self.MinX) and (p.X <= self.MaxX) \
               and (p.Y >= self.MinY) and (p.Y <= self.MaxY) \
               and (p.Z >= self.MinZ) and (p.Z <= self.MaxZ)

    # ----------------------------------------------------------------------------------------------

    def is_potential_intersect_with_segment(self, s):
        """
        Check if there is intersection with segment.

        :param s: Segment.
        :return:
            True - if there is potential intersection with segment,
            False - there is not intersection for sure.
        """

        # Get segment box.
        sb = Box.from_points(s.Points)

        # Check no intersection along coordinates.
        is_no_x = (sb.MinX > self.MaxX) or (sb.MaxX < self.MinX)
        is_no_y = (sb.MinY > self.MaxY) or (sb.MaxY < self.MinY)
        is_no_z = (sb.MinZ > self.MaxZ) or (sb.MaxZ < self.MinZ)

        if is_no_x or is_no_y or is_no_z:
            # There is definitely no intersection.
            return False
        else:
            # Maybe there is intersection, but maybe not.
            # Maybe there is some sense to check intersection in all rest cases.
            return True

    # ----------------------------------------------------------------------------------------------

    def is_potential_intersect_with_box(self, other_box):
        """
        Intersection of two boxes.

        :param other_box: The second box.
        :return: Intersection of boxes.
        """

        # Check no intersection along coordinates.
        is_no_x = (other_box.MinX > self.MaxX) or (other_box.MaxX < self.MinX)
        is_no_y = (other_box.MinY > self.MaxY) or (other_box.MaxY < self.MinY)
        is_no_z = (other_box.MinZ > self.MaxZ) or (other_box.MaxZ < self.MinZ)

        if is_no_x or is_no_y or is_no_z:
            # There is definitely no intersection.
            return False
        else:
            # Maybe there is intersection, but maybe not.
            # Maybe there is some sense to check intersection in all rest cases.
            return True

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
