"""
Triangles cloud realization.
"""

import math
from geom.segment import Segment
from geom.triangle import Triangle

# ==================================================================================================


class TrianglesCloud:
    """
    Triangles cloud realization.
    Triangles cloud is container for triangles objects.
    It can contain triangles in self.Triangles field or subclouds in self.Subclouds list
      (BUT NOT IN BOTH).
    Each subcloud in self.Subclouds list is an instance of class TrianglesCloud.

    Example 1:
      One solid cloud of 4 triangles and no subclouds.

      TrianglesCloud:
        Triangles = [t1, t2, t3, t4]
        Subclouds = []

    Example 2:
      Cloud is separated into binary tree of subclouds.

                                          TrianglesCloud([])
                                                | |
                         *----------------------* *----------------------*
                         |                                               |
                         V                                               V
                  TrianglesCloud([])                              TrianglesCloud([])
                        | |                                             | |
             *----------* *----------*                       *----------* *----------*
             |                       |                       |                       |
             V                       V                       V                       V
      TrianglesCloud([t1])    TrianglesCloud([t2])    TrianglesCloud([t3])    TrianglesCloud([t4])
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, triangles_list):
        """
        Constructor.
        :param triangles_list: List of triangles.
        """

        # List of triangles.
        self.Triangles = triangles_list

        # List of children clouds.
        self.Subclouds = []

    # ----------------------------------------------------------------------------------------------

    def is_list(self):
        """
        Check cloud is list.
        :return: True - if cloud is list, False - otherwise.
        """

        is_triangles = (self.Triangles != [])
        is_subclouds = (self.Subclouds != [])

        if is_triangles and is_subclouds:
            raise Exception('internal error')

        return is_triangles

    # ----------------------------------------------------------------------------------------------

    def first_intersection_with_segment(self, s):
        """
        Find first intersection with segment.
        :param s: Segment.
        :return:  Triangle (if there is intersection) or None.
        """

        # Processing of the list.
        # If cloud contains triangles in Triangles field, just check all of them.
        if self.is_list():
            for t in self.Triangles:
                if t.intersection_with_segment(s) != []:
                    return t
            return None

        # TODO:
        #   process subclouds.

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
