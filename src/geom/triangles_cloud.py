"""
Triangles cloud realization.
"""

import numpy as np
import itertools
from geom.box import Box

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

    # Maximum count of triangles in list.
    MaxListTrianglesCount = 1

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

        # Box.
        self.Box = Box.from_triangles(self.Triangles)

        # Build subclouds tree.
        self.build_subclouds_tree()

    # ----------------------------------------------------------------------------------------------

    def build_subclouds_tree(self):
        """
        Build subclouds tree.
        """

        if len(self.Triangles) > TrianglesCloud.MaxListTrianglesCount:

            # Separate triangles list and buld new subtrees.
            # self.Triangles must be cleaned up.
            new_tri_lists = self.separate_triangles_list(self.Triangles)
            self.Triangles = []
            self.Subclouds = [TrianglesCloud(li) for li in new_tri_lists]

        else:

            # Do nothings.
            # Triangles stay triangels.
            pass

    # ----------------------------------------------------------------------------------------------

    def separate_triangles_list(self, mash):
        """
        Separete triangles list into pair of lists.
        :param mash: List of triangles.
        :return: A list of two lists of triangles.
        """

        assert len(mash) > 1, 'internal error'

        # Edge points box.
        xmax, ymax, zmax = self.Box.MaxX, self.Box.MaxY, self.Box.MaxZ
        xmin, ymin, zmin = self.Box.MinX, self.Box.MinY, self.Box.MinZ

        # Checking the long side.
        lenxyz = [xmax - xmin, ymax - ymin, zmax - zmin]
        indxyz = lenxyz.index(np.amax(lenxyz))

        # Calculating the center of the long side.
        sumxyz = [xmax + xmin, ymax + ymin, zmax + zmin]
        mid_surf = sumxyz[indxyz] / 2

        # Binary partitioning of an array relative to the center point.
        arr_left = [t for t in mash if t.centroid()[indxyz] < mid_surf]
        arr_right = [t for t in mash if t.centroid()[indxyz] >= mid_surf]

        # Checking the correctness of the split.
        if len(arr_left) == 0:
            arr_left = [arr_right[0]]
            arr_right = arr_right[1:]
        elif len(arr_right) == 0:
            arr_right = [arr_left[0]]
            arr_left = arr_left[1:]

        res_arr = [arr_left, arr_right]
        return res_arr

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
        Find intersections with segment.
        :param s: Segment.
        :return:  Triangle (if there is intersection) or None.
        """

        # Cold check.
        if not self.Box.is_potential_intersect_with_segment(s):
            return None

        if self.is_list():
            for t in self.Triangles:
                if t.intersection_with_segment(s) != []:
                    return t
        else:
            for sc in self.Subclouds:
                res = sc.first_intersection_with_segment(s)
                if res != None:
                    return res

        # No intersection is found.
        return None

    # ----------------------------------------------------------------------------------------------

    def intersection_with_triangles_cloud(self, tc):
        """
        Find intersections with another triangle cloud.
        :param tc: triangle cloud
        :return: [] or list of triangle pairs - [t1, t2]
        """

        # Cold check.
        if not self.Box.is_potential_intersect_with_box(tc.Box):
            return []

        if self.Subclouds != []:
            return list(itertools.chain(*[a.intersection_with_triangles_cloud(tc) for a in self.Subclouds]))

        elif tc.Subclouds != []:
            return list(itertools.chain(*[self.intersection_with_triangles_cloud(b) for b in tc.Subclouds]))

        else:
            return [[t1, t2]
                   for t1 in self.Triangles
                   for t2 in tc.Triangles
                   if t1.intersection_with_triangle(t2) != []]


# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
