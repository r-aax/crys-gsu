"""
Triangles cloud realization.
"""

# import math
import numpy as np
# from geom.segment import Segment
# from geom.triangle import Triangle
# from geom.vect import Vect
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

    # ----------------------------------------------------------------------------------------------

    def intersection_with_segment(self, s):
        """
        Find first intersection with segment.
        :param s: Segment.
        :return: Triangle (if there is intersection) or None.
        """

        # предподгоотовка: основа дерева
        self.Subclouds = [self.Triangles]

        # основной процесс: проверяются ветви до тех пор, пока они есть
        while self.Subclouds:

            arr = self.Subclouds[0]

            # формирование куба массива
            self.boxarr = Box.from_triangles(arr)

            # проверка вхождений вектора в куб массива
            if self.boxarr.is_potential_intersect_with_segment(s):

                # проверка на тип элемента дерва
                if len(arr) > 1:

                    # ветка
                    # отращиваем новые ветви
                    self.get_branches(arr)

                else:

                    # лист
                    # проверка на пересечение вектора и листа
                    if arr[0].intersection_with_segment(s) != []:
                        return arr[0]

            # уничтожение ветви
            self.Subclouds.pop(0)

    # ----------------------------------------------------------------------------------------------

    def get_branches(self, arr):
        """
        Make new branches

        :param arr: cloud coordinates of triangles
        :return: append self.Subclouds
        """

        # получаем новые ветви
        branches = self.tree(arr)

        # формируем левую ветвь
        if len(branches[0]) != 0:
            self.Subclouds.append(branches[0])
        else:
            pass

        # формируем правую ветвь
        if len(branches[1]) != 0:
            self.Subclouds.append(branches[1])
        else:
            pass

    # ----------------------------------------------------------------------------------------------

    def tree(self, mass):
        """

        :param mass: список треугольников
        :return: список из двух списков треугольников
        """

        if len(mass) == 2:
            arr_left = [mass[0]]
            arr_right = [mass[1]]
            res_arr = [arr_left, arr_right]
            return res_arr

        else:
            # краевые точки box
            xmax, ymax, zmax = self.boxarr.MaxX, self.boxarr.MaxY, self.boxarr.MaxZ
            xmin, ymin, zmin = self.boxarr.MinX, self.boxarr.MinY, self.boxarr.MinZ

            # проверка длинной стороны
            lenxyz = [xmax - xmin, ymax - ymin, zmax - zmin]
            indxyz = lenxyz.index(np.amax(lenxyz))

            # вычисление центра длинной стороны
            sumxyz = [xmax + xmin, ymax + ymin, zmax + zmin]
            mid_surf = sumxyz[indxyz] / 2

            # бинарное разбиение массива относительно центральной точки
            arr_left = [t for t in mass if t.centroid()[indxyz] < mid_surf]
            arr_right = [t for t in mass if t.centroid()[indxyz] >= mid_surf]

            # проверка корректности разбиения
            if len(arr_left) == 0 and len(arr_right) > 2:
                arr_left = [arr_right[0]]
                arr_right = arr_right[1:]
            elif len(arr_right) == 0 and len(arr_left) > 2:
                arr_right = [arr_left[0]]
                arr_left = arr_left[1:]

            res_arr = [arr_left, arr_right]
            return res_arr

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
