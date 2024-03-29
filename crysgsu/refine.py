"""
refine grid from file
"""

import os
import time
import functools
import math
from gsu.gsu import Grid
from gsu.node import Node
from gsu.edge import Edge
from gsu.face import Face
from geom.vect import Vect
from geom.triangle import Triangle
from geom.triangles_cloud import TrianglesCloud


def refine_grid(grid_file, out_grid_file):
    """

    Parameters
    ----------
    grid_file: file with grid
    out_grid_file: out \*.dat file

    Returns: saved a copy of the grid file
    -------

    """

    # TODO добвить функционал перестроения плохих треугольников
    #     (с углами либо тремящимися к 0, либо к 180 градусам)

    def point_on_edge_of_triangle(point, tri):
        """
        The method finds the edge to which the point belongs and the opposite vertex of the triangle.

        Parameters
        ----------
        point: Vect object
        tri: Triangle object

        Returns: [[Vect, Vect], Vect]
        -------

        """

        # Get all points.
        a, b, c = tri.a(), tri.b(), tri.c()

        if not Triangle(a, b, point).area():
            return [[a, b], c]
        elif not Triangle(point, b, c).area():
            return [[b, c], a]
        else:
            return [[c, a], b]

    # ----------------------------------------------------------------------------------------------

    def in_point_class(arg, points_list):
        """

        Parameters
        ----------
        arg: class number (0, 1, 2 or 3)
        points_list: list of point lists and their class numbers

        Returns: True or False
        -------

        """

        for i in points_list:
            if arg == i[0]:
                return True
        return False

    # ----------------------------------------------------------------------------------------------

    def get_node_from_list_of_nodes(vect, nodes):
        """

        Parameters
        ----------
        vect: Vect
        nodes: List of Nodes

        Returns: Node or None
        -------

        """

        node = Node(vect)
        if node in nodes:
            return nodes[nodes.index(node)]

        return None

    # ----------------------------------------------------------------------------------------------

    def refine(g, intersections_list):
        """

        Данный метод перестраивает сетку на основе списка пересечений.

        Parameters
        ----------
        g: Сетка
        intersections_list: Список пересечений,
        элементы которого предстваленны в виде списка пар [треугольник, [точки пересечений]]

        Returns: Перестроенная сетка g
        -------

        """

        # Поочереди перестраиваем все треугольники (если у них есть точки для перестроения)
        for tri, intersection_points in intersections_list:

            # список точек в треугольнике
            rebuilding_points = []

            # объект сетки для перестроения
            face = tri.BackRef
            assert (face)

            # analyze the intersection points for the triangle
            for point in intersection_points:
                position_in_the_triangle = tri.point_in_triangle(point)
                if position_in_the_triangle != 0 and position_in_the_triangle != 1:
                    rebuilding_points.append([position_in_the_triangle, point])

            if not rebuilding_points:
                # do not rebuild this triangle
                pass

            else:

                # перестроения треугольника при 1 точке
                if len(rebuilding_points) == 1:

                    if rebuilding_points[0][0] == 2:
                        # перестроение треугольника по одной точке Р2 на два новых треугольника
                        # чтобы воспользоваться атомарным методом "Разбиение ребра", необходимо:
                        # 1) определить какому ребру принадлежит точка
                        # 2) есть ли с обратной стороны ребра еще один face
                        res_pos = point_on_edge_of_triangle(rebuilding_points[0][1], tri)
                        # вершины ребра
                        n1 = res_pos[0][0]
                        n2 = res_pos[0][1]
                        node1 = Node(n1)
                        node2 = Node(n2)
                        ed = [ed for ed in face.Edges if (ed.Nodes[0].P == node1.P and ed.Nodes[1].P == node2.P) or
                              (ed.Nodes[1].P == node1.P and ed.Nodes[0].P == node2.P)][0]

                        if len(ed.Faces) == 2:
                            g.cut_edge(ed, rebuilding_points[0][1])
                        else:
                            g.cut_single_edge(ed, rebuilding_points[0][1])

                    else:
                        # перестроение треугольника по одной точке Р3 на 3 новых треугольника
                        g.divide_face(face, rebuilding_points[0][1])
                        pass

                # перестроения треугольника при 2 точках
                elif len(rebuilding_points) == 2:

                    if rebuilding_points[0][0] == 2 and rebuilding_points[1][0] == 2:
                        res_pos1 = point_on_edge_of_triangle(rebuilding_points[0][1], tri)
                        res_pos2 = point_on_edge_of_triangle(rebuilding_points[1][1], tri)

                        # belong the same edge
                        if functools.reduce(lambda x, y: x and y, map(lambda p, q: p == q, res_pos1[0], res_pos2[0]),
                                            True):

                            # вершины ребра
                            n1 = res_pos1[0][0]
                            n2 = res_pos1[0][1]
                            node1 = Node(n1)
                            node2 = Node(n2)
                            ed = [ed for ed in face.Edges if (ed.Nodes[0].P == node1.P and ed.Nodes[1].P == node2.P) or
                                  (ed.Nodes[1].P == node1.P and ed.Nodes[0].P == node2.P)][0]
                            g.cut_edge_with_two_nodes(ed, rebuilding_points[0][1], rebuilding_points[1][1])

                        else:
                            # TODO построить три новых треугольника в сетке на месте старого (две точки на двух ребрах)
                            pass

                    elif rebuilding_points[0][0] == 3 and rebuilding_points[1][0] == 3:
                        # TODO построить три новых треугольника в сетке на месте старого (обе точки в центре)
                        pass

                    else:
                        # TODO построить три новых треугольника в сетке на месте старого (одна в центре, вторая на ребре)
                        pass
                    pass

                # # перестроения треугольника при 3 точках
                # elif len(rebuilding_points) == 3:
                #     # TODO перестроение треугольника при 3 точках
                #     pass

                # перестроения треугольника при прочих количествах точек
                else:
                    # TODO перестроение треугольника при прочих количествах точек
                    if in_point_class(3, rebuilding_points):
                        # если len(rebuilding_points) > 3, то разбиваем треугольник в сетке на 3 треугольника по точке,
                        # ближайшей к центру треугольника
                        # выводим их в виде объектов класса Треугольник в конец списка res_intersect_new_list как,
                        # res_intersect_new_list.append([tri1, rebuilding_points],
                        #                               [tri2, rebuilding_points],
                        #                               [tri3, rebuilding_points])
                        # где rebuilding_points - набор точек, которые однозначно пересекали исходный треугольник
                        # так новые треугольники попадут в конец цикла и тоже будут рассмотрены и разбиты на части
                        pass
                    else:
                        # нет ни одной точки в центре треугольника => есть много точек на ребрах (больше трех)
                        # нужно разработать алгоритм действий в этом случае
                        pass

    # =================================================================================================

    # Check for grid file.
    if not os.path.isfile(grid_file):
        raise Exception(f'crys-gsu-refine : no such file ({grid_file})')

    print(f'crys-gsu-refine : start, grid_file = {grid_file}, '
          f'out_grid_file = {out_grid_file}')

    start_time = time.time()

    # Load grid.
    g = Grid()
    g.load(grid_file)

    # finding intersections of grid elements (triangles)
    triangles_list = g.get_triangles_list()
    triangles_cloud = TrianglesCloud(triangles_list)
    res_intersect = triangles_cloud.intersection_with_triangles_cloud(triangles_cloud)

    # необходимо отсортировать по треугольникам и объединить точки пересечений по одному треугольнику в один список
    # сортировка списка треугольников по id
    res_intersect.sort(key=lambda list_of_tri_and_points: id(list_of_tri_and_points[0]))
    # первый элемент просто записываем в новый список
    res_intersect_new_list = [res_intersect[0]]
    # остальные элементы сравниваем с предыдущими
    # если треугольники одинаковвые, то к предыдущему дописывем список точек
    # если треугольники разные - дописываем в новый список треугольник с точками
    for i in range(1, len(res_intersect)):
        if id(res_intersect[i][0]) != id(res_intersect[i-1][0]):
            res_intersect_new_list.append(res_intersect[i])
        else:
            res_intersect_new_list[-1][1] += res_intersect[i][1]

    # одинаковые точки нужно отсеять
    for j in res_intersect_new_list:
        j[1].sort()
        j[1] = [j[1][0]] + [j[1][i] for i in range(1, len(j[1])) if
                           not (math.isclose(j[1][i].X, j[1][i - 1].X) and
                                math.isclose(j[1][i].Y, j[1][i - 1].Y) and
                                math.isclose(j[1][i].Z, j[1][i - 1].Z))]

    # Поочереди перестраиваем все треугольники (если у них есть точки для перестроения)
    refine(g, res_intersect_new_list)

    g.store(out_grid_file)

    print(f'crys-gsu-refine : done (time estimated = {time.time() - start_time} s)')

# ========================================================================================================


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='refine',
                                     description='The tool for improving/rebuilding the grid.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('grid_file', help='grid file name')
    parser.add_argument('out_grid_file', help='out file with result grid')
    parser.add_argument('-v', '--verbosity', action='count',
                        help='increase output verbosity', default=0)
    args = parser.parse_args()

    refine_grid(args.grid_file, args.out_grid_file)
