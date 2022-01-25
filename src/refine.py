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

    def del_face_from_grid_and_get_data_from_face(tri_in_grid, g):
        """

        Parameters
        ----------
        face: Face in Grid
        g: Grid object

        Returns: data, nodes, edges, zone form Face in Grid
        -------

        """
        # 1) временно сохранить Data
        data = tri_in_grid.Data.copy()

        # 2) открыть Edges и в каждом Edge удалить Face этого треугольника
        edges = []
        for edge in tri_in_grid.Edges:
            edge.Faces.remove(tri_in_grid)
            edges.append(edge)

        # 3) открыть Nodes и в каждом Node удалить Face этого треугольника
        nodes = []
        for node in tri_in_grid.Nodes:
            node.Faces.remove(tri_in_grid)
            nodes.append(node)

        # 4) открыть Zone и удалить Face этого треугольника
        zone = tri_in_grid.Zone
        zone.Faces.remove(tri_in_grid)

        # 5) удаление треугольника из сетки
        g.Faces.remove(tri_in_grid)

        return data, nodes, edges, zone

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

    for tri, intersection_points in res_intersect_new_list:

        # список точек в треугольнике
        points_in_triangle = []

        # объект сетки для перестроения
        tri_in_grid = tri.BackRef
        assert(tri_in_grid)

        # analyze the intersection points for the triangle
        for point in intersection_points:
            position_in_the_triangle = tri.point_in_triangle(point)
            if position_in_the_triangle != 0 and position_in_the_triangle != 1:
                points_in_triangle.append([position_in_the_triangle, point])

        if not points_in_triangle:
            # do not rebuild this triangle
            pass

        else:

            # перестроения треугольника при 1 точке
            if len(points_in_triangle) == 1:

                if points_in_triangle[0][0] == 2:
                    # точка находится на ребре у этого треугольника и у еще одного
                    # в этом проходе цикла необходимо перестроить только этот треугольник,
                    # потому что второй или уже перестроен или в очереди на перестроение
                    # с использованием этой же точки

                    # перестроение треугольника по одной точке Р2 на два новых треугольника
                    # TODO перестроить этот треугольник
                    # return [Triangle(А, B, P2), Triangle(B, C, P2)]

                    # удалить Face из Grid и вынуть все необходимые данные
                    data, nodes, edges, zone = del_face_from_grid_and_get_data_from_face(tri_in_grid, g)

                    # создать фейсы
                    f1 = Face(list(data.keys()), list(data.values()))
                    f2 = Face(list(data.keys()), list(data.values()))

                    # Nodes
                    # n4 - новый нод, который необходимо создать
                    # не известно между какими нодами лежит n4
                    n4 = Node(points_in_triangle[0][1])

                    # необходимо это определить
                    list_of_point_for_node, opposite_point = point_on_edge_of_triangle(points_in_triangle[0][1], tri)

                    # зададим однозначные ссылки на конкретные ноды (ноды существующие в сетке)
                    # n1 и n2 - ноды ребра для разбиения
                    # n3 - противоположный нод
                    for node in nodes:
                        if node.P == list_of_point_for_node[0]:
                            n1 = node
                        elif node.P == list_of_point_for_node[1]:
                            n2 = node
                        else:
                            n3 = node

                    # Edges
                    # ищем соответствующие ребра
                    # e1 и е2 - целые ребра
                    # е3 - ребро для разбиения и последующего удаления
                    # е4 и е5 - новые ребра
                    for edge in edges:
                        if n1 in edge.Nodes and n3 in edge.Nodes:
                            e1 = edge
                        elif n2 in edge.Nodes and n3 in edge.Nodes:
                            e2 = edge
                        else:
                            e3 = edge

                    # новые ребра
                    e4 = Edge()
                    e5 = Edge()

                    # TODO если удалить е3 из Grid, то у противоположного треугольника не будет ребра е3,
                    #  т.е. треугольник будет состоять из 4 ребер, а нет трех, что поломает все алгоритмы

                    pass
                else:
                    # перестроение треугольника по одной точке Р3 на три новых треугольника
                    # удалить Face из Grid и вынуть все необходимые данные
                    data, nodes, edges, zone = del_face_from_grid_and_get_data_from_face(tri_in_grid, g)

                    # добавить новые фейсы
                    # создать фейс
                    f1 = Face(list(data.keys()), list(data.values()))
                    f2 = Face(list(data.keys()), list(data.values()))
                    f3 = Face(list(data.keys()), list(data.values()))

                    # добавить и залинковать к ним Node и Edge
                    # Nodes
                    n1 = nodes[0]
                    n2 = nodes[1]
                    n3 = nodes[2]
                    n4 = Node(points_in_triangle[0][1])

                    f1.Nodes += [n1, n2, n4]
                    f2.Nodes += [n2, n3, n4]
                    f3.Nodes += [n1, n3, n4]

                    n1.Faces += [f1, f3]
                    n2.Faces += [f1, f2]
                    n3.Faces += [f2, f3]
                    n4.Faces += [f1, f2, f3]

                    # Edges
                    for edge in edges:
                        if n1 in edge.Nodes and n2 in edge.Nodes:
                            e1 = edge
                        elif n2 in edge.Nodes and n3 in edge.Nodes:
                            e2 = edge
                        else:
                            e3 = edge

                    e4 = Edge()
                    e5 = Edge()
                    e6 = Edge()

                    n1.Edges += [e4]
                    n2.Edges += [e5]
                    n3.Edges += [e6]
                    n4.Edges += [e4, e5, e6]

                    # это исходные ребра, их точки не поменялись
                    # e1.Nodes += []
                    # e2.Nodes += []
                    # e3.Nodes += []
                    # новые ребра
                    e4.Nodes += [n4, n1]
                    e5.Nodes += [n4, n2]
                    e6.Nodes += [n4, n3]

                    e1.Faces += [f1]
                    e2.Faces += [f2]
                    e3.Faces += [f3]
                    e4.Faces += [f1, f3]
                    e5.Faces += [f1, f2]
                    e6.Faces += [f3, f2]

                    # добавить данные в Grid
                    g.Faces += [f1, f2, f3]
                    g.Edges += [e4, e5, e6]
                    # TODO а надо ли добавлять Node в Grid? ведь это чья-то вершина, которая уже там есть
                    # может все к существующему объекту привязывать, а не к новому?
                    if not n4 in g.Nodes:
                        g.Nodes += [n4]
                        g.RoundedCoordsBag.add(n4.RoundedCoords)

                    # линк Zone
                    f1.Zone = zone
                    f2.Zone = zone
                    f3.Zone = zone

                    zone.Faces += [f1]
                    zone.Faces += [f2]
                    zone.Faces += [f3]

                    zone.Edges += [e4, e5, e6]

                    zone.Nodes += [n4]

                    del tri_in_grid

            # перестроения треугольника при 2 точках
            elif len(points_in_triangle) == 2:

                if points_in_triangle[0][0] == 2 and points_in_triangle[1][0] == 2:
                    res_pos1 = point_on_edge_of_triangle(points_in_triangle[0][1], tri)
                    res_pos2 = point_on_edge_of_triangle(points_in_triangle[1][1], tri)

                    # belong the same edge
                    if functools.reduce(lambda x, y : x and y, map(lambda p, q: p == q, res_pos1[0], res_pos2[0]),
                                        True):
                        # TODO построить три новых треугольника в сетке на месте старого (две точки на одном ребре)
                        pass
                    else:
                        # TODO построить три новых треугольника в сетке на месте старого (две точки на двух ребрах)
                        pass

                elif points_in_triangle[0][0] == 3 and points_in_triangle[1][0] == 3:
                    # TODO построить три новых треугольника в сетке на месте старого (обе точки в центре)
                    pass

                else:
                    # TODO построить три новых треугольника в сетке на месте старого (одна в центре, вторая на ребре)
                    pass

            # # перестроения треугольника при 3 точках
            # elif len(points_in_triangle) == 3:
            #     # TODO перестроение треугольника при 3 точках
            #     pass

            # перестроения треугольника при прочих количествах точек
            else:
                # TODO перестроение треугольника при прочих количествах точек
                if in_point_class(3, points_in_triangle):
                    # если len(points_in_triangle) > 3, то разбиваем треугольник в сетке на 3 треугольника по точке,
                    # ближайшей к центру треугольника
                    # выводим их в виде объектов класса Треугольник в конец списка res_intersect_new_list как,
                    # res_intersect_new_list.append([tri1, points_in_triangle],
                    #                               [tri2, points_in_triangle],
                    #                               [tri3, points_in_triangle])
                    # где points_in_triangle - набор точек, которые однозначно пересекали исходный треугольник
                    # так новые треугольники попадут в конец цикла и тоже будут рассмотрены и разбиты на части
                    pass
                else:
                    # нет ни одной точки в центре треугольника => есть много точек на ребрах (больше двух)
                    # нужно разработать алгоритм действий в этом случае
                    pass

    # Save grid.
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
