"""
refine grid from file
"""

import os
import time
from gsu.gsu import Grid
from geom.triangle import Triangle
from geom.triangles_cloud import TrianglesCloud

# TODO реализовать метод перестроения сетки - перестроение пересекающихся треугольников


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

    # Check for grid file.
    if not os.path.isfile(grid_file):
        raise Exception(f'crys-gsu-refine : no such file ({grid_file})')

    print(f'crys-gsu-refine : start, grid_file = {grid_file}, '
          f'out_grid_file = {out_grid_file}')

    start_time = time.time()

    # Load grid.
    g = Grid()
    g.load(grid_file)

    # поиск пересечений элементов сетки (треугольников)
    triangles_list = g.get_triangles_list()
    triangles_cloud = TrianglesCloud(triangles_list)
    res_intersect = triangles_cloud.intersection_with_triangles_cloud(triangles_cloud)

    for tr1, tr2 in res_intersect:

        # intersection points
        intersection_points = tr1.intersection_with_triangle(tr2)

        # analysis of each triangle
        for tri in [tr1, tr2]:
            points_in_triangle = []

            # analyze the intersection points for the triangle
            for point in intersection_points:
                position_in_the_triangle = tri.point_in_triangle(point)
                if position_in_the_triangle != 0 and position_in_the_triangle != 1:
                    points_in_triangle.append([position_in_the_triangle, point])

            if points_in_triangle == []:
                # do not rebuild this triangle
                pass

            else:

                if len(points_in_triangle) == 1:

                    if points_in_triangle[0][0] == 2:
                        res_pos = point_on_edge_of_triangle(points_in_triangle[0][1], tri)
                        # TODO перестроить этот треугольник внутри сетки
                        # return [Triangle(А, B, P2), Triangle(B, C, P2)]
                        pass
                    else:
                        # TODO построить три новых треугольника в сетке на месте старого
                        # return [Triangle(А, B, P3), Triangle(B, C, P3) , Triangle(А, C, P3)]
                        pass

                else:

                    if points_in_triangle[0][0] == 2 and points_in_triangle[1][0] == 2:
                        res_pos1 = point_on_edge_of_triangle(points_in_triangle[0][1], tri)
                        res_pos2 = point_on_edge_of_triangle(points_in_triangle[0][1], tri)

                        # belong the same edge
                        if functools.reduce(lambda x, y : x and y, map(lambda p, q: p == q, res_pos1[0], res_pos2[0]),
                                            True):
                            # TODO построить три новых треугольника в сетке на месте старого
                            pass
                        else:
                            # TODO построить три новых треугольника в сетке на месте старого
                            pass

                    elif points_in_triangle[0][0] == 3 and points_in_triangle[1][0] == 3:
                        # TODO построить три новых треугольника в сетке на месте старого
                        pass

                    else:
                        # TODO построить три новых треугольника в сетке на месте старого
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