import argparse
from os.path import isfile
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from time import time
from scipy.spatial import KDTree
from numpy import array
import gsu


def check_argument(filename):
    """Check existance and naming of meshes.

    Parameters
    ----------
    filename : filename
    """
    if not isfile(filename):
        print('File {} does not exist'.format(filename))
        exit(1)
    else:
        if not filename[-4:] == '.dat':
            print('File {} should be .dat file'.format(filename))
            exit(1)

def all_coordinates_of_aux_nodes(mesh):
    """Compose array of coordinates of aux nodes of size [N, 3]

    Parameters
    ----------
      mesh : mesh with faces having aux_node field (Vect)

    Returns
    -------
      array
    """
    x, y, z = list(), list(), list()

    for f in mesh.Faces:
        x.append(f.aux_node.X)
        y.append(f.aux_node.Y)
        z.append(f.aux_node.Z)

    return array([x, y, z]).T


def dual_mesh(mesh):
    """Compose dual mesh.

    Parameters
    ----------
      mesh : triangular mesh
    """
    for f in mesh.Faces:
        f.aux_node = f.get_triangle().centroid()


def face_centered_interpolation(source_mesh, target_mesh):
    """Interpolate between meshes by 1-nearest neigbhour search.
    
    Parameters
    ----------
      source_mesh : source mesh
      target_mesh : target mesh
    """
    dual_mesh(source_mesh)
    dual_mesh(target_mesh)
    old_aux_nodes = all_coordinates_of_aux_nodes(source_mesh)

    kdtree = KDTree(old_aux_nodes)
    for f in target_mesh.Faces:
        i = kdtree.query(f.aux_node.coords_as_np())[1]
        f['T'] = source_mesh.Faces[i]['T']
        f['Hw'] = source_mesh.Faces[i]['Hw']


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source', help='source mesh .dat file to interpolate from')
    parser.add_argument('target', help='target_mesh .dat file to interpolate to')
    parser.add_argument('-res', '--result_grid', help='interpolated grid. if not provided than the name of the '
                                                    'result file is \"new_grid\" + \"_interpolated\"')
    parser.add_argument("-v", "--verbosity", action="count",
                        help="increase output verbosity", default=0)
    args = parser.parse_args()

    source_mesh = args.source
    target_mesh = args.target
    result_grid = args.result_grid

    check_argument(source_mesh)
    check_argument(target_mesh)
    if result_grid:
        if not result_grid[-4:] == '.dat':
            print('File {} should be .dat file'.format(result_grid))
            exit(1)

    start = time()
    grid1 = gsu.Grid()
    grid2 = gsu.Grid()
    grid1.load(source_mesh)

    if args.verbosity > 0:
        print('Old grid read')

    grid2.load(target_mesh)

    if args.verbosity > 0:
        print('New grid read')

    face_centered_interpolation(grid1, grid2)

    if args.verbosity > 0:
        print('Interpolation made')

    if args.result_grid:
        grid2.store(result_grid)
    else:
        grid2.store(target_mesh[:-4] + '_interpolated.dat')

    if args.verbosity > 0:
        print('Result grid was written')
        print('Total time:', time() - start)
