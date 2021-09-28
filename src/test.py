"""
Test module.
"""
from datetime import datetime
import time
import os
import string

import gsu
import split


# --------------------------------------------------------------------------------------------------


def case_001_node_face_data(test='wing_1_mz'):
    """
    Load small grid and store it without any changes.
    Test objective:
      To make sure that node data and face data are correct.
      Node data is just a point.
      Face data is arbitrary array.
    :param test: test name
    """

    print('case_001_node_face_data({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.store('grids/{0}_original.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_002_decompose(test='bunny'):
    """
    Test all decompose algorithms.
    Test objective:
      To make sure that all grid decomposition algorithms work.
    :param test: test name
    """

    print('case_002_decompose({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))

    # mono
    g.decompose_mono(new_name=test + ' mono')
    g.store('grids/{0}_mono.dat'.format(test))
    # random
    g.decompose_random(new_name=test + ' random')
    g.store('grids/{0}_random.dat'.format(test))
    # linear
    g.decompose_linear(new_name=test + ' linear')
    g.store('grids/{0}_linear.dat'.format(test))

    # hierarchical
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                 gsu.fun_face_cy(),
                                                 gsu.fun_face_cz()],
                             levels=6,
                             new_name=test + ' hierarchical')
    g.store('grids/{0}_hierarchical.dat'.format(test))

    # pressure
    g.decompose_pressure(new_name=test + ' pressure')
    g.store('grids/{0}_pressure.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_003_fixed_zones(test='bunny_pos'):
    """
    Fixed zones mechanism.
    Test objective:
      To make sure that zone marked with flag 'Fixed' stay unchanged.
      This feature is implemented only for hierarchical algorithm.
    :param test: test name
    """

    print('case_003_fixed_zones({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                 gsu.fun_face_cy(),
                                                 gsu.fun_face_cz()],
                             levels=6,
                             new_name=test + ' hierarchical',
                             fixed_zones=['POS1', 'POS2'])
    g.store('grids/{0}_hierarchical.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_004_load_store_load(test='bunny'):
    """
    Load grid, save it immediately, and load again.
    Test objective:
      To make sure that operation load-store doesn't corrupt the grid.
    :param test: test name
    """

    print('case_004_load_store_load({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.store('grids/{0}_c1.dat'.format(test))
    g.load('grids/{0}_c1.dat'.format(test))
    g.store('grids/{0}_c2.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_005_explode_bunny(test='bunny_pos'):
    """
    Visual bunny explosion (illusion of different zones run away from each other).
    Test objective:
      To check how we can move zones in different directions.
    :param test: test name
    """

    print('case_005_explode_bunny({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                 gsu.fun_face_cy(),
                                                 gsu.fun_face_cz()],
                             levels=5,
                             new_name=test + ' hierarchical',
                             fixed_zones=['POS1', 'POS2'])
    g.store('grids/{0}_hierarchical.dat'.format(test))
    g.load('grids/{0}_hierarchical.dat'.format(test), is_merge_same_nodes=False)
    g.move_from_mean_point(0.25)
    g.store('grids/{0}_explode.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_006_store_faces_t_hw_hi(test='bunny'):
    """
    Store faces T, Hw, Hi data.
    Test objective:
      Storing data ot T, Hw, Hi to file (emulation the result of MPI program).
    :param test: test name
    """

    print('case_006_store_faces_t_hw_hi({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.store_faces_calc_data('grids/{0}.txt'.format(test))


# --------------------------------------------------------------------------------------------------


def case_007_load_faces_t_hw_hi(test='bunny'):
    """
    Load faces T, Hw, Hi from file.
    Test objective:
      Check how we can replace data T, Hw, Hi in grid from file.
    :param test: test name
    """

    print('case_007_load_faces_t_hw_hi({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.store_faces_calc_data('grids/{0}.txt'.format(test))
    g.load_faces_calc_data('grids/{0}.txt'.format(test))
    g.store('grids/{0}_data.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_008_decompose_pressure(test='lucy'):
    """
    Load grid and decompose with pressure algorithm.
    Test objective:
      Check decompose pressure algorithm.
    :param test: test name
    """

    print('case_008_decompose_pressure({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.decompose_pressure(new_name=test + ' pressure', fz_names=['POS1', 'POS2'])
    g.print_info(True, True, True)
    g.store('grids/{0}_pressure.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_009_store_mpi(test='bunny'):
    """
    Load grid, decompose it and store in MPI.
    Test objective:
      Check we can store data for multiprocessing mode calculations.
    :param test: test name
    """

    print('case_009_store_mpi({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                 gsu.fun_face_cy(),
                                                 gsu.fun_face_cz()],
                             levels=3,
                             new_name=test + ' hierarchical')
    g.store_mpi('grids/mpi', '2021-08-11-12-00-00')


# --------------------------------------------------------------------------------------------------

def case_010_decompose_incremental(test='bunny'):
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))

    g.decompose_incremental(32)
    g.print_info(True, True)
    print(g.get_metric())
    g.store('grids/{0}_incremental.dat'.format(test))


def case_011_decompose_incremental_vs_local_refinement(test='bunny'):
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))

    g.decompose_incremental(32)
    g.print_info(True, True)
    print(g.get_metric())
    g.store('grids/{0}_incremental.dat'.format(test))

    g.local_refinement()
    g.print_info(True, True)
    print(g.get_metric())

    g.store('grids/{0}_incremental_refined.dat'.format(test))


def case_012_decompose_pressure_vs_local_refinement(test='bunny'):
    print('case_008_decompose_pressure({0}):'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.decompose_pressure(new_name=test + ' pressure', fz_names=['POS1', 'POS2'])
    g.print_info(True, True)
    print(g.get_metric())
    g.store('grids/{0}_pressure.dat'.format(test))

    g.local_refinement()
    g.print_info(True, True)
    print(g.get_metric())

    g.store('grids/{0}_pressure_refined.dat'.format(test))


def case_013_get_notebook_stat():
    grids = ['bunny']  # , 'dragon', 'lucy']
    for i in range(0, 100):
        for grid in grids:
            process_grid(grid, i)


def case_014_decompose_incremental_consistently(test='bunny'):
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))

    g.decompose_incremental_consistently(32)
    g.print_info(True, True)
    print(g.get_metric())
    g.store('grids/{0}_incremental.dat'.format(test))


def process_grid(grid, ii):
    cross_edges = []
    max_border = []
    fd = []
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(grid))
    seconds = datetime.now()

    for i in range(2, 64):
        print('--------------------' + datetime.now().strftime("%H:%M:%S") + ' grid ' + grid + ' zones ' + str(
            i) + ' attempt ' + str(ii) + ', prev take {0}  s --------------------'.format(datetime.now() - seconds))
        g.decompose_incremental(i)

        zam = gsu.ZonesAdjacencyMatrix(g.Edges, g.Zones)
        cross_edges.append(zam.edges_statistics()[6])
        max_border.append(zam.max_cross_border_len())
        faces_distr_dev, max_zone_faces_count = g.get_fd_deviation()
        fd.append(faces_distr_dev)

        g.local_refinement()

        zam = gsu.ZonesAdjacencyMatrix(g.Edges, g.Zones)
        cross_edges.append(zam.edges_statistics()[6])
        max_border.append(zam.max_cross_border_len())
        fd.append(g.get_fd_deviation())
        seconds = datetime.now()

    print(grid)
    print(grid + '_inc_ce=[' + ','.join(str(a) for a in cross_edges) + ']')
    print(grid + '_inc_mb=[' + ','.join(str(a) for a in max_border) + ']')
    print(grid + '_inc_fd=[' + ','.join(str(a) for a in fd) + ']')

    f = open("{0}-{1}.txt".format(grid, ii), "a")
    f.write(grid + '_inc_ce=[' + ','.join(str(a) for a in cross_edges) + ']')
    f.write('\r\n')
    f.write(grid + '_inc_mb=[' + ','.join(str(a) for a in max_border) + ']')
    f.write('\r\n')
    f.write(grid + '_inc_fd=[' + ','.join(str(a) for a in fd) + ']')
    f.close()


if __name__ == '__main__':
    # case_001_node_face_data()
    # case_002_decompose()
    # case_003_fixed_zones()
    # case_004_load_store_load()
    # case_005_explode_bunny()
    # case_006_store_faces_t_hw_hi()
    # case_007_load_faces_t_hw_hi()
    # case_008_decompose_pressure()
    # case_009_store_mpi()
    # case_010_decompose_incremental()
    # case_011_decompose_incremental_vs_local_refinement()
    # case_012_decompose_pressure_vs_local_refinement()
    # case_013_decompose_hierarchical_vs_local_refinement()
    # case_013_get_notebook_stat()
    case_014_decompose_incremental_consistently()
    pass

# --------------------------------------------------------------------------------------------------
