"""
Test module.
"""

import geom
import gsu
import split

# ==================================================================================================


def case_001_load_store(test='wing_1'):
    """
    Load and store small grid.
    Test objective:

      To make sure the grid can be loaded and stored
      and loading and storing does not corrupt the data.

    :param test: test name
    """

    print('case_001_load_store({0})'.format(test))
    g = gsu.Grid()
    g.load('cases/grids/{0}.dat'.format(test))
    g.store('cases/grids/{0}_store1.dat'.format(test))
    g.load('cases/grids/{0}_store1.dat'.format(test))
    g.store('cases/grids/{0}_store2.dat'.format(test))

# --------------------------------------------------------------------------------------------------


def case_002_decompose(test='bunny_pos',
                       methods=['mono', 'random', 'linear', 'hierarchical', 'farhat']):
    """
    Test all decompose algorithms.
    Test objective:

      To make sure that all grid decomposition algorithms work.

    :param test: test name
    """

    print('case_002_decompose({0}): {1}'.format(test, methods))
    g = gsu.Grid()

    # Create new name of grid and file
    new_grid_name = lambda m: '{0} {1}'.format(test, m)
    new_file_name = lambda m: 'cases/grids/{0}_{1}.dat'.format(test, m)

    if 'mono' in methods:
        g.load('cases/grids/{0}.dat'.format(test))
        g.decompose_mono(new_name=new_grid_name('mono'))
        g.store(new_file_name('mono'))

    if 'random' in methods:
        g.load('cases/grids/{0}.dat'.format(test))
        g.decompose_random(new_name=new_grid_name('random'))
        g.store(new_file_name('random'))

    if 'linear' in methods:
        g.load('cases/grids/{0}.dat'.format(test))
        g.decompose_linear(new_name=new_grid_name('linear'))
        g.store(new_file_name('linear'))

    if 'hierarchical' in methods:
        g.load('cases/grids/{0}.dat'.format(test))
        g.decompose_hierarchical(new_name=new_grid_name('hierarchical'),
                                 fixed_zones=['POS1', 'POS2'])
        g.store(new_file_name('hierarchical'))

    if 'farhat' in methods:
        g.load('cases/grids/{0}.dat'.format(test))
        g.decompose_farhat(new_name=new_grid_name('farhat'),
                           fz_names=['POS1', 'POS2'])
        g.store(new_file_name('farhat'))


# --------------------------------------------------------------------------------------------------


def case_005_explode_bunny(test='bunny'):
    """
    Visual bunny explosion (illusion of different zones run away from each other).
    Test objective:

      To check how we can move zones in different directions.

    :param test: test name
    """

    print('case_005_explode_bunny({0})'.format(test))
    g = gsu.Grid()
    g.load('cases/grids/{0}.dat'.format(test))
    g.decompose_hierarchical(levels=5, new_name=test + ' hierarchical')
    g.store('cases/grids/{0}_hierarchical.dat'.format(test))
    g.load('cases/grids/{0}_hierarchical.dat'.format(test), is_merge_same_nodes=False)
    g.move_from_mean_point(0.25)
    g.store('cases/grids/{0}_explode.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_007_store_load_faces_calc_data(test='bunny'):
    """
    Store and load faces calc data.
    Test objective:

      Check how we can replace data T, Hw, Hi and other in grid from file.

    :param test: test name
    """

    print('case_007_store_load_faces_calc_data({0})'.format(test))
    g = gsu.Grid()
    g.load('cases/grids/{0}.dat'.format(test))
    g.store_faces_calc_data('cases/grids/{0}.txt'.format(test))
    g.load_faces_calc_data('cases/grids/{0}.txt'.format(test))
    g.store('cases/grids/{0}_data.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_009_store_mpi(test='bunny'):
    """
    Load grid, decompose it and store in MPI.
    Test objective:

      Check we can store data for multiprocessing mode calculations.

    :param test: test name
    """

    print('case_009_store_mpi({0})'.format(test))
    g = gsu.Grid()
    g.load('cases/grids/{0}.dat'.format(test))
    g.decompose_hierarchical(levels=3, new_name=test + ' hierarchical')
    g.store_mpi('cases/grids/{0}_mpi'.format(test), '000000000100')


# --------------------------------------------------------------------------------------------------


def case_014_convert_grid_stall_to_check_point():
    """
    Load grid, convert and store.
    Test objective:

      Verification grids conversion.
    """

    print('case_014_convert_grid_stall_to_check_point')
    g = gsu.Grid()
    g.load('cases/drops/cyl/cyl_stall.dat')
    g.convert_grid_stall_to_check_point()
    g.store('cases/drops/cyl/cyl_stall_case_014.dat')

# --------------------------------------------------------------------------------------------------


def case_015_GloId_in_grid_for_divide_face():

    g = gsu.Grid()
    g.load('cases/grids/wing_1.dat')
    print('face id')
    print([f.GloId for f in g.Faces])
    print('edge id')
    print([f.GloId for f in g.Edges])
    print('node id')
    print([f.GloId for f in g.Nodes])

    f = g.Faces[0]
    g.divide_face(f, f.get_triangle().centroid())

    print('face id after')
    print([f.GloId for f in g.Faces])
    print('edge id after')
    print([f.GloId for f in g.Edges])
    print('node id after')
    print([f.GloId for f in g.Nodes])

# --------------------------------------------------------------------------------------------------


def case_016_GloId_in_grid_for_collapse_face():

    g = gsu.Grid()
    g.load('cases/grids/wing_1.dat')
    print('face id')
    print([f.GloId for f in g.Faces])
    print('edge id')
    print([f.GloId for f in g.Edges])
    print('node id')
    print([f.GloId for f in g.Nodes])

    f = g.Faces[0]
    g.collapse_face(f)

    print('face id after')
    print([f.GloId for f in g.Faces])
    print('edge id after')
    print([f.GloId for f in g.Edges])
    print('node id after')
    print([f.GloId for f in g.Nodes])

# --------------------------------------------------------------------------------------------------


def case_017_GloId_in_grid_for_cut_edge():

    g = gsu.Grid()
    g.load('cases/grids/wing_1.dat')
    print('face id')
    print([f.GloId for f in g.Faces])
    print('edge id')
    print([f.GloId for f in g.Edges])
    print('node id')
    print([f.GloId for f in g.Nodes])

    n = 0
    len(g.Faces[n].Edges[0].Faces)
    w = True
    while w:
        w = not len(g.Faces[n].Edges[0].Faces) == 2
        if not w:
            g.cut_edge(g.Faces[n].Edges[0], g.Faces[n].get_triangle().centroid())
        n += 1

    print('face id after')
    print([f.GloId for f in g.Faces])
    print('edge id after')
    print([f.GloId for f in g.Edges])
    print('node id after')
    print([f.GloId for f in g.Nodes])

# --------------------------------------------------------------------------------------------------


def case_018_GloId_in_grid_for_collapse_edge():

    g = gsu.Grid()
    g.load('cases/grids/wing_1.dat')
    print('face id')
    print([f.GloId for f in g.Faces])
    print('edge id')
    print([f.GloId for f in g.Edges])
    print('node id')
    print([f.GloId for f in g.Nodes])

    n = 0
    len(g.Faces[n].Edges[0].Faces)
    w = True
    while w:
        w = not len(g.Faces[n].Edges[0].Faces) == 2
        if not w:
            g.collapse_edge(g.Faces[n].Edges[0])
        n += 1

    print('face id after')
    print([f.GloId for f in g.Faces])
    print('edge id after')
    print([f.GloId for f in g.Edges])
    print('node id after')
    print([f.GloId for f in g.Nodes])

# --------------------------------------------------------------------------------------------------


def case_019_GloId_in_grid_for_cut_single_edge():

    g = gsu.Grid()
    g.load('cases/grids/wing_1.dat')
    print('face id')
    print([f.GloId for f in g.Faces])
    print('edge id')
    print([f.GloId for f in g.Edges])
    print('node id')
    print([f.GloId for f in g.Nodes])

    n = 0
    len(g.Faces[n].Edges[0].Faces)
    w = True
    while w:
        w = not len(g.Faces[n].Edges[0].Faces) == 1
        if not w:
            g.cut_single_edge(g.Faces[n].Edges[0], g.Faces[n].get_triangle().centroid())
        n += 1

    print('face id after')
    print([f.GloId for f in g.Faces])
    print('edge id after')
    print([f.GloId for f in g.Edges])
    print('node id after')
    print([f.GloId for f in g.Nodes])

# --------------------------------------------------------------------------------------------------


def case_020_GloId_in_grid_for_cut_edge_with_two_nodes():

    g = gsu.Grid()
    g.load('cases/grids/wing_1.dat')
    print('face id')
    print([f.GloId for f in g.Faces])
    print('edge id')
    print([f.GloId for f in g.Edges])
    print('node id')
    print([f.GloId for f in g.Nodes])

    n = 0
    len(g.Faces[n].Edges[0].Faces)
    w = True
    while w:
        w = not len(g.Faces[n].Edges[0].Faces) == 1
        if not w:
            g.cut_edge_with_two_nodes(g.Faces[n].Edges[0],
                                      g.Faces[n].get_triangle().centroid() - geom.Vect(0.01, 0.01, 0.01),
                                      g.Faces[n].get_triangle().centroid() + geom.Vect(0.01, 0.01, 0.01))
        n += 1

    print('face id after')
    print([f.GloId for f in g.Faces])
    print('edge id after')
    print([f.GloId for f in g.Edges])
    print('node id after')
    print([f.GloId for f in g.Nodes])

    n = 0
    len(g.Faces[n].Edges[0].Faces)
    w = True
    while w:
        w = not len(g.Faces[n].Edges[0].Faces) == 2
        if not w:
            g.cut_edge_with_two_nodes(g.Faces[n].Edges[0],
                                      g.Faces[n].get_triangle().centroid() - geom.Vect(0.01, 0.01, 0.01),
                                      g.Faces[n].get_triangle().centroid() + geom.Vect(0.01, 0.01, 0.01))
        n += 1

    print('face id after')
    print([f.GloId for f in g.Faces])
    print('edge id after')
    print([f.GloId for f in g.Edges])
    print('node id after')
    print([f.GloId for f in g.Nodes])

# ==================================================================================================


if __name__ == '__main__':
    # case_001_load_store()
    # case_002_decompose()
    # case_005_explode_bunny()
    # case_007_store_load_faces_calc_data()
    # case_009_store_mpi()
    # case_014_convert_grid_stall_to_check_point()
    # case_015_GloId_in_grid_for_divide_face()
    # case_016_GloId_in_grid_for_collapse_face()
    # case_017_GloId_in_grid_for_cut_edge()
    # case_018_GloId_in_grid_for_collapse_edge()
    # case_019_GloId_in_grid_for_cut_single_edge()
    # case_020_GloId_in_grid_for_cut_edge_with_two_nodes()

    pass

# ==================================================================================================
