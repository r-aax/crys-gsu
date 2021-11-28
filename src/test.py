"""
Test module.
"""

from gsu import gsu
import gsu_geom
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
    g.load('grids/{0}.dat'.format(test))
    g.store('grids/{0}_store1.dat'.format(test))
    g.load('grids/{0}_store1.dat'.format(test))
    g.store('grids/{0}_store2.dat'.format(test))

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
    new_file_name = lambda m: 'grids/{0}_{1}.dat'.format(test, m)

    if 'mono' in methods:
        g.load('grids/{0}.dat'.format(test))
        g.decompose_mono(new_name=new_grid_name('mono'))
        g.store(new_file_name('mono'))

    if 'random' in methods:
        g.load('grids/{0}.dat'.format(test))
        g.decompose_random(new_name=new_grid_name('random'))
        g.store(new_file_name('random'))

    if 'linear' in methods:
        g.load('grids/{0}.dat'.format(test))
        g.decompose_linear(new_name=new_grid_name('linear'))
        g.store(new_file_name('linear'))

    if 'hierarchical' in methods:
        g.load('grids/{0}.dat'.format(test))
        g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                     gsu.fun_face_cy(),
                                                     gsu.fun_face_cz()],
                                 new_name=new_grid_name('hierarchical'),
                                 fixed_zones=['POS1', 'POS2'])
        g.store(new_file_name('hierarchical'))

    if 'farhat' in methods:
        g.load('grids/{0}.dat'.format(test))
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
    g.load('grids/{0}.dat'.format(test))
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                 gsu.fun_face_cy(),
                                                 gsu.fun_face_cz()],
                             levels=5,
                             new_name=test + ' hierarchical')
    g.store('grids/{0}_hierarchical.dat'.format(test))
    g.load('grids/{0}_hierarchical.dat'.format(test), is_merge_same_nodes=False)
    g.move_from_mean_point(0.25)
    g.store('grids/{0}_explode.dat'.format(test))


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
    g.load('grids/{0}.dat'.format(test))
    g.store_faces_calc_data('grids/{0}.txt'.format(test))
    g.load_faces_calc_data('grids/{0}.txt'.format(test))
    g.store('grids/{0}_data.dat'.format(test))


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
    g.load('grids/{0}.dat'.format(test))
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                 gsu.fun_face_cy(),
                                                 gsu.fun_face_cz()],
                             levels=3,
                             new_name=test + ' hierarchical')
    g.store_mpi('grids/{0}_mpi'.format(test), '000000000100')


# --------------------------------------------------------------------------------------------------


def case_014_convert_grid_stall_to_check_point(test='cyl/cyl_stall'):
    """
    Load grid, convert and store.
    Test objective:
      Verification grids conversion.
    :param test: test name
    """

    print('case_014_convert_grid_stall_to_check_point({0})'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.convert_grid_stall_to_check_point()
    g.store('grids/{0}_case_014.dat'.format(test))


# ==================================================================================================


if __name__ == '__main__':
    # case_001_load_store()
    # case_002_decompose()
    # case_005_explode_bunny()
    # case_007_store_load_faces_calc_data()
    # case_009_store_mpi()
    # case_014_convert_grid_stall_to_check_point()

    pass

# ==================================================================================================
