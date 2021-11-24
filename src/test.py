"""
Test module.
"""

import gsu
import gsu_geom
import split


# --------------------------------------------------------------------------------------------------


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
    g.load('grids/{0}.dat'.format(test))

    # Create new name of grid and file
    new_grid_name = lambda m: '{0} {1}'.format(test, m)
    new_file_name = lambda m: 'grids/{0}_{1}.dat'.format(test, m)

    if 'mono' in methods:
        g.decompose_mono(new_name=new_grid_name('mono'))
        g.store(new_file_name('mono'))

    if 'random' in methods:
        g.decompose_random(new_name=new_grid_name('random'))
        g.store(new_file_name('random'))

    if 'linear' in methods:
        g.decompose_linear(new_name=new_grid_name('linear'))
        g.store(new_file_name('linear'))

    if 'hierarchical' in methods:
        g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                     gsu.fun_face_cy(),
                                                     gsu.fun_face_cz()],
                                 new_name=new_grid_name('hierarchical'),
                                 fixed_zones=['POS1', 'POS2'])
        g.store(new_file_name('hierarchical'))

    if 'farhat' in methods:
        g.decompose_farhat(new_name=new_grid_name('farhat'),
                           fz_names=['POS1', 'POS2'])
        g.store(new_file_name('farhat'))


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


def case_010_add_mimp2_vd2(test='bunny'):
    """
    Load grid without MImp2, Vd2 fields and save it with MImp2, Vd2 fields (all zeros).
    Test objective:
      Check that program supports work with grids of both formats.
    :param test: test name
    """

    print('case_010_add_mimp2_vd2({0})'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.store('grids/{0}_mimp2_vd2.dat'.format(test))


# --------------------------------------------------------------------------------------------------


def case_013_clean_grid(test='cyl'):
    """
    Load grid, clean and store.
    Test objective:
      Verification of cleaning mechanism.
    :param test: test name
    """

    print('case_013_clean_grid({0})'.format(test))
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    for f in g.Faces:
        f.set_t(0.0)
        f.set_hw(0.0)
        f.set_hi(0.0)
        f.set_mimp2(0.0)
        f.set_vd2(0.0)
    g.store('grids/{0}_case_013.dat'.format(test))

# --------------------------------------------------------------------------------------------------


def case_014_convert_grid_stall_to_check_point(test='cyl/cyl_1_stall'):
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

# --------------------------------------------------------------------------------------------------


def case_015_self_intersection():
    """
    Self intersection of the grid.
    Test objective:
      Self intersection algorithm.
    """

    m = gsu_geom.Mesh()
    m.load('meshes/cyl_100.dat')
    m.shred()
    m.store('meshes/cyl_100_sh.dat')

# --------------------------------------------------------------------------------------------------


def case_016_wrapping():
    """
    Wrapping mechanism test.
    Test objective:
      Wrapping algorithm.
    """

    m = gsu_geom.Mesh()
    m.load('meshes/cyl_100_sh.dat')
    m.wrap()
    m.store('meshes/cyl_100_sh_wr.dat')

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    # case_001_load_store()
    # case_002_decompose()
    # case_005_explode_bunny()
    # case_006_store_faces_t_hw_hi()
    # case_007_load_faces_t_hw_hi()
    # case_009_store_mpi()
    # case_010_add_mimp2_vd2()
    # case_013_clean_grid()
    # case_014_convert_grid_stall_to_check_point()
    # case_015_self_intersection()
    # case_016_wrapping()

    pass

# --------------------------------------------------------------------------------------------------
