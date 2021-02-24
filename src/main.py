"""
Main module.
"""

import gsu

# ----------------------------------------------------------------------------------------------------------------------


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

# ----------------------------------------------------------------------------------------------------------------------


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
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(), gsu.fun_face_cy(), gsu.fun_face_cz()],
                             levels=6,
                             new_name=test + ' hierarchical')
    g.store('grids/{0}_hierarchical.dat'.format(test))

# ----------------------------------------------------------------------------------------------------------------------


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
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(), gsu.fun_face_cy(), gsu.fun_face_cz()],
                             levels=6,
                             new_name=test + ' hierarchical',
                             fixed_zones=['POS1', 'POS2'])
    g.store('grids/{0}_hierarchical.dat'.format(test))

# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    # case_001_node_face_data()
    # case_002_decompose()
    case_003_fixed_zones()

# ----------------------------------------------------------------------------------------------------------------------
