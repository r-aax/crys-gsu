"""
Main module.
"""

import gsu

# ----------------------------------------------------------------------------------------------------------------------


def case_001_node_face_data():
    """
    Load small grid and store in without any changes.
    Test objective:
      To make sure that node data and face data are correct.
      Node data is just a point.
    """

    print('case_001_node_face_data:\n')
    test = 'wing_1_mz'
    g = gsu.Grid()
    g.load('grids/{0}.dat'.format(test))
    g.store('grids/{0}_original.dat'.format(test))

# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    case_001_node_face_data()

# ----------------------------------------------------------------------------------------------------------------------
