"""
Split grid for MPI.
"""

import pathlib
import gsu

# ----------------------------------------------------------------------------------------------------------------------


def split_2_power_n(grid_file, n, fixed_zones=[]):
    """
    Split grid on 2^n pieces (in additional to fixed zones).
    So if there are m fixed zones, total count of zones is 2^n + m.
    Splitting is produced by hierarchical method.
    :param grid_file: file with grid
    :param n: power of 2
    :param fixed_zones: list of fixed zones
    :return:
    """

    print('split_2_power_n : grid file = {0}'.format(grid_file))
    print('split_2_power_n : total zones count = '
          '{0} (2^{1} + {2} fixed zones)'.format(2 ** n + len(fixed_zones), n, len(fixed_zones)))

    pp = pathlib.PurePath(grid_file)
    base, suff = '{0}/{1}'.format(pp.parents[0], pp.stem), pp.suffix

    # Check file name.
    if suff != '.dat':
        raise Exception('Wrong name of grid file (extension must be *.dat)')

    # Load grid.
    g = gsu.Grid()
    g.load(grid_file)

    # Decompose grid.
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(), gsu.fun_face_cy(), gsu.fun_face_cz()],
                             levels=n + 1,
                             fixed_zones=fixed_zones)

    # Store for MPI.
    g.store_mpi(base)


# ----------------------------------------------------------------------------------------------------------------------
