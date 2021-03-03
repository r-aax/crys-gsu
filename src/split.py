"""
Split grid for MPI.
"""

import pathlib
import gsu
import utils
import sys

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

    pp = pathlib.PurePath(grid_file)
    # Get characteristics of file:
    #   bs - base of path
    #   nm - name without timestamp
    #   ts - timestamp string
    #   sf - suffix
    bs, sf = str(format(pp.parents[0])), pp.suffix
    nm, ts = utils.get_filename_and_timestamp_pair(pp.stem)

    print('split_2_power_n : grid file = '
          '{0} (bs {1}, nm {2}, ts {3}, sf {4})'.format(grid_file, bs, nm, ts, sf))
    print('split_2_power_n : total zones count = '
          '{0} (2^{1} + {2} fixed zones)'.format(2 ** n + len(fixed_zones), n, len(fixed_zones)))

    # Check file name.
    if sf != '.dat':
        raise Exception('Wrong name of grid file (extension must be *.dat)')

    # Load grid.
    g = gsu.Grid()
    g.load(grid_file)

    # Decompose grid.
    g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(), gsu.fun_face_cy(), gsu.fun_face_cz()],
                             levels=n + 1,
                             fixed_zones=fixed_zones)

    # Store for MPI.
    g.store_mpi('{0}/{1}'.format(bs, nm), ts)

    print('split_2_power_n : done')

# ----------------------------------------------------------------------------------------------------------------------


# split.py should be called from shell script in the following manner:
#     py "grids/bunny.dat" 2 "POS1" "POS2"
if __name__ == '__main__':

    split_2_power_n(grid_file=sys.argv[1],
                    n=int(sys.argv[2]),
                    fixed_zones=sys.argv[3:])

# ----------------------------------------------------------------------------------------------------------------------
