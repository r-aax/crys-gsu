"""
Split grid for MPI.
"""

import pathlib
import gsu
import utils
import sys
import os

# ----------------------------------------------------------------------------------------------------------------------


def split(grid_file, cry_dir, split_policy, fixed_zones=[]):
    """
    Split grid.
    :param grid_file: file with grid
    :param cry_dir: directory for *.cry files
    :param split_policy: policy for grid splitting
    :param fixed_zones: list of fixed zones
    :return: actual count of zones after splitting
    """

    pp = pathlib.PurePath(grid_file)
    # Get characteristics of file:
    #   bs - base of path
    #   nm - name without timestamp
    #   ts - timestamp string
    #   sf - suffix
    bs, sf = str(format(pp.parents[0])), pp.suffix
    nm, ts = utils.get_filename_and_timestamp_pair(pp.stem)

    print('crys-gsu-split : grid file='
          '{0} (bs={1}, nm={2}, ts={3}, sf={4})'.format(grid_file, bs, nm, ts, sf))

    # Check file name.
    if sf != '.dat':
        raise Exception('Wrong name of grid file (extension must be *.dat)')

    # Load grid.
    g = gsu.Grid()
    g.load(grid_file)

    # Decompose grid.
    if split_policy[0] == 'h':
        # Hierarchical split.
        n = int(split_policy[1:])
        target_zones_count = 2 ** n + len(fixed_zones)
        print('crys-gsu-split : hierarchical, target_zones_count={0} '
              '(2^{1} + {2})'.format(target_zones_count, n, len(fixed_zones)))
        g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(), gsu.fun_face_cy(), gsu.fun_face_cz()],
                                 levels=n + 1,
                                 fixed_zones=fixed_zones)
        actual_zones_count = len(g.Zones)
        print('crys-gsu-split : hierarchical, actual_zones_count={0}'.format(actual_zones_count))
    else:
        raise Exception('unknown split policy')

    # Store for MPI.
    try:
        os.makedirs(cry_dir)
    except FileExistsError:
        # If directory already exists - there is nothing to do.
        pass
    g.store_mpi('{0}/{1}'.format(cry_dir, nm), ts)

    print('crys-gsu-split : done')

    return actual_zones_count

# ----------------------------------------------------------------------------------------------------------------------


# split.py should be called from shell script in the following manner:
#     split.py grids/bunny.dat grids/cry h2 "POS1" "POS2"
if __name__ == '__main__':

    if sys.argv[1] == '-h':
        print('split.py <grid-file> <output-dir> <split-strategy> <ais-zone-name-1> ... <ais-zone-name-n>')
        print('    <split-strategy>:')
        print('        h<n> - hierarchical split into 2^n zones in addition to AIS zones')
        exit(0)

    split(grid_file=sys.argv[1],
          cry_dir=sys.argv[2],
          split_policy=sys.argv[3],
          fixed_zones=sys.argv[4:])

# ----------------------------------------------------------------------------------------------------------------------
