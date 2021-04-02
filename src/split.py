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

    # Check for grid file.
    if not os.path.isfile(grid_file):
        raise Exception('crys-gsu-split : no such file ({0})'.format(grid_file))

    # Start splitting.
    print('crys-gsu-split : grid-file='
          '{0} (bs={1}, nm={2}, ts={3}, sf={4})'.format(grid_file, bs, nm, ts, sf))

    # Check file name.
    if sf != '.dat':
        raise Exception('crys-gsu-split : wrong name of grid file (extension must be *.dat)')

    # Load grid.
    g = gsu.Grid()
    g.load(grid_file)

    # Decompose grid.
    if split_policy[0] == 'h':
        # Hierarchical split.
        n = int(split_policy[1:])
        target_zones_count = 2 ** n + len(fixed_zones)
        print('crys-gsu-split : hierarchical, target-zones-count={0} '
              '(2^{1} + {2}) // {3}'.format(target_zones_count,
                                            n,
                                            len(fixed_zones),
                                            fixed_zones))
        g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(), gsu.fun_face_cy(), gsu.fun_face_cz()],
                                 levels=n + 1,
                                 fixed_zones=fixed_zones)
        actual_zones_count = len(g.Zones)
        print('crys-gsu-split : hierarchical, actual-zones-count={0}'.format(actual_zones_count))
    else:
        raise Exception('crys-gsu-split : unknown split policy ({0})'.format(split_policy))

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


def print_help():
    """
    Print help.
    """

    print('[Overview]:')
    print('split.py script reads *.dat file with surface grid and split it into several *.cry files')
    print('')
    print('[Usage]:')
    print('split.py <grid-file> <output-dir> <split-strategy> <ais-zone-name-1> ... <ais-zone-name-n>')
    print('    <grid-file> - full name of grid file')
    print('    <output-dir> - directory name for output files')
    print('        if there is no such directory, it will be created')
    print('    <split-strategy>:')
    print('        h<n> - hierarchical split into 2^n zones in addition to AIS zones')
    print('    <json> - json file from which names of ones should be extracted')
    print('')
    print('[Examples]:')
    print('bunny.dat -> bunny_<mpi_i>_000000000000.cry')
    print('bunny_000000000100.dat -> bunny_<mpi_i>_000000000100.cry')


# ----------------------------------------------------------------------------------------------------------------------

def extract_fixed_zones_from_json(filename):
    """
    Extract fixed zones from json.
    :param filename: file name
    :return:
    """

    f = open(filename, 'r')
    ll = f.readlines()
    f.close()

    # Find zones for AIS.
    ll = [s.split(':')[1] for s in ll if '"zone"' in s]
    ll = [s[s.index('"') + 1: s.rindex('"')] for s in ll]

    return ll

# ----------------------------------------------------------------------------------------------------------------------


# split.py should be called from shell script in the following manner:
#     split.py grids/bunny_pos.dat grids/cry h2 grids/bunny.json
if __name__ == '__main__':

    c = len(sys.argv)

    if c == 1:
        print_help()
        exit(0)

    if (sys.argv[1] == '-h') or (sys.argv[1] == '--help'):
        print_help()
        exit(0)

    if c < 4:
        raise Exception('crys-gsu-split : not enough arguments')
    elif c == 4:
        fixed_zones = []
    elif c == 5:
        fixed_zones = extract_fixed_zones_from_json(sys.argv[4])
    else:
        raise Exception('crys-gsu-split : too many arguments')

    split(grid_file=sys.argv[1],
          cry_dir=sys.argv[2],
          split_policy=sys.argv[3],
          fixed_zones=fixed_zones)

# ----------------------------------------------------------------------------------------------------------------------
