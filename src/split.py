"""
Split grid for MPI.
"""

import pathlib
import gsu
import utils
import sys
import os
import json

# --------------------------------------------------------------------------------------------------


def split(grid_file, cry_dir, split_strategy, fixed_zones=[]):
    """
    Split grid.
    :param grid_file: file with grid
    :param cry_dir: directory for *.cry files
    :param split_strategy: strategy for grid splitting
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
    if split_strategy[0] == 'h':
        # Hierarchical split.
        n = int(split_strategy[1:])
        target_zones_count = 2 ** n + len(fixed_zones)
        print('crys-gsu-split : hierarchical, target-zones-count={0} '
              '(2^{1} + {2}) // {3}'.format(target_zones_count,
                                            n,
                                            len(fixed_zones),
                                            fixed_zones))
        g.decompose_hierarchical(extract_signs_funs=[gsu.fun_face_cx(),
                                                     gsu.fun_face_cy(),
                                                     gsu.fun_face_cz()],
                                 levels=n + 1,
                                 fixed_zones=fixed_zones)
        actual_zones_count = len(g.Zones)
        print('crys-gsu-split : hierarchical, actual-zones-count={0}'.format(actual_zones_count))
    elif split_strategy[0] == 'n':
        # Split with exact n domains defined.
        n = int(split_strategy[1:])
        fz = len(fixed_zones)
        if n == 0:
            raise Exception('crys-gsu-split : zero zones count to split')
        if fz >= n:
            raise Exception('crys-gsu-split : fixed zones count is more or equal '
                            'than target zones count ({0} zones, {1} fixed zones)'.format(n, fz))
        print('crys-gsu-split : {0} zones ({1} fixed zones) : begin'.format(n, fz))
        g.decompose_pressure(count=n - fz, fz_names=fixed_zones)
        actual_zones_count = len(g.Zones)
        print('crys-gsu-split : {0} zones : done'.format(actual_zones_count))
    else:
        raise Exception('crys-gsu-split : unknown split strategy ({0})'.format(split_strategy))

    # Store for MPI.
    try:
        os.makedirs(cry_dir)
    except FileExistsError:
        # If directory already exists - there is nothing to do.
        pass
    g.store_mpi('{0}/{1}'.format(cry_dir, nm), ts)

    print('crys-gsu-split : done')

    return actual_zones_count

# --------------------------------------------------------------------------------------------------


def print_help():
    """
    Print help.
    """

    print('[Overview]:')
    print('split.py script reads *.dat file with surface grid and '
          'splits it into several *.cry files.')
    print('All additional parameters are contained in *.json file.')
    print('')
    print('[Usage]:')
    print('split.py <grid-file> <output-dir> <json-file>')
    print('    <grid-file> - full name of grid file')
    print('    <json-file> - json file with parameters')
    print('    <output-dir> - directory name for output files')
    print('        if there is no such directory, it will be created')
    print('')
    print('[Examples]:')
    print('bunny.dat -> bunny_<mpi_i>_000000000000.cry')
    print('bunny_000000000100.dat -> bunny_<mpi_i>_000000000100.cry')

# --------------------------------------------------------------------------------------------------


def extract_strategy_from_json(filename):
    """
    Extract strategy from json file.
    :param filename: file name
    :return: strategy
    """

    if not os.path.isfile(filename):
        raise Exception('crys-gsu-split : no such file ({0})'.format(filename))

    with open(filename, 'r') as jf:
        data = json.load(jf)
        jf.close()

    try:
        return data['split_strategy']
    except KeyError:
        return 'n1'


# --------------------------------------------------------------------------------------------------


def extract_fixed_zones_from_json(filename):
    """
    Extract fixed zones from json.
    :param filename: file name
    :return: list of zones' names
    """

    if not os.path.isfile(filename):
        raise Exception('crys-gsu-split : no such file ({0})'.format(filename))

    with open(filename, 'r') as jf:
        data = json.load(jf)
        jf.close()

    try:
        els = data['heating_system']['elements']
        return [e['zone'] for e in els]
    except KeyError:
        return []

# --------------------------------------------------------------------------------------------------


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

    if c != 4:
        print('crys-gsu-split : script receives exactly 3 parameters')
        print_help()
        exit(0)

    # Extract parameters.
    strategy = extract_strategy_from_json(sys.argv[2])
    fixed_zones = extract_fixed_zones_from_json(sys.argv[2])

    split(grid_file=sys.argv[1],
          cry_dir=sys.argv[3],
          split_strategy=strategy,
          fixed_zones=fixed_zones)

# --------------------------------------------------------------------------------------------------
