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
    grid_zones_names = [z.Name for z in g.Zones]
    for fz in fixed_zones:
        if not fz in grid_zones_names:
            raise Exception('crys-gsu-split : no zone with name "{0}" in the grid'.format(fz))

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


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='split',
                                     description='GSU split script into separate domains.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('grid_file', help='grid file name')
    parser.add_argument('json_file', help='JSON file name')
    parser.add_argument('cry_dir', help='out dir for *.cry files')
    args = parser.parse_args()

    # Extract parameters.
    strategy = extract_strategy_from_json(args.json_file)
    fixed_zones = extract_fixed_zones_from_json(args.json_file)

    split(grid_file=args.grid_file,
          cry_dir=args.cry_dir,
          split_strategy=strategy,
          fixed_zones=fixed_zones)

# --------------------------------------------------------------------------------------------------
