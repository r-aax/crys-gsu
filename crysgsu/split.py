"""
Split grid for MPI.
"""

import pathlib
import utils
import sys
import os
import json
import time
import gsu

# --------------------------------------------------------------------------------------------------


def split(grid_file,
          cry_dir, out_dat_file,
          split_strategy, fixed_zones=[], recovery_factor=1.0):
    """
    Split grid.

    :param grid_file: file with grid
    :param cry_dir: directory for \*.cry files
    :param out_dat_file: out \*.dat file
    :param split_strategy: strategy for grid splitting
    :param fixed_zones: list of fixed zones
    :return: actual count of zones after splitting
    """

    start_time = time.time()

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
          '{0} (bs={1}, nm={2}, ts={3}, sf={4}), cry-dir={5}, '
          'out-grid-file={6}, split-strategy={7}, '
          'fixed-zones={8}'.format(grid_file, bs, nm, ts, sf,
                                   cry_dir, out_dat_file, split_strategy, fixed_zones))

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

    # Issue#12 recovery_factor support
    # Check recovery factor in grid faces data.
    # If there is no recovery factor in face, then add it.
    recovery_factor_name = 'RecoveryFactor'
    for f in g.Faces:
        if not recovery_factor_name in f.Data.keys():
            f.Data[recovery_factor_name] = recovery_factor

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
        g.decompose_hierarchical(levels=n + 1,
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
        g.decompose_farhat(count=n - fz, fz_names=fixed_zones)
        actual_zones_count = len(g.Zones)
        print('crys-gsu-split : {0} zones : done'.format(actual_zones_count))
    else:
        raise Exception('crys-gsu-split : unknown split strategy ({0})'.format(split_strategy))

    # Store for MPI.
    if not out_dat_file is None:
        g.store(out_dat_file)
    elif not cry_dir is None:
        try:
            os.makedirs(cry_dir)
        except FileExistsError:
            # If directory already exists - there is nothing to do.
            pass
        g.store_mpi('{0}/{1}'.format(cry_dir, nm), ts)
    else:
        raise Exception('crys-gsu-split : internal error')


    time_stat = 'crys-gsu-split : done, ' \
                '{0:.3f} seconds spent)'.format(time.time() - start_time)
    print(time_stat)
    if not cry_dir is None:
        with open(cry_dir + '/report.txt', 'w') as f:
            print(time_stat, file=f)

    return g

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


def extract_recovery_factor_from_json(filename):
    """
    Extract recovery factor from json.

    :param filename: file name
    :return: recovery factor or 1.0
    """

    if not os.path.isfile(filename):
        raise Exception('crys-gsu-split : no such file ({0})'.format(filename))

    with open(filename, 'r') as jf:
        data = json.load(jf)
        jf.close()

    try:
        return data['recovery_factor']
    except KeyError:
        return 1.0

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='split',
                                     description='GSU split script into separate domains.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('grid_file', help='grid file name')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--cry_dir', help='result in *.cry format (directory)')
    group.add_argument('-o', '--out_dat_file', help='result in *.dat format (file)')
    parser.add_argument('-j', '--json_file', help='JSON file name')
    parser.add_argument('-s', '--strategy', help='strategy of decomposition (z, n<n> or h<h>)')
    parser.add_argument('-v', '--verbosity', action='count',
                        help='increase output verbosity', default=0)
    args = parser.parse_args()

    # Extract split_strategy parameter.
    if not args.strategy is None:
        split_strategy = args.strategy
    elif not args.json_file is None:
        split_strategy = extract_strategy_from_json(args.json_file)
    else:
        split_strategy = 'n1'

    # Extract fixed zones parameters.
    if not args.json_file is None:
        fixed_zones = extract_fixed_zones_from_json(args.json_file)
    else:
        fixed_zones = []

    # Extract recovery_factor from json.
    if not args.json_file is None:
        recovery_factor = extract_recovery_factor_from_json(args.json_file)
    else:
        recovery_factor = 1.0

    split(grid_file=args.grid_file,
          cry_dir=args.cry_dir,
          out_dat_file=args.out_dat_file,
          split_strategy=split_strategy,
          fixed_zones=fixed_zones,
          recovery_factor=recovery_factor)

# --------------------------------------------------------------------------------------------------
