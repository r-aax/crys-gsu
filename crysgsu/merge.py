"""
Merge \*.txt files with original grid.
"""

import sys
import os
import pathlib
import utils
from gsu import gsu
import time
import numpy as np
from numpy import genfromtxt

# --------------------------------------------------------------------------------------------------


def merge(grid_file, src_dir, dst_dir):
    """
    Merge txt files with original grid.
    Merge files separately for each timestamp and produce out _r_ file for crys-remesh tool.

    :param grid_file: file with original grid
    :param src_dir: directory with sources
    :param dst_dir: destination directory
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
        raise Exception('crys-gsu-merge : no such file ({0})'.format(grid_file))

    print('crys-gsu-merge : grid-file={0} (bs={1}, nm={2}, ts={3}, sf={4}), '
          'src_dir={5}, dst_dir={6}'.format(grid_file, bs, nm, ts, sf,
                                            src_dir, dst_dir))

    all_files = os.listdir(src_dir)

    # We have to analyze only our *.txt files.
    match_files = [fn for fn in all_files if utils.is_filename_correct_crys_txt_file(fn, nm)]

    # Group match files.
    gmf = utils.group_txt_files_by_timestamps(match_files)

    # Make dirs
    try:
        os.makedirs(dst_dir)
    except FileExistsError:
        # Ignore if directories already exist.
        pass

    # Add all faces data from cry files.
    g = gsu.Grid()
    for tm in gmf:
        fs = gmf.get(tm)
        g.load(grid_file)
        for f in fs:
            g.load_faces_calc_data(src_dir + '/' + f)
        g.store('{0}/{1}_r_{2}.dat'.format(dst_dir, nm, tm))

    # Function for reading 2D array from *.csv file.
    # #172 : Make 2D data matrix manually if csv file contains just one row.
    def read_csv2d(file_name):
        data = genfromtxt(file_name, delimiter=';', skip_header=True)
        if len(data.shape) == 1:
            data = np.array([data])
        return data

    # Finally put total *.csv file to dst_dir/../ directory.
    csv_files = [fn for fn in all_files if fn[-4:] == '.csv']
    if len(csv_files) > 0:
        d = read_csv2d(src_dir + '/' + csv_files[0])
        t = d[:, :1]
        d = d[:, 1:]
        for csv_file in csv_files[1:]:
            d = d + read_csv2d(src_dir + '/' + csv_file)[:, 1:]
        d = np.hstack((t, d))
        np.savetxt('{0}/../{1}.csv'.format(dst_dir, nm), d, delimiter=";",
                   header='time;M_imp;M_es;M_stall;M_ini_water;M_water;M_ini_ice;M_ice;Q_conv;Q_es',
                   comments='')

    time_stat = 'crys-gsu-merge : done ({0} _r_ files generated, ' \
                '{1:.3f} seconds spent)'.format(len(gmf.keys()), time.time() - start_time)
    print(time_stat)
    with open(dst_dir + '/report.txt', 'w') as f:
        print(time_stat, file=f)

# --------------------------------------------------------------------------------------------------


# Example of running merge.py script:
#     merge.py grids/bunny.dat dir_to_txt_files
# Script produces out files for crys-remesh tool.
if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='merge',
                                     description='Merge calc data.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('grid_file', help='grid file name')
    parser.add_argument('src_dir', help='directory with source files')
    parser.add_argument('dst_dir', help='destination directory')
    parser.add_argument('-v', '--verbosity', action='count',
                        help='increase output verbosity', default=0)
    args = parser.parse_args()

    merge(grid_file=args.grid_file,
          src_dir=args.src_dir,
          dst_dir=args.dst_dir)

# --------------------------------------------------------------------------------------------------
