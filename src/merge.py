"""
Merge *.txt files with original grid.
"""

import sys
import os
import pathlib
import utils
import gsu
import time

# --------------------------------------------------------------------------------------------------


def merge(grid_file, txt_files_dir, r_files_dir):
    """
    Merge txt files with original grid.
    Merge files separately for each timestamp and produce out _r_ file for crys-remesh tool.
    :param grid_file: file with original grid
    :param txt_files_dir: dir with txt files
    :param r_files_dir: out dir for _r_ files
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
          'txt_files_dir={5}, r_files_dir={6}'.format(grid_file, bs, nm, ts, sf,
                                                      txt_files_dir, r_files_dir))

    all_files = os.listdir(txt_files_dir)

    # We have to analyze only our *.txt files.
    match_files = [fn for fn in all_files if utils.is_filename_correct_crys_txt_file(fn, nm)]

    # Group match files.
    gmf = utils.group_txt_files_by_timestamps(match_files)

    # Make dirs
    try:
        os.makedirs(r_files_dir)
    except FileExistsError:
        # Ignore if directories already exist.
        pass

    # Add all faces data from cry files.
    g = gsu.Grid()
    for tm in gmf:
        fs = gmf.get(tm)
        g.load(grid_file)
        for f in fs:
            g.load_faces_calc_data(txt_files_dir + '/' + f)
        g.store('{0}/{1}_r_{2}.dat'.format(r_files_dir, nm, tm))

    print('crys-gsu-merge : done ({0} _r_ files generated, '
          '{1:.3f} seconds spent)'.format(len(gmf.keys()), time.time() - start_time))

# --------------------------------------------------------------------------------------------------


def print_help():
    """
    Print help.
    """

    print('[Overview]:')
    print('merge.py script merge initial *.dat grid file with *.txt data files')
    print('         and produces *.dat files with timestamps')
    print('')
    print('[Usage]:')
    print('merge.py <grid-file> <txt-files-dir> <r-files-dir>')
    print('    <grid-file> - initial file name')
    print('    <txt-files-dir> - name of directory with data files')
    print('    <r-files-dir> - name of directory with output *_r_* files')
    print('')
    print('[Examples]:')
    print('            [bunny_00000_000000000100.txt bunny_00001_000000000100.txt] '
          '-> [bunny_r_000000000100.dat]')
    print('bunny.dat + [bunny_00000_000000000200.txt bunny_00001_000000000200.txt] '
          '-> [bunny_r_000000000200.dat]')
    print('            [bunny_00000_000000000300.txt bunny_00001_000000000300.txt] '
          '-> [bunny_r_000000000300.dat]')

# --------------------------------------------------------------------------------------------------


# Example of running merge.py script:
#     merge.py grids/bunny.dat dir_to_txt_files
# Script produces out files for crys-remesh tool.
if __name__ == '__main__':

    if len(sys.argv) == 1:
        print_help()
        exit(0)

    if (sys.argv[1] == '-h') or (sys.argv[1] == '--help'):
        print_help()
        exit(0)

    if len(sys.argv) < 4:
        raise Exception('crys-gsu-merge : not enough arguments')

    merge(grid_file=sys.argv[1],
          txt_files_dir=sys.argv[2],
          r_files_dir=sys.argv[3])

# --------------------------------------------------------------------------------------------------
