"""
Merge *.txt files with original grid.
"""

import sys
import os
import pathlib
import utils
import gsu

# ----------------------------------------------------------------------------------------------------------------------


def merge_txt_from_mpi(grid_file, txt_files_dir):
    """
    Merge txt files with original grid.
    Merge files separately for each timestamp and produce out _r_ file for crys-remesh tool.
    :param grid_file: file with original grid
    :param txt_files_dir: dir with txt files
    """

    print('crys-gsu : merge_txt_from_mpi : '
          'grid file = {0}, txt_files_dir = {1}'.format(grid_file, txt_files_dir))

    pp = pathlib.PurePath(grid_file)
    stem = pp.stem
    all_files = os.listdir(txt_files_dir)

    # We have to analyze only our *.txt files.
    match_files = [fn for fn in all_files if utils.is_filename_correct_crys_txt_file(fn, stem)]

    # Add all faces data from cry files.
    g = gsu.Grid()
    g.load(grid_file)
    for fn in match_files:
        g.load_faces_t_hw_hi(txt_files_dir + '/' + fn)
    nm, ts = utils.get_filename_and_timestamp_pair(stem)
    g.store('{0}_r_000000000100.dat'.format(stem))

    print('crys-gsu : merge_txt_from_mpi : done (data from {0} files is added)'.format(len(match_files)))

# ----------------------------------------------------------------------------------------------------------------------


# Example of running merge.py script:
#     merge.py grid/bunny.dat dir_to_txt_files
# Script produces out files for crys-remesh tool.
if __name__ == '__main__':

    merge_txt_from_mpi(grid_file=sys.argv[1], txt_files_dir=sys.argv[2])

# ----------------------------------------------------------------------------------------------------------------------
