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
    print(all_files)

    # We have to analyze only our *.txt files.
    match_files = [fn for fn in all_files if utils.is_filename_correct_crys_txt_file(fn, stem)]
    print(match_files)

    # Add all faces data from cry files.
    g = gsu.Grid()
    g.load(grid_file)
    for fn in match_files:
        g.load_faces_t_hw_hi(txt_files_dir + '/' + fn)
    g.store(grid_file + '.2')

    print('crys-gsu : merge_txt_from_mpi : done')

# ----------------------------------------------------------------------------------------------------------------------


# Example of running merge.py script:
#     merge.py grid/bunny.dat dir_to_txt_files
# Script produces out files for crys-remesh tool.
if __name__ == '__main__':

    merge_txt_from_mpi(grid_file=sys.argv[1], txt_files_dir=sys.argv[2])

# ----------------------------------------------------------------------------------------------------------------------
