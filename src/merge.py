"""
Merge *.txt files with original grid.
"""

import sys
import os
import pathlib
import utils
import gsu

# ----------------------------------------------------------------------------------------------------------------------


def merge(grid_file, txt_files_dir, r_files_dir):
    """
    Merge txt files with original grid.
    Merge files separately for each timestamp and produce out _r_ file for crys-remesh tool.
    :param grid_file: file with original grid
    :param txt_files_dir: dir with txt files
    :param r_files_dir: out dir for _r_ files
    """

    pp = pathlib.PurePath(grid_file)
    # Get characteristics of file:
    #   bs - base of path
    #   nm - name without timestamp
    #   ts - timestamp string
    #   sf - suffix
    bs, sf = str(format(pp.parents[0])), pp.suffix
    nm, ts = utils.get_filename_and_timestamp_pair(pp.stem)

    print('crys-gsu-merge : grid file = {0} (bs={1}, nm={2}, ts={3}, sf={4}), '
          'txt_files_dir={5}, r_files_dir={6}'.format(grid_file, bs, nm, ts, sf, txt_files_dir, r_files_dir))

    all_files = os.listdir(txt_files_dir)

    # We have to analyze only our *.txt files.
    match_files = [fn for fn in all_files if utils.is_filename_correct_crys_txt_file(fn, nm)]

    # Group match files.
    gmf = utils.group_txt_files_by_timestamps(match_files)

    # Add all faces data from cry files.
    os.makedirs(r_files_dir)
    g = gsu.Grid()
    for tm in gmf:
        fs = gmf.get(tm)
        g.load(grid_file)
        for f in fs:
            g.load_faces_t_hw_hi(txt_files_dir + '/' + f)
        g.store('{0}/{1}_r_{2}.dat'.format(r_files_dir, nm, tm))

    print('crys-gsu-merge : done ({0} _r_ files generated)'.format(len(gmf.keys())))

# ----------------------------------------------------------------------------------------------------------------------


# Example of running merge.py script:
#     merge.py grid/bunny.dat dir_to_txt_files
# Script produces out files for crys-remesh tool.
if __name__ == '__main__':

    if sys.argv[1] == '-h':
        print('merge.py <grid-file> <txt-files-dir> <r-files-dir>')
        exit(0)

    merge(grid_file=sys.argv[1],
          txt_files_dir=sys.argv[2],
          r_files_dir=sys.argv[3])

# ----------------------------------------------------------------------------------------------------------------------
