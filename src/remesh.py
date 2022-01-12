"""
Module describing Tong remesher.
"""

import argparse
import array
import os
from gsu.gsu import Grid
from remeshing import TongRemesher
import time
import struct


# ==================================================================================================


def check_argument(name):
    """
    Check argument.

    Parameters
    ----------
        name : string
            name of argument
    """

    if not os.path.isfile(name):
        print('File {} does not exist'.format(name))
        exit(1)
    else:
        if not name[-4:] == '.dat':
            print('File {} should be .dat file'.format(name))
            exit(1)


# --------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-g', '--grid', help='grid in .dat format')
group.add_argument('-d', '--dir', help='dir with .dat grids to read all')
parser.add_argument('-j', '--json', help='json file')
parser.add_argument('-o', '--outdir', help='all result grids will be stored here', default='.')
parser.add_argument('-v', '--verbosity', action='count',
                    help='increase output verbosity', default=0)
args = parser.parse_args()

if args.dir:
    source_dir = args.dir
    meshes = [source_dir + '/' + path for path in os.listdir(source_dir)]
    if args.verbosity > 0:
        if len(meshes) == 0:
            print('WARNING: no meshes found')

if args.grid:
    grid = args.grid
    meshes = [grid]

if args.json is not None:
    if not os.path.isfile(args.json):
        raise Exception('No such json-file {0}'.format(args.json))

if args.outdir:
    outdir = args.outdir

start = time.time()

# Create out dir.
try:
    os.makedirs(outdir)
except FileExistsError:
    pass

meshes = [mesh for mesh in meshes if mesh.endswith('.dat')]

for i, mesh in enumerate(meshes):
    grid = Grid()
    grid.load(mesh)

    if args.verbosity > 1:
        print('Mesh #{} is read'.format(i))

    if mesh.find('/') != -1:
        filename = mesh[mesh.rfind('/'):]
    else:
        filename = mesh

    remesher = TongRemesher(filename, args.json, outdir, verbose=True if args.verbosity > 1 else False)
    remesher(grid)

    outputfilename = outdir + '/' + filename.replace('_r_', '_')
    grid.store(outputfilename)

    if args.verbosity > 1:
        print('Mesh #{} was written into {}'.format(i, outputfilename))

if args.verbosity > 1:
    print('Total time: {} s'.format(time.time() - start))

# ==================================================================================================
