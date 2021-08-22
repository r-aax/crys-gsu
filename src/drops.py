"""
Drops realization.
"""

import sys
import os
import gsu
import numpy as np
import utils
import time

# ==================================================================================================


class SpacePartition:
    """
    Single space partitio.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, ds):
        """
        Constructor.
        :param ds: data array
        """

        self.Ds = ds

        # Construct box.
        xs = [d[0][0] for d in ds]
        ys = [d[0][1] for d in ds]
        zs = [d[0][2] for d in ds]
        self.Box = (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))

# ==================================================================================================


class SpaceSeparator:
    """
    Class for space separation into small partitions.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, ds):
        """
        Constructor.
        :param ds: data array
        """

        # By default we make one single partition.

        self.Partitions = [SpacePartition(ds)]

    # ----------------------------------------------------------------------------------------------

    def print_info(self):
        """
        Print info.
        """

        print('SpaceSeparator : {0} partitions'.format(len(self.Partitions)))

        for i, p in enumerate(self.Partitions):
            print('    : part {0} : {1} points'.format(i, len(p.Ds)))
            print('      box {0}'.format(p.Box))

    # ----------------------------------------------------------------------------------------------

    def find_nearest(self, p):
        """
        Find nearest data.
        :param p: point
        :return: nearest data
        """

        #
        # Warning! This works only if one partition is present.
        #

        m = np.array([utils.dist2(p, pi) for (pi, _) in self.Partitions[0].Ds])

        return self.Partitions[0].Ds[m.argmin()]

    # ----------------------------------------------------------------------------------------------

    def inside(self, p):
        """
        Check if point is inside.
        :param p: point
        :return: True - if point in box, False - otherwise.
        """

        #
        # Warning! This works only if one partition is present.
        #

        (min_x, min_y, min_z, max_x, max_y, max_z) = self.Partitions[0].Box
        (px, py, pz) = p

        return (px >= min_x) and (px <= max_x)\
               and (py >= min_y) and (py <= max_y)\
               and (pz >= min_z) and (pz <= max_z)

    # ----------------------------------------------------------------------------------------------

    def fly(self, p, dt, g, max_steps):
        """
        Flying of a point.
        :param p: point
        :param dt: time step
        :param g: grid (surface)
        :param max_steps: max count of fly steps
        :return: tuple of 3 elements.
                 first element - diagnostic
                     'S' - stop on place
                     'C' - cross surface
                     'O' - left the box out
                     'N' - too long flying
                 second element - face if intersect (None otherwise)
                 third element - trajectory
        """

        # Trajectory is initilized with start point.
        cp = p
        tr = [cp]
        i = 0

        while self.inside(cp):

            np = self.find_nearest(cp)
            v = np[1]
            new_cp = utils.a_kb(cp, dt, v)

            # If point stay on one place we can exit.
            if utils.dist2(cp, new_cp) < 1e-10:
                return ('S', None, tr)

            tr.append(new_cp)

            # Check intersection.
            for f in g.Faces:
                if utils.is_triangle_and_segment_intersect(f.Nodes[0].P,
                                                           f.Nodes[1].P,
                                                           f.Nodes[2].P,
                                                           cp, new_cp):
                    return ('C', f, tr)

            cp = new_cp
            i = i + 1
            if i > max_steps:
                return ('N', None, tr)

        # We left box out.
        return ('O', None, tr)

# ==================================================================================================


def read_vel_field_from_file(grid_air_file):
    """
    Read velocity field from 3D grid.
    :param grid_air_file: file with air grid
    """

    ds = []

    with open(grid_air_file, 'r') as f:
        l = f.readline()
        while l:

            if 'TITLE=' in l:
                pass
            elif 'VARIABLES=' in l:
                pass
            elif 'ZONE' in l:
                ss = l.split()[1:]
                if 'FEBRICK' in ss[0]:
                    # Air data.
                    ns, es = int(ss[5]), int(ss[8])
                    for i in range(ns):
                        l = f.readline()
                        d = l.split()
                        p = float(d[0]), float(d[1]), float(d[2])
                        v = float(d[5]), float(d[6]), float(d[7])
                        ds.append((p, v))
                    # Ignore all febricks links.
                    for i in range(es):
                        l = f.readline()
                else:
                    # Not air data. Ignore all.
                    ns, es = int(ss[4]), int(ss[7])
                    for i in range(ns + es):
                        l = f.readline()
            else:
                print('unexpected line ({0})'.format(l))

            l = f.readline()

    f.close()

    # Now create space separator for search points.
    sep = SpaceSeparator(ds)

    return sep

# --------------------------------------------------------------------------------------------------


def drops(grid_file, grid_air_file, out_grid_file,
          d=0.0001, dt=0.00001, wet_thr=0.0004, max_fly_steps=10):
    """
    Calculate drops.
    :param grid_file: file with grid
    :param grid_air_file: file with air grid
    :param out_grid_file: out file
    :param d: distance from face for flying point start
    :param dt: time step
    :param wet_thr: thresjold for hw, to make decision about water stall
    :param max_fly_steps: max fly steps
    """

    # Check for grid file.
    if not os.path.isfile(grid_file):
        raise Exception('crys-gsu-drops : no such file ({0})'.format(grid_file))

    # Check for grid file.
    if not os.path.isfile(grid_air_file):
        raise Exception('crys-gsu-drops : no such air file ({0})'.format(grid_file))

    print('crys-gsu-drops : start, grid_file = {0}, grid_air_file = {1}, '
          'out_grid_file = {2}, d = {3}, dt = {4}, '
          'wet_thr = {5}, max_fly_steps = {6}'.format(grid_file, grid_air_file, out_grid_file,
                                                      d, dt, wet_thr, max_fly_steps))
    start_time = time.time()

    # Load grid.
    g = gsu.Grid()
    g.load(grid_file)

    # Clean all ice.
    for f in g.Faces:
        f.set_hi(0.0)

    # Read air from file.
    air = read_vel_field_from_file(grid_air_file)
    air.print_info()

    # Check all faces.
    for f in g.Faces:
        if f.get_hw() > wet_thr:
            # print('... face {0} is wet. Start flying.'.format(f.GloId))
            res = air.fly(f.get_point_above(d), dt, g, max_fly_steps)
            if res[0] == 'C':
                print('... secondary impingement '
                      'from face {0} to face {1}'.format(f.GloId, res[1].GloId))
                res[1].set_hi(1.0)

    # Save grid back.
    g.store(out_grid_file)

    print('crys-gsu-drops : done (time estimated = {0} s)'.format(time.time() - start_time))

# --------------------------------------------------------------------------------------------------


def print_help():
    """
    Print help.
    """

    print('[Overview]:')
    print('drops.py script calculates drops trajectories and secondary impingement')
    print('')
    print('[Usage]:')
    print('merge.py <grid-file> <grid-air-file> <out-grid-file> ')
    print('         [<d> <dt> <wet_thr> <max_fly_steps>]')
    print('    <grid-file> - grid file name')
    print('    <grid-air-file> - name of file with air grid')
    print('    <out-grid-file> - out file with result grid')
    print('    <d> - distance above face surface for start point of trajectory (default = 0.0001)')
    print('    <dt> - time step (default = 0.00001)')
    print('    <wet_thr> - threshold for stall faces (default = 0.0004)')
    print('    <max_fly_steps> - maximum steps count for drops flying (default = 10)')

# ==================================================================================================


# Example of running drops.py script:
#     drops.py grids/cyl.dat grids/cyl_air.dat grids/out_cyl.dat
if __name__ == '__main__':

    if len(sys.argv) == 1:
        print_help()
        exit(0)

    if (sys.argv[1] == '-h') or (sys.argv[1] == '--help'):
        print_help()
        exit(0)

    if len(sys.argv) < 4:
        raise Exception('crys-gsu-drops : not enough arguments')

    if len(sys.argv) == 4:
        drops(grid_file=sys.argv[1],
              grid_air_file=sys.argv[2],
              out_grid_file=sys.argv[3])
    elif len(sys.argv) == 8:
        drops(grid_file=sys.argv[1],
              grid_air_file=sys.argv[2],
              out_grid_file=sys.argv[3],
              d=float(sys.argv[4]),
              dt=float(sys.argv[5]),
              wet_thr=float(sys.argv[6]),
              max_fly_steps=int(sys.argv[7]))
    else:
        raise Exception('crys-gsu-drops : wrong parameters count')

# ==================================================================================================
