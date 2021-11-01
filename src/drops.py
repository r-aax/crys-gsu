"""
Drops realization.
"""

import sys
import os
import gsu
import numpy as np
import utils
import time
import math
from geom.vect import Vect
from geom.segment import Segment
from geom.triangle import Triangle
from geom.trajectory import Trajectory
from geom.box import Box

# ==================================================================================================


def face_triangle(f):
    """
    Get face triangle.
    :param f: Face.
    :return: Triangle.
    """

    return Triangle(Vect.from_iterable(f.Nodes[0].P),
                    Vect.from_iterable(f.Nodes[1].P),
                    Vect.from_iterable(f.Nodes[2].P))

# ==================================================================================================


class SpacePartition:
    """
    Single space partitio.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, ds):
        """
        Constructor.
        :param ds: Data array.
        """

        self.Ds = ds
        self.Bx = Box([d[0] for d in ds])

# ==================================================================================================


class SpaceSeparator:
    """
    Class for space separation into small partitions.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, ds):
        """
        Constructor.
        :param ds: Data array.
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
            print('      {0}'.format(p.Bx))

    # ----------------------------------------------------------------------------------------------

    def find_nearest(self, p):
        """
        Find nearest data.
        :param p: Point.
        :return: Nearest data.
        """

        #
        # Warning! This works only if one partition is present.
        #

        m = np.array([(p - pi).mod2() for (pi, _) in self.Partitions[0].Ds])

        return self.Partitions[0].Ds[m.argmin()]

    # ----------------------------------------------------------------------------------------------

    def inside(self, p):
        """
        Check if point is inside.
        :param p: Point.
        :return: True - if point in box, False - otherwise.
        """

        return self.Partitions[0].Bx.is_inside(p)

    # ----------------------------------------------------------------------------------------------

    def fly_step(self, p, v, v_air, d, dt):
        """
        Step of flying.
        :param p: Point.
        :param v: Velocity.
        :param v_air: Air velocity.
        :param d: Diameter.
        :param dt: Time step.
        :return: New point position and new velocity.
        """

        # Calculate new point with old velocity.
        new_p = p + v * dt

        # Calculate new velocity.
        vis = 1.4607 * 0.00001
        re = (v - v_air).mod() * d / vis
        if re <= 350.0:
            cd = (24.0 / re) * (1 + 0.166 * math.pow(re, 0.33))
        else:
            cd = 0.178 * math.pow(re, 0.217)
        k = (3.0 / 4.0) * ((cd * 1.3) / (d * 1000.0)) * (v - v_air).mod() * dt
        vv = v_air - v
        new_v = v + vv * k

        return new_p, new_v

    # ----------------------------------------------------------------------------------------------

    def fly(self, p, vel, d, dt, g, max_steps):
        """
        Flying of a point.
        :param p: Point.
        :param vel: Velocity.
        :param d: Diameter.
        :param dt: Time step.
        :param g: Grid (surface).
        :param max_steps: Max count of fly steps.
        :return: Tuple of 3 elements.
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
        tr = Trajectory(cp)
        v = vel
        i = 0

        while self.inside(cp):

            np = self.find_nearest(cp)
            new_cp, v = self.fly_step(cp, v, np[1], d, dt)

            # If point stay on one place we can exit.
            if (cp - new_cp).mod2() < 1e-15:
                return ('S', None, tr)

            tr.add_point(new_cp)

            # Check intersection.
            for f in g.Faces:
                tri = face_triangle(f)
                if tri.intersection_with_segment(Segment(cp, new_cp)) != []:
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
    :param grid_air_file: File with air grid.
    :return: Space separator.
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
                        p = Vect(float(d[0]), float(d[1]), float(d[2]))
                        v = Vect(float(d[5]), float(d[6]), float(d[7]))
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
          d=1.0e-4, dt=1.0e-6, stall_thr=1.0e-6, max_fly_steps=200):
    """
    Calculate drops.
    :param grid_file: File with grid.
    :param grid_air_file: File with air grid.
    :param out_grid_file: Out file.
    :param d: Distance from face for flying point start.
    :param dt: Time step.
    :param stall_thr: Threshold for Stall value, to make decision about water stall.
    :param max_fly_steps: Max fly steps.
    """

    # Check for grid file.
    if not os.path.isfile(grid_file):
        raise Exception('crys-gsu-drops : no such file ({0})'.format(grid_file))

    # Check for grid file.
    if not os.path.isfile(grid_air_file):
        raise Exception('crys-gsu-drops : no such air file ({0})'.format(grid_file))

    print('crys-gsu-drops : start, grid_file = {0}, grid_air_file = {1}, '
          'out_grid_file = {2}, d = {3}, dt = {4}, '
          'stall_thr = {5}, max_fly_steps = {6}'.format(grid_file, grid_air_file, out_grid_file,
                                                        d, dt, stall_thr, max_fly_steps))
    start_time = time.time()

    # Load grid.
    g = gsu.Grid()
    g.load(grid_file)

    # Indexes of Stall and MImp2 fields.
    stall_ind = g.get_variable_index('Stall')
    stall_d_ind = g.get_variable_index('StallD')
    stall_vx_ind = g.get_variable_index('StallVX')
    stall_vy_ind = g.get_variable_index('StallVY')
    stall_vz_ind = g.get_variable_index('StallVZ')
    mimp2_ind = g.get_variable_index('MImp2')

    # Read air from file.
    air = read_vel_field_from_file(grid_air_file)
    air.print_info()

    tr_f = open(out_grid_file + '.tr.dat', 'w')
    tr_f.write('# Trajectories file.\n')
    tr_f.write('TITLE="Trajectories"\n')
    tr_f.write('VARIABLES="X", "Y", "Z"\n')

    # Check all faces.
    for f in g.Faces:
        stall_value = f.Data[stall_ind - 3]
        if stall_value > stall_thr:
            # print('... face {0} is wet. Start flying.'.format(f.GloId))
            stall_d = f.Data[stall_d_ind - 3]
            stall_vel = Vect(f.Data[stall_vx_ind - 3],
                             f.Data[stall_vy_ind - 3],
                             f.Data[stall_vz_ind - 3])
            tri = face_triangle(f)
            start_point = tri.centroid() + tri.normal_orth() * d
            res = air.fly(start_point, stall_vel, stall_d, dt, g, max_fly_steps)
            traj = res[2]
            print(traj)
            traj.dump(tr_f, 'Trajectory-{0}'.format(f.GloId))

            if res[0] == 'C':
                print('... secondary impingement '
                      'from face {0} to face {1}'.format(f.GloId, res[1].GloId))
                tri2 = face_triangle(res[1])
                res[1].Data[mimp2_ind - 3] = (stall_value / tri2.area())

    # Save grid back.
    g.convert_grid_stall_to_check_point()
    g.clean_t_hw_hi()
    g.store(out_grid_file)

    tr_f.close()

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
    print('    <d> - distance above face surface for start point of trajectory (default = 1.0e-4)')
    print('    <dt> - time step (default = 1.0e-6)')
    print('    <stall_thr> - threshold for stall faces (default = 1.0e-6)')
    print('    <max_fly_steps> - maximum steps count for drops flying (default = 200)')

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
              stall_thr=float(sys.argv[6]),
              max_fly_steps=int(sys.argv[7]))
    else:
        raise Exception('crys-gsu-drops : wrong parameters count')

# ==================================================================================================
