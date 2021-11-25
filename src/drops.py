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
from geom.triangles_cloud import TrianglesCloud
from geom.trajectory import Trajectory
from geom.box import Box

# ==================================================================================================


class SpacePartition:
    """
    Single space partition.
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

    def __init__(self, ds, ds_tuple):
        """
        Constructor.
        :param ds: Data array.
        """

        # By default we make one single partition.

        self.Partitions = [SpacePartition(ds)]
        self.PartitionsTuple = np.array(ds_tuple)

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
        :return:  Nearest data.
        """

        #
        # Warning! This works only if one partition is present.
        #

        # m = np.array([(p - pi).mod2() for (pi, _) in self.Partitions[0].Ds])
        #
        # return self.Partitions[0].Ds[m.argmin()]

        find_point = np.array(p.coords_tuple())
        res = self.PartitionsTuple - find_point
        distances = np.linalg.norm(res, axis=1)
        min_index = np.argmin(distances)

        return self.Partitions[0].Ds[min_index]

    # ----------------------------------------------------------------------------------------------

    def inside(self, p):
        """
        Check if point is inside.
        :param p: Point.
        :return:  True - if point in box, False - otherwise.
        """

        return self.Partitions[0].Bx.is_inside(p)

    # ----------------------------------------------------------------------------------------------

    def fly_step(self, p, v, v_air, d, dt):
        """
        Step of flying.
        :param p:     Point.
        :param v:     Velocity.
        :param v_air: Air velocity.
        :param d:     Diameter.
        :param dt:    Time step.
        :return:      New point position and new velocity.
        """

        # Air viscosity and density for zero temperature.
        air_dns = 1.292
        wtr_dns = 1000.0
        air_vis = 1.736e-5

        # Reynolds number and C_d.
        vsm = (v - v_air).mod()
        re = vsm * d / air_vis
        if re <= 350.0:
            cd = (24.0 / re) * (1.0 + 0.166 * math.pow(re, 0.33))
        else:
            cd = 0.178 * math.pow(re, 0.217)

        # Speedup.
        a = (v_air - v) * 0.75 * ((cd * air_dns) / (d * wtr_dns)) * vsm

        # Calculate new position through velocity,
        # and new velocity through speedup.
        new_p = p + v * dt
        new_v = v + a * dt

        return new_p, new_v

    # ----------------------------------------------------------------------------------------------

    def fly(self, p, vel, d, dt, triangles_cloud, max_steps):
        """
        Flying of a point.
        :param p:               Point.
        :param vel:             Velocity.
        :param d:               Diameter.
        :param dt:              Time step.
        :param triangles_cloud: Triangles cloud.
        :param max_steps:       Max count of fly steps.
        :return:                Tuple of 3 elements.
                                  first element - diagnostic
                                    'N' - too long flying
                                    'O' - left the box out
                                    'S' - stop on place
                                    'C' - cross surface
                                  second element - face if intersect (None otherwise)
                                  third element - trajectory
        """

        tr = Trajectory(p)
        v = vel

        # Infinite loop.
        while True:

            # Get last point (lp) in trajectory.
            lp = tr.last_point()

            # Check if there is too many points in trajectory.
            if tr.points_count() > max_steps:
                return ('N', None, tr)

            # Check if we still in our space.
            if not self.inside(lp):
                return ('O', None, tr)

            # Find air velocity.
            v_air = self.find_nearest(lp)[1]

            # Calculate future point (fp) and new velocity value.
            fp, v = self.fly_step(lp, v, v_air, d, dt)

            # If point stay on one place we can exit.
            if fp.is_near(lp, 1e-15):
                return ('S', None, tr)

            # Add point to trajectory.
            tr.add_point(fp)

            # Everything is OK.
            # Check intersection.
            tri = triangles_cloud.first_intersection_with_segment(Segment(lp, fp))
            if not tri is None:
                return ('C', tri.BackRef, tr)

# ==================================================================================================


def read_vel_field_from_file(grid_air_file):
    """
    Read velocity field from 3D grid.
    :param   grid_air_file: File with air grid.
    :return: Space separator.
    """

    ds = []
    ds_tuple = []

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
                        ds_tuple.append(p.coords_tuple())
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
    sep = SpaceSeparator(ds, ds_tuple)

    return sep

# --------------------------------------------------------------------------------------------------


def drops(grid_stall_file, grid_air_file, out_grid_file,
          d, dt, stall_thr, max_fly_steps):
    """
    Calculate drops.
    :param grid_stall_file: File with grid.
    :param grid_air_file:   File with air grid.
    :param out_grid_file:   Out file.
    :param d:               Distance from face for flying point start.
    :param dt:              Time step.
    :param stall_thr:       Threshold for Stall value, to make decision about water stall.
    :param max_fly_steps:   Max fly steps.
    """

    # Check for grid file.
    if not os.path.isfile(grid_stall_file):
        raise Exception('crys-gsu-drops : no such file ({0})'.format(grid_stall_file))

    # Check for grid file.
    if not os.path.isfile(grid_air_file):
        raise Exception('crys-gsu-drops : no such air file ({0})'.format(grid_air_file))

    print('crys-gsu-drops : start, grid_stall_file = {0}, grid_air_file = {1}, '
          'out_grid_file = {2}, d = {3}, dt = {4}, stall_thr = {5}, '
          'max_fly_steps = {6}'.format(grid_stall_file, grid_air_file, out_grid_file,
                                       d, dt, stall_thr, max_fly_steps))
    start_time = time.time()

    # Load grid.
    g = gsu.Grid()
    g.load(grid_stall_file)
    triangles_cloud = TrianglesCloud(g.get_triangles_list())

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
            tri = f.get_triangle()
            start_point = tri.centroid() + tri.normal_orth() * d
            res = air.fly(start_point, stall_vel, stall_d, dt, triangles_cloud, max_fly_steps)
            traj = res[2]
            print(traj)
            traj.dump(tr_f, 'Trajectory-{0}'.format(f.GloId))

            if res[0] == 'C':
                print('... secondary impingement '
                      'from face {0} to face {1}'.format(f.GloId, res[1].GloId))
                tri2 = res[1].get_triangle()
                res[1].Data[mimp2_ind - 3] = (stall_value / tri2.area())

    # Save grid back.
    g.convert_grid_stall_to_check_point()
    g.clean_t_hw_hi()
    g.store(out_grid_file)

    tr_f.close()

    print('crys-gsu-drops : done (time estimated = {0} s)'.format(time.time() - start_time))

# ==================================================================================================


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='drops',
                                     description='Drops trajectories and secondary impingement calculation.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('grid_stall_file', help='grid file name in STALL format')
    parser.add_argument('grid_air_file', help='name of file with air grid')
    parser.add_argument('out_grid_file', help='out file with result grid, trajectories are stored in <out_grid_file>.tr.dat')
    parser.add_argument('-d', '--distance', dest='distance', type=float, default=1.0e-4,
                        help='distance above face surface for start point of trajectory (m)')
    parser.add_argument('-t', '--time_delta', dest='time_delta', type=float, default=1.0e-5,
                        help='time step (s)')
    parser.add_argument('-s', '--stall_threshold', dest='stall_threshold', type=float, default=1.0e-6,
                        help='threshold for stall faces (kg / s)')
    parser.add_argument('-m', '--max_fly_steps', dest='max_fly_steps', type=int, default=200,
                        help='maximum points in droplets trajectory')
    args = parser.parse_args()

    # Run.
    drops(args.grid_stall_file, args.grid_air_file, args.out_grid_file,
          args.distance, args.time_delta, args.stall_threshold, args.max_fly_steps)

# ==================================================================================================
