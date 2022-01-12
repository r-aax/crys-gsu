"""
This module describes classes of Remeshers.

Remesher allows to perform rebuilding of surface based on
heights of ice obtained during the computer modeling.

Shumilin Sergei noisd@yandex.ru
"""

from gsu.gsu import Grid
from gsu.node import Node
from geom.vect import Vect
from geom.utils import prizmatoid_volume_coefs, displaced_triangle_volume
from smoothing import NullSpaceSmoothing
import json
import statistics
import time
import math
from copy import deepcopy

# ==================================================================================================


def zipwith(a, b, f):
    """
    Zip two lists with a function.

    Parameters
    ----------
        a : list
            first list
        b : list
            second list
        f : function
            zip function

    Return
    ------
        New list.
    """

    return [f(ai, bi) for (ai, bi) in zip(a, b)]


# --------------------------------------------------------------------------------------------------


def wsum(a, w):
    """
    Weighted sum.

    Parameters
    ----------
        a : list
        w : weights

    Return
    ------
        Weighted sum.
    """

    return sum(zipwith(a, w, lambda ai, wi: ai * wi)) / sum(w)

# ==================================================================================================


class Remesher:
    """
    Base remesher of triangular surface grid.
    """

# --------------------------------------------------------------------------------------------------

    def __init__(self, n_rebuild_steps=1):
        if n_rebuild_steps < 1:
            raise ValueError('Wrong argument n_rebuild_steps')

        self.n_rebuild_steps = n_rebuild_steps

# --------------------------------------------------------------------------------------------------

    def __call__(self, grid : Grid):
        """
        Remesh the grid.

        Parameters
        ----------
            grid : Grid obj, triangular grid with Nodes, Faces lists
                   each face containes Hi param - height of ice.
        """
        raise NotImplementedError

# --------------------------------------------------------------------------------------------------

    def distance(self, node : Node) -> float:
        """
        Magnitude of shift.

        Parameters
        ----------
            node : Node obj
                the node being shifted

        Return
        ------
            distance : float
                magnitude of shift
        """
        raise NotImplementedError

# --------------------------------------------------------------------------------------------------

    def direction(self, node : Node) -> Vect:
        """
        Direction of shift.

        Parameter
        ---------
            node : Node obj
                the node being shifted

        Return
        ------
            direction : Vector obj  
                direction of shift, unit vector
        """
        raise NotImplementedError


# ==================================================================================================


class TongRemesher:
    """
    Remesher based on Tong's article.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, mesh_filename, json_file, outdir, verbose=False):
        """
        Init Tong remesher.

        Parameters
        ----------
            mesh_filename : string
                name of mesh file
            json_file : string
                file with parameters (in json format)
            outdir : string
                output directory
        """

        self.mesh_filename = mesh_filename
        self.json_file = json_file
        self.outdir = outdir
        self.set_default_options()

        if json_file is not None:

            with open(json_file, 'r') as jf:
                data = json.load(jf)
                jf.close()

            self.try_to_read_options(data)

        if verbose:
            self.verbose = True
    # ----------------------------------------------------------------------------------------------

    def set_default_options(self):
        """
        Set default options.
        """

        self.export_stl = False
        self.verbose = False
        self.nsteps_min = 0
        self.nsteps_max = 10
        self.nsteps_hi_side_fact = 0.1
        self.nsteps_ignore_part = 0.05
        self.nsmooth_steps = 0
        self.nsmooth_s = 10.0
        self.nsmooth_k = 0.15
        self.hsmooth_steps = 0
        self.hsmooth_alfa = 0.2
        self.hsmooth_beta = 0.1
        self.rest_volume_redistr = True
        self.nss_steps = 0
        self.nss_epsilon = 0.01
        self.nss_st = 0.2
        self.nss_fix_nodes = True
        self.pseudo_2d_oz = False

    # ----------------------------------------------------------------------------------------------

    def try_to_read_options(self, data):
        """
        Try to read options from json data.

        Parameters
        ----------
            data : json-data
                structure with json-data
        """

        # Read option export_stl.
        try:
            self.export_stl = data['export_stl']
        except KeyError:
            pass

        # Read verbose option.
        try:
            self.verbose = data['remesh_verbose']
        except KeyError:
            pass

        # Read nsteps_min option.
        try:
            self.nsteps_min = data['remesh_nsteps_min']
        except KeyError:
            pass

        # Read nsteps_max option.
        try:
            self.nsteps_max = data['remesh_nsteps_max']
        except KeyError:
            pass

        # Read nsteps_hi_side_fact option.
        try:
            self.nsteps_hi_side_fact = data['remesh_nsteps_hi_side_fact']
        except KeyError:
            pass

        # Read nstepd_ignore_part option.
        try:
            self.nsteps_ignore_part = data['remesh_nsteps_ignore_part']
        except KeyError:
            pass

        # Read nsmooth_steps.
        try:
            self.nsmooth_steps = data['remesh_nsmooth_steps']
        except KeyError:
            pass

        # Read nsmooth_s coefficient.
        try:
            self.nsmooth_s = data['remesh_nsmooth_s']
        except KeyError:
            pass

        # Read nsmooth_k coefficient.
        try:
            self.nsmooth_k = data['remesh_nsmooth_k']
        except KeyError:
            pass

        # Read hsmooth_spets parameter.
        try:
            self.hsmooth_steps = data['remesh_hsmooth_steps']
        except KeyError:
            pass

        # Read hsmooth_alfa coefficient.
        try:
            self.hsmooth_alfa = data['remesh_hsmooth_alfa']
        except KeyError:
            pass

        # Read hsmooth_beta parameter.
        try:
            self.hsmooth_beta = data['remesh_hsmooth_beta']
        except KeyError:
            pass

        # Read rest_volume_redistr parameter.
        try:
            self.rest_volume_redistr = data['remesh_rest_volume_redistr']
        except KeyError:
            pass

        # Read remesh_nss_steps parameter.
        try:
            self.nss_steps = data['remesh_nss_steps']
        except KeyError:
            pass

        # Read remesh_nss_epsilon parameter.
        try:
            self.nss_epsilon = data['remesh_nss_epsilon']
        except KeyError:
            pass

        # Read remesh_nss_st parameter.
        try:
            self.nss_st = data['remesh_nss_st']
        except KeyError:
            pass

        # Read remesh_nss_st parameter.
        try:
            self.nss_fix_nodes = data['remesh_nss_fix_nodes']
        except KeyError:
            pass

        # Read remesh_pseudo_2d_oz parameter.
        try:
            self.pseudo_2d_oz = data['remesh_pseudo_2d_oz']
        except:
            pass

    # ----------------------------------------------------------------------------------------------

    def __call__(self, grid: Grid):
        """
        Remesh the grid.

        Parameters
        ----------
            grid : Grid obj, triangular grid with Nodes, Faces lists
                   each face contains Hi param - height of ice.
        """

        self.grid = grid

        # Global parameters of accreted ice.
        for face in self.grid.Faces:

            # Ice we plan to grow.
            face.target_ice = face.get_triangle().area() * face['Hi']

            # Rest ice to grow.
            face.rest_ice = face.target_ice

            # Actual ice that has been accreted.
            face.actual_ice = 0.0

        if self.verbose:
            self.print_options()

        total_time = 0.0

        if self.nss_fix_nodes:
            grid.mark_all_fixed_nodes()

        # Define nsteps.
        nsteps = self.calculate_nsteps()

        # Mark nodes classes for option pseudo_2d_oz.
        if self.pseudo_2d_oz:
            self.clusterize_nodes_for_pseudo_2d_oz()
            self.freeze_nodes_z_coord()

        # Remesh steps.
        for i in range(nsteps):

            start_step_time = time.time()
            self.calc_ice_chunks(nsteps - i)
            self.remesh_step()
            finish_step_time = time.time()
            local_time = finish_step_time - start_step_time
            total_time += local_time

            if self.verbose:
                print('..... step {0} -- {1:.4} seconds'.format(i, local_time))

            # Debug export intermediate mesh.
            if self.verbose:
                self.grid.store(self.outdir + self.mesh_filename.replace('.dat',
                                                                       '_{0:05d}.dat'.format(i)))

        if self.verbose:
            print('Tong remesh total time : {0:.4} seconds'.format(total_time))
            total_target_ice = sum([face.target_ice for face in self.grid.Faces])
            total_ice_diff = sum(abs(face.actual_ice - face.target_ice) for face in self.grid.Faces)
            if total_target_ice != 0.0:
                perc = 100.0 * (total_ice_diff / total_target_ice)
            elif total_ice_diff == 0.0:
                perc = 0.0
            else:
                perc = math.inf
            print('Tong remesh accuracy   : {0} %'.format(perc))
            print('                   ice : '
                  '{0:.4} - target, {1:.4} - diff'.format(total_target_ice,
                                                          total_ice_diff))

        # Zero Hi.
        for face in self.grid.Faces:
            face['Hi'] = 0.0

    # ----------------------------------------------------------------------------------------------

    def print_options(self):
        """
        Print options.
        """

        print('Tong remesh with parameters:')
        print('    mesh_filename              = {0}'.format(self.mesh_filename))
        print('    json_file                  = {0}'.format(self.json_file))
        print('    outdir                     = {0}'.format(self.outdir))
        print('    export_stl                 = {0}'.format(self.export_stl))
        print('    remesh_verbose             = {0}'.format(self.verbose))
        print('    remesh_nsteps_min          = {0}'.format(self.nsteps_min))
        print('    remesh_nsteps_max          = {0}'.format(self.nsteps_max))
        print('    remesh_nsteps_hi_side_fact = {0}'.format(self.nsteps_hi_side_fact))
        print('    remesh_nsteps_ignore_part  = {0}'.format(self.nsteps_ignore_part))
        print('    remesh_nsmooth_steps       = {0}'.format(self.nsmooth_steps))
        print('    remesh_nsmooth_s           = {0}'.format(self.nsmooth_s))
        print('    remesh_nsmooth_k           = {0}'.format(self.nsmooth_k))
        print('    remesh_hsmooth_steps       = {0}'.format(self.hsmooth_steps))
        print('    remesh_hsmooth_alfa        = {0}'.format(self.hsmooth_alfa))
        print('    remesh_hsmooth_beta        = {0}'.format(self.hsmooth_beta))
        print('    remesh_rest_volume_redistr = {0}'.format(self.rest_volume_redistr))
        print('    remesh_nss_steps           = {0}'.format(self.nss_steps))
        print('    remesh_nss_epsilon         = {0}'.format(self.nss_epsilon))
        print('    remesh_nss_st              = {0}'.format(self.nss_st))
        print('    remesh_nss_fix_nodes       = {0}'.format(self.nss_fix_nodes))
        print('    remesh_pseudo_2d_oz        = {0}'.format(self.pseudo_2d_oz))

    # ----------------------------------------------------------------------------------------------

    def calculate_nsteps(self):
        """
        Calculate count of steps.

        Return
        ------
            Count of remesher steps.
        """

        if self.nsteps_min > self.nsteps_max:
            raise Exception('nsteps_max ({0}) can not be less '
                            'than nsteps_min ({1})'.format(self.nsteps_max,
                                                           self.nsteps_min))

        # We want to satisfy condition
        #   hi / (nsteps * min_side) < hi_side_fac
        #   hi / min_side < hi_side_fact * nsteps
        #   nsteps = [hi / (min_side * hi_side_fact)] + 1

        # Local nsteps array.
        lna = [int(f['Hi'] / (f.get_triangle().min_side() * self.nsteps_hi_side_fact)) + 1
               for f in self.grid.Faces]
        lna.sort()
        lna = list(reversed(lna))

        # Calculate count of faces we can ignore.
        lna_len = len(lna)
        may_ignore = int(lna_len * self.nsteps_ignore_part)

        if self.verbose:
            print('Calc nsteps: nsteps_min={0}, nsteps_max={1}, '
                  'nsteps_hi_side_fact={2}, '
                  'nsteps_ignore_part={3}'.format(self.nsteps_min,
                                                  self.nsteps_max,
                                                  self.nsteps_hi_side_fact,
                                                  self.nsteps_ignore_part))
            print('    may ignore {0} of {1} faces'.format(may_ignore, lna_len))

        # If after ignore we have too big value of nsteps, just use nsteps_max.
        nsteps = self.nsteps_max
        if lna[may_ignore] > self.nsteps_max:
            if self.verbose:
                print('    ignore bad faces did not help: '
                      'first value after ignore = {0}'.format(lna[may_ignore]))
        else:
            # Find first value that no more than nsteps_max.
            for i in range(may_ignore):
                if lna[i] <= self.nsteps_max:
                    nsteps = lna[i]
                    if self.verbose:
                        print('    take {0}-th value ({1}) '
                              'while ignoring bad faces'.format(i, nsteps))
                    break

        if nsteps < self.nsteps_min:
            print('    low nsteps threshold fires ({0})'.format(self.nsteps_min))
            nsteps = self.nsteps_min

        if self.verbose:
            print('    nsteps = {0}'.format(nsteps))

        return nsteps

    # ----------------------------------------------------------------------------------------------

    def clusterize_nodes_for_pseudo_2d_oz(self):
        """
        Clusterize nodes when using option pseudo_2d_oz.
        """

        big_value = 10 ** 20
        d = {}

        # Add nodes to appropriate classes with the same x, y coords.
        for node in self.grid.Nodes:
            c = (int(node.P.x * big_value), int(node.P.y * big_value))
            v = d.get(c)
            if v is not None:
                v.append(node)
            else:
                d[c] = [node]

        # Calc global parameters for all columns.
        vs = list(d.values())
        glo_len = len(vs[0])

        # Mark pseudo 2d oz mode columns and check by the fly.
        for i, v in enumerate(vs):
            if len(v) != glo_len:
                raise Exception('pseudo_2d_oz mode is unable for pure 3d mesh (diff columns lens)')
            if len(v) == 1:
                raise Exception('pseudo_2d_oz mode is unable for pure 3d mesh (column len = 1)')
            v.sort(key=lambda n: n.P.z)

        # Set information to grid.
        self.grid.pseudo_2d_oz_columns = vs

    # ----------------------------------------------------------------------------------------------

    def freeze_nodes_z_coord(self):
        """
        Freeze z coords for all nodes.
        """

        for n in self.grid.Nodes:
            n.frozen_z = n.P.z

    # ----------------------------------------------------------------------------------------------

    def calc_ice_chunks(self, rest_iterations):
        """
        Calculate ice chunks in faces.

        Parameters
        ----------
            rest_iterations : int
                iteration rest to calculate
        """

        # Calculate ice_chunk - ice to grow on current iteration.
        # Calculate ice chunk from rest ice and rest count of iterations.
        # Requirement FR.RM.MT.TN.04.
        for face in self.grid.Faces:
            face.ice_chunk = face.rest_ice / rest_iterations

    # ----------------------------------------------------------------------------------------------

    def null_space_smoothing(self):
        """Smooth the grid with null-space smoothing algorithm.
        """
        NullSpaceSmoothing(self.grid,
                           self.nss_steps,
                           fix_boundary_nodes=self.nss_fix_nodes,
                           epsilon=self.nss_epsilon,
                           st=self.nss_st).smoothing()

    def remember_current_areas(self):
        """Remember current area of the face""" 
        for f in self.grid.Faces:
            f.saved_area = f.get_triangle().area()

    def correction_by_area(self):
        """Correct height of ice inside a face based on:
                - area of the face before smoothing
                - current area of the face.
        """
        for face in self.grid.Faces:
            assert face.saved_area != 0.0
            face.rest_ice = face.rest_ice * face.get_triangle().area() / face.saved_area

    # ----------------------------------------------------------------------------------------------

    def remesh_step(self):
        """
        One step of the remeshing process.
        """

        # While deforming mesh we work with nods' and faces' ice_dirs.
        self.init_ice_dirs()

        # Local normals (ice_directions) smoothing.
        self.normals_smoothing()

        # Define heights field.
        self.define_ice_shifts()

        # Smoothing of heights.
        self.heights_smoothing()

        # Move nodes.
        self.move_nodes()

        # Remember squares of faces before nss.
        self.remember_current_areas()

        # Null-space smoothing.
        self.null_space_smoothing()

        # Correction by area.
        self.correction_by_area()

        # Correct coords for pseudo_2d_oz crutch.
        if self.pseudo_2d_oz:
            self.correct_coords_for_pseudo_2d_oz()

    # ----------------------------------------------------------------------------------------------

    def correct_coords_for_pseudo_2d_oz(self):
        """
        Correct coords for pseudo_2d_oz.
        """

        cols = self.grid.pseudo_2d_oz_columns
        mid_index = len(cols[0]) // 2

        # Correct each column.
        for col in cols:
            for i, n in enumerate(col):
                if i != mid_index:
                    n.P.X, n.P.Y, n.P.Z = col[mid_index].x, col[mid_index].y, n.frozen_z

    # ----------------------------------------------------------------------------------------------

    def init_ice_dirs(self):
        """
        Initialize ice directions.
        """

        # By default ice direction for face is its normal.
        for face in self.grid.Faces:
            face.ice_dir = deepcopy(face.get_triangle().normal_orth())

        # Requirement FR.RM.MT.TN.02.
        # Direction of ice growth in node is mean of all adjacent faces normals.
        for node in self.grid.Nodes:
            node.ice_dir = sum([f.ice_dir for f in node.Faces], Vect())
            node.ice_dir = node.ice_dir.orth()

    # ----------------------------------------------------------------------------------------------

    def normals_smoothing(self):
        """
        Smoothing of normals.
        Requirement FR.RM.MT.TN.03.
        """

        s = self.nsmooth_s
        k = self.nsmooth_k

        for _ in range(self.nsmooth_steps):

            # First, we smooth faces' ice directions through nodes' ice directions.
            for face in self.grid.Faces:
                ws = [max(s * (1.0 - (face.ice_dir * node.ice_dir)), k) for node in face.Nodes]
                face.ice_dir = wsum([node.ice_dir for node in face.Nodes], ws)

            # Second, we smooth nodes' ice directions through faces' ice directions.
            for node in self.grid.Nodes:
                ws = [1.0 / face.get_triangle().area() for face in node.Faces]
                node.ice_dir = wsum([face.ice_dir for face in node.Faces], ws)

    # ----------------------------------------------------------------------------------------------

    def define_ice_shifts(self):
        """
        Define ice shifts in cells and nodes.
        For calculation ice_shifts we use cur_ice_chunks.
        Requirement FR.RM.MT.TN.05.
        """

        eps = 1.0e-10

        # First define ice shifts in faces.
        # We use formula from Tong's article for obtaining volume of prizmatoid from its height.
        # V = a * h + b * h * h + c * h * h * h
        # We ignore c * h * h * h member of equation (it is o-small of overall value).
        # V = a * h + b * h * h
        # b * h * h + a * h - V = 0
        for f in self.grid.Faces:

            # Classic calculation.
            f.ice_shift = f.ice_chunk / f.get_triangle().area()

            # Trying to calculate more accurately.
            if f.ice_chunk > eps:
                a, b, _ = prizmatoid_volume_coefs(f.Nodes[0].P,
                                                  f.Nodes[1].P,
                                                  f.Nodes[2].P,
                                                  f.Nodes[0].ice_dir,
                                                  f.Nodes[1].ice_dir,
                                                  f.Nodes[2].ice_dir,
                                                  f.get_triangle().normal_orth())

                # Check if it is needed to solve square equation.
                if abs(b) > eps:
                    d = a * a + 4.0 * b * f.ice_chunk
                    if d >= 0.0:
                        d = math.sqrt(d)
                        h1, h2 = (-a + d) / (2.0 * b), (-a - d) / (2.0 * b)
                        if (h1 >= 0.0) and (h2 >= 0.0):
                            f.ice_shift = min(h1, h2)
                        elif h1 >= 0.0:
                            f.ice_shift = h1
                        elif h2 >= 0.0:
                            f.ice_shift = h2

        # Define ice shifts for nodes.
        for n in self.grid.Nodes:
            n.ice_shift = statistics.mean([f.ice_shift for f in n.Faces])
            if n.ice_shift < 0.0:
                raise Exception('Negative node ice shift.')

# --------------------------------------------------------------------------------------------------

    def heights_smoothing(self):
        """
        Heights smoothing.
        Requirement FR.RM.MT.TN.06.
        """

        # Smoothing iterations.
        for _ in range(self.hsmooth_steps):

            # While redistributing we work with local ice chunks.
            for face in self.grid.Faces:
                face.loc_ice_chunk = face.ice_chunk

            # Calculate max h.
            max_h = max([face.ice_shift for face in self.grid.Faces])

            # Process all inner edges.
            for edge in self.grid.Edges:
                if len(edge.Faces) == 2:
                    f1, f2 = edge.Faces[0], edge.Faces[1]

                    # Suppose h1 >= h2.
                    if f1.ice_shift < f2.ice_shift:
                        f1, f2 = f2, f1

                    # Calculate delta volume.
                    if f1.ice_shift == 0.0:
                        mid_area = f1.get_triangle().area()
                    else:
                        mid_area = f1.ice_chunk / f1.ice_shift
                    delta_v = min(f1.ice_shift - f2.ice_shift,
                                  self.hsmooth_alfa * max_h) * mid_area

                    # Redistribute volume.
                    f1.loc_ice_chunk -= self.hsmooth_beta * delta_v
                    f2.loc_ice_chunk += self.hsmooth_beta * delta_v

            # Put current ice chunks back.
            for face in self.grid.Faces:
                face.ice_chunk = face.loc_ice_chunk

            # Recalc ice_shifts.
            self.define_ice_shifts()

    # ----------------------------------------------------------------------------------------------

    def move_nodes(self):
        """
        Move nodes
        Requirement FR.RM.MT.TN.07.
        """

        # First we have to calculate all shifts and after that we can move nodes.
        for node in self.grid.Nodes:
            node.shift = node.ice_dir * node.ice_shift

        # Update actual ice value.
        for face in self.grid.Faces:
            n1, n2, n3 = face.Nodes[0], face.Nodes[1], face.Nodes[2]
            p1, p2, p3 = n1.P, n2.P, n3.P
            np1 = p1 + n1.shift
            np2 = p2 + n2.shift
            np3 = p3 + n3.shift
            v = displaced_triangle_volume(p1, p2, p3, np1, np2, np3)

            # Requirement FR.RM.MT.TN.08
            # Save volume difference that we will try to accrete on the next step.
            # Apply this step only when option in on.
            if self.rest_volume_redistr:
                face.rest_ice -= v
                if face.rest_ice < 0.0:
                    face.rest_ice = 0.0

            face.actual_ice += v

        # Move nodes.
        for node in self.grid.Nodes:
            node.P = node.P + node.shift

# ==================================================================================================
