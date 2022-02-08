from scipy.linalg import eig, det
from numpy import argmax, array, vstack, diag, dot, abs, cumprod, sum, full, arccos, exp, real
from copy import deepcopy
from collections import deque
import geom
import gsu

DETERMINANT_ACCURACY = 10e-5
EPSILON = 10e-5

# ==================================================================================================


class Smoothing:
    __name__ = ''

    def __init__(self, grid, num_interations=0, fix_boundary_nodes=False):
        """Constructor.

        Parameters
        ----------
            grid : Grid object
              grid to smooth
            num_itarations : int
              number of smoothing iterations
            fix_boundary_nodes : bool
              whether to keep boundary nodes fixed
        """
        self.grid = grid
        self.num_iterations = num_interations
        self.fix_boundary_nodes = fix_boundary_nodes

# --------------------------------------------------------------------------------------------------

    @staticmethod
    def name_of_iteration(i):
        """Compose name for iteration."""
        iteration_name = str(i)
        zeros = '000'
        return zeros[:-len(iteration_name)] + iteration_name

# --------------------------------------------------------------------------------------------------

    def smoothing(self):
        """Perdorm smoothing."""
        raise NotImplementedError

# --------------------------------------------------------------------------------------------------

    def move_node(self, node, shift: geom.Vect):
        """Move node for a given shift.

        Parameters
        ----------
            node : Node object
              node
            shift : Vector object
              shift
        """
        if self.fix_boundary_nodes:
            if node.border:
                pass
            else:
                node.P = node.P + shift
        else:
             node.P = node.P + shift

# --------------------------------------------------------------------------------------------------

    def write_grid_and_print_info(self, iteration):
        """Write info and save current grid.
        
        Parameters
        ----------
            iteration : int
              number of iteration
        """
        print('{}th iteration of {} smoothing'.format(iteration, self.__name__))
        self.grid.store('{}_smoothing_{}.dat'.format(self.__name__,
                                                              self.name_of_iteration(iteration)))

# --------------------------------------------------------------------------------------------------

    def apply_shifts(self, shifts):
        """Apply precomputed list of laplacians to the nodes.
        
        Parameters
        ----------
            laplacians : list of Vector objects
              shifts for nodes
        """
        assert len(self.grid.Nodes) == len(shifts)
        for n, l in zip(self.grid.Nodes, shifts):
            self.move_node(n, l)

# ==================================================================================================


class NullSpaceSmoothing(Smoothing):
    __name__ = "NullSpace"

    def __init__(self, grid: gsu.Grid,
                 num_iterations=0, st=0.2, epsilon=1e-2,
                 print_intermediate_steps=False,
                 fix_boundary_nodes=False):
        """Constructor.

        Parameters
        ----------
            grid : Grid object
              grid to smooth
            num_itarations : int
              number of smoothing iterations
            st : float
              control of magnitude of node's shift
            epsilon : float
              control of degree of smoothing
            print_intermediate_steps : bool
              print intermediate steps of smoothing
            fix_boundary_nodes : bool
              whether to keep boundary nodes fixed
        """
        assert len(grid.Nodes) > 0, 'the grid is empty'
        assert len(grid.Faces) > 0, 'the grid is empty'
        assert len(grid.Edges) > 0, 'the grid is empty'
        self.st = st
        self.epsilon = epsilon
        self.print_intermediate_steps = print_intermediate_steps
        Smoothing.__init__(self, grid, num_iterations, fix_boundary_nodes)

# --------------------------------------------------------------------------------------------------

    @staticmethod
    def create_matrix_of_normals(node) -> array:
        """Create matrix consisting of adjacent faces' normals.
        
        Parameters
        ----------
            node : Node object
              node
        
        Returns
        -------
            numpy array (N_FACES, 3)
        """
        N = node.Faces[0].get_triangle().normal_orth().coords_as_numpy()
        N = N.reshape((1, 3))
        for f in node.Faces[1:]:
            normal = f.get_triangle().normal_orth().coords_as_numpy()
            assert N.shape[1] == normal.shape[1], print(N.shape, normal.shape)
            N = vstack((N, normal))

        return N

# --------------------------------------------------------------------------------------------------

    def get_laplacian(self, node):
        """Calculate the vector to the sum of centroids from a node.

        Parameters
        ----------
            node : Node object
              target node
        
        Returns
        -------
            Vector obj
        """
        dv = Vect()
        p = node.P
        weights = 0

        neighbours = node.Faces

        assert len(neighbours) > 0

        for f in neighbours:
            c = f.get_triangle().centroid()
            c = c - p
            weight = 1.0

            assert isinstance(c, Vect)

            c = c * weight
            dv = dv + c
            weights += weight

        dv = dv / weights
        return dv

# --------------------------------------------------------------------------------------------------

    def smoothing(self):
        """Perform smoothing."""
        for i in range(0, self.num_iterations):
            shifts = []
            for n in self.grid.Nodes:
                m = len(n.faces)
                w = [f.get_triangle().area() for f in n.faces]
                N = self.create_matrix_of_normals(n)
                dv = self.get_laplacian(n)

                assert N.shape == (m, 3)
                assert len(w) == m
                assert isinstance(dv, Vect)

                W = diag(w)
                NTW = dot(N.T, W)
                A = dot(NTW, N)

                assert W.shape == (m, m)
                assert NTW.shape == (3, m)
                assert A.shape == (3, 3)

                eigenValues, eigenVectors = eig(A)
                idx = eigenValues.argsort()[::-1]
                eigenValues = eigenValues[idx]
                eigenVectors = eigenVectors[:, idx]

                assert abs(det(A) - cumprod(eigenValues)[-1]) < DETERMINANT_ACCURACY, \
                    print(det(A), cumprod(eigenValues)[-1])

                k = sum((eigenValues > self.epsilon * eigenValues[0]))
                ns = eigenVectors[:, k:]
                dv = dv.coords_as_np().reshape(3, 1)

                assert ns.shape == (3, 3 - k), 'Wrong eigenvectors shape'
                assert dv.shape == (3, 1), 'Wrong shape'

                if k < 3:
                    ttT = dot(ns, ns.T)
                    t = self.st * dot(ttT, dv)
                    assert t.shape == (3, 1), 'Wrong shift shape'
                    shift = Vect(real(t[0, 0]), real(t[1, 0]), real(t[2, 0]))
                else:
                    shift = Vect(0, 0, 0)
                shifts.append(shift)

            Smoothing.apply_shifts(self, shifts)

            if self.print_intermediate_steps:
                Smoothing.write_grid_and_print_info(self, i)

# ==================================================================================================
