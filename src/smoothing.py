from crysremesh.geom import *
from scipy.linalg import eig, det
from numpy import argmax, array, vstack, diag, dot, abs, cumprod, sum, full, arccos, exp, real
from crysremesh.geom import Vector
from crysremesh.io import write_tecplot
from crysremesh.triangular_grid import Grid
from copy import deepcopy
from collections import deque

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

    def move_node(self, node, shift: Vector):
        """Move node for a given shift.

        Parameters
        ----------
            node : Node object
              node
            shift : Vector object
              shift
        """
        if self.fix_boundary_nodes:
            if node.fixed:
                pass
            else:
                node.move(shift)
        else:
            node.move(shift)

# --------------------------------------------------------------------------------------------------

    def write_grid_and_print_info(self, iteration):
        """Write info and save current grid.
        
        Parameters
        ----------
            iteration : int
              number of iteration
        """
        print('{}th iteration of {} smoothing'.format(iteration, self.__name__))
        write_tecplot(self.grid, '{}_smoothing_{}.dat'.format(self.__name__,
                                                              self.name_of_iteration(iteration)))

# --------------------------------------------------------------------------------------------------

    def apply_laplacians(self, laplacians):
        """Apply precomputed list of laplacians to the nodes.
        
        Parameters
        ----------
            laplacians : list of Vector objects
              shifts for nodes
        """
        assert len(self.grid.Nodes) == len(laplacians)
        for n, l in zip(self.grid.Nodes, laplacians):
            self.move_node(n, l)


# ==================================================================================================

class NullSpaceSmoothing(Smoothing):
    __name__ = "NullSpace"

    def __init__(self, grid: Grid, num_iterations=0, st=0.2, epsilon=1e-2, print_intermediate_steps=False, 
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
        N = node.faces[0].normal().coords_np_array()
        N = N.reshape((1, 3))
        for f in node.faces[1:]:
            normal = f.normal().coords_np_array()
            assert N.shape[1] == normal.shape[1], print(N.shape, normal.shape)
            N = vstack((N, normal))

        return N

# --------------------------------------------------------------------------------------------------

    def vector_to_avg_of_centroids(self, node):
        """Calculate the vector to the sum of centroids from a node.

        Parameters
        ----------
            node : Node object
              target node
        
        Returns
        -------
            Vector obj
        """
        dv = Vector()
        p = node.as_vector()
        weights = 0

        neighbours = node.faces

        assert len(neighbours) > 0
        assert isinstance(node.Id, int)

        for f in neighbours:
            c = f.centroid()
            c -= p
            weight = 1.0

            assert isinstance(c, Vector)

            c *= weight
            dv += c
            weights += weight

        dv /= weights
        return dv

# --------------------------------------------------------------------------------------------------

    def smoothing(self):
        """Perform smoothing."""
        for i in range(0, self.num_iterations):
            laplacians = []
            for n in self.grid.Nodes:
                m = len(n.faces)
                w = [f.area() for f in n.faces]
                N = self.create_matrix_of_normals(n)
                dv = self.vector_to_avg_of_centroids(n)

                assert N.shape == (m, 3)
                assert len(w) == m
                assert isinstance(dv, Vector)

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
                dv = dv.coords_np_array().reshape(3, 1)

                assert ns.shape == (3, 3 - k), 'Wrong eigenvectors shape'
                assert dv.shape == (3, 1), 'Wrong shape'

                if k < 3:
                    ttT = dot(ns, ns.T)
                    t = self.st * dot(ttT, dv)
                    assert t.shape == (3, 1), 'Wrong shift shape'
                    laplacian = Vector(real(t[0, 0]), real(t[1, 0]), real(t[2, 0]))
                else:
                    laplacian = Vector(0, 0, 0)
                laplacians.append(laplacian)

            Smoothing.apply_laplacians(self, laplacians)

            if self.print_intermediate_steps:
                Smoothing.write_grid_and_print_info(self, i)
