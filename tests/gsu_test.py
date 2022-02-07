import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../crysgsu/')
import geom
import gsu


class TestUtils(unittest.TestCase):

    def setUp(self):
        self.g = gsu.Grid()
        self.g.load('cases/grids/wing_1.dat')

    def test_atomic_grid_transformation(self):

        def check_GloId(list_old_id, list_new_id, number_new_items):
            len_old = len(list_old_id)
            len_new = len(list_new_id)
            if len_new == len_old + number_new_items:
                list_new_id.sort()
                for n, i in enumerate(list_new_id):
                    if n != i:
                        return False
                return True
            return False

        # divide_face
        list_old_id = [f.GloId for f in self.g.Faces]
        f = self.g.Faces[0]
        self.g.divide_face(f, f.get_triangle().centroid())
        list_new_id =[f.GloId for f in self.g.Faces]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 3 - 1)) # плюс 3 face и минус 1 face
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Edges]
        f = self.g.Faces[0]
        self.g.divide_face(f, f.get_triangle().centroid())
        list_new_id = [f.GloId for f in self.g.Edges]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 3)) # плюс 3 edge
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Nodes]
        f = self.g.Faces[0]
        self.g.divide_face(f, f.get_triangle().centroid())
        list_new_id = [f.GloId for f in self.g.Nodes]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 1)) # плюс 1 node
        self.setUp()

        # collapse_face
        list_old_id = [f.GloId for f in self.g.Faces]
        f = self.g.Faces[60]
        self.g.collapse_face(f)
        list_new_id = [f.GloId for f in self.g.Faces]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 0))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Edges]
        f = self.g.Faces[60]
        self.g.collapse_face(f)
        list_new_id = [f.GloId for f in self.g.Edges]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 0))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Nodes]
        f = self.g.Faces[60]
        self.g.collapse_face(f)
        list_new_id = [f.GloId for f in self.g.Nodes]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 0))
        self.setUp()

        # cut_edge
        list_old_id = [f.GloId for f in self.g.Faces]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.cut_edge(self.g.Faces[n].Edges[0], self.g.Faces[n].get_triangle().centroid())
            n += 1
        list_new_id = [f.GloId for f in self.g.Faces]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 2))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Edges]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.cut_edge(self.g.Faces[n].Edges[0], self.g.Faces[n].get_triangle().centroid())
            n += 1
        list_new_id = [f.GloId for f in self.g.Edges]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 3))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Nodes]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.cut_edge(self.g.Faces[n].Edges[0], self.g.Faces[n].get_triangle().centroid())
            n += 1
        list_new_id = [f.GloId for f in self.g.Nodes]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 1))
        self.setUp()

        # collapse_edge
        list_old_id = [f.GloId for f in self.g.Faces]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.collapse_edge(self.g.Faces[n].Edges[0])
            n += 1
        list_new_id = [f.GloId for f in self.g.Faces]
        self.assertTrue(check_GloId(list_old_id, list_new_id, -2))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Edges]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.collapse_edge(self.g.Faces[n].Edges[0])
            n += 1
        list_new_id = [f.GloId for f in self.g.Edges]
        self.assertTrue(check_GloId(list_old_id, list_new_id, -3))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Nodes]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.collapse_edge(self.g.Faces[n].Edges[0])
            n += 1
        list_new_id = [f.GloId for f in self.g.Nodes]
        self.assertTrue(check_GloId(list_old_id, list_new_id, -1))
        self.setUp()

        # cut_single_edge
        list_old_id = [f.GloId for f in self.g.Faces]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 1
            if not w:
                self.g.cut_single_edge(self.g.Faces[n].Edges[0], self.g.Faces[n].get_triangle().centroid())
            n += 1
        list_new_id = [f.GloId for f in self.g.Faces]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 2-1))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Edges]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 1
            if not w:
                self.g.cut_single_edge(self.g.Faces[n].Edges[0], self.g.Faces[n].get_triangle().centroid())
            n += 1
        list_new_id = [f.GloId for f in self.g.Edges]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 3-1))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Nodes]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 1
            if not w:
                self.g.cut_single_edge(self.g.Faces[n].Edges[0], self.g.Faces[n].get_triangle().centroid())
            n += 1
        list_new_id = [f.GloId for f in self.g.Nodes]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 1))
        self.setUp()

        # cut_edge_with_two_nodes
        # 1 face in cut edge
        list_old_id = [f.GloId for f in self.g.Faces]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 1
            if not w:
                self.g.cut_edge_with_two_nodes(self.g.Faces[n].Edges[0],
                                          self.g.Faces[n].get_triangle().centroid() - geom.Vect(0.01, 0.01, 0.01),
                                          self.g.Faces[n].get_triangle().centroid() + geom.Vect(0.01, 0.01, 0.01))
            n += 1
        list_new_id = [f.GloId for f in self.g.Faces]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 3-1))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Edges]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 1
            if not w:
                self.g.cut_edge_with_two_nodes(self.g.Faces[n].Edges[0],
                                               self.g.Faces[n].get_triangle().centroid() - geom.Vect(0.01, 0.01, 0.01),
                                               self.g.Faces[n].get_triangle().centroid() + geom.Vect(0.01, 0.01, 0.01))
            n += 1
        list_new_id = [f.GloId for f in self.g.Edges]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 5 - 1))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Nodes]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 1
            if not w:
                self.g.cut_edge_with_two_nodes(self.g.Faces[n].Edges[0],
                                               self.g.Faces[n].get_triangle().centroid() - geom.Vect(0.01, 0.01, 0.01),
                                               self.g.Faces[n].get_triangle().centroid() + geom.Vect(0.01, 0.01, 0.01))
            n += 1
        list_new_id = [f.GloId for f in self.g.Nodes]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 2))
        self.setUp()

        # 2 faces in cut edge
        list_old_id = [f.GloId for f in self.g.Faces]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.cut_edge_with_two_nodes(self.g.Faces[n].Edges[0],
                                               self.g.Faces[n].get_triangle().centroid() - geom.Vect(0.01, 0.01, 0.01),
                                               self.g.Faces[n].get_triangle().centroid() + geom.Vect(0.01, 0.01, 0.01))
            n += 1
        list_new_id = [f.GloId for f in self.g.Faces]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 6 - 2))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Edges]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.cut_edge_with_two_nodes(self.g.Faces[n].Edges[0],
                                               self.g.Faces[n].get_triangle().centroid() - geom.Vect(0.01, 0.01, 0.01),
                                               self.g.Faces[n].get_triangle().centroid() + geom.Vect(0.01, 0.01, 0.01))
            n += 1
        list_new_id = [f.GloId for f in self.g.Edges]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 7 - 1))
        self.setUp()

        list_old_id = [f.GloId for f in self.g.Nodes]
        n = 0
        len(self.g.Faces[n].Edges[0].Faces)
        w = True
        while w:
            w = not len(self.g.Faces[n].Edges[0].Faces) == 2
            if not w:
                self.g.cut_edge_with_two_nodes(self.g.Faces[n].Edges[0],
                                               self.g.Faces[n].get_triangle().centroid() - geom.Vect(0.01, 0.01, 0.01),
                                               self.g.Faces[n].get_triangle().centroid() + geom.Vect(0.01, 0.01, 0.01))
            n += 1
        list_new_id = [f.GloId for f in self.g.Nodes]
        self.assertTrue(check_GloId(list_old_id, list_new_id, 2))
        self.setUp()


t = TestUtils()
t.setUp()
t.test_atomic_grid_transformation()