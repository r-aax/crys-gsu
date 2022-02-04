"""
GSU.
Geometry aspects of Grid Surface Unstructured.
"""

import math
from gsu import gsu
import numpy as np
import scipy as sp
import scipy.linalg as lalg
from geom.vect import Vect
from geom.segment import Segment

# ==================================================================================================

eps = 1.0e-10

# ==================================================================================================


class Triangle:
    """
    Triangle.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, a, b, c):
        """
        Constructor.

        :param a: A point.
        :param b: B point.
        :param c: C point.
        """

        self.Points = [a, b, c]

    # ----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        :return: String.
        """

        return 'T: {0}, {1}, {2}'.format(self.Points[0], self.Points[1], self.Points[2])

    # ----------------------------------------------------------------------------------------------

    def a(self):
        """
        A point.

        :return: A point.
        """

        return self.Points[0]

    # ----------------------------------------------------------------------------------------------

    def b(self):
        """
        B point.

        :return: B point.
        """

        return self.Points[1]

    # ----------------------------------------------------------------------------------------------

    def c(self):
        """
        C point.

        :return: C point.
        """

        return self.Points[2]

    # ----------------------------------------------------------------------------------------------

    def center(self):
        """
        Center point.

        :return: Center.
        """

        return (self.a() + self.b() + self.c()) / 3.0

    # ----------------------------------------------------------------------------------------------

    def ab(self):
        """
        AB segment.

        :return: AB segment.
        """

        return Segment(self.a(), self.b())

    # ----------------------------------------------------------------------------------------------

    def bc(self):
        """
        BC segment.

        :return: BC segment.
        """

        return Segment(self.b(), self.c())

    # ----------------------------------------------------------------------------------------------

    def ac(self):
        """
        AC segment.

        :return: AC segment.
        """

        return Segment(self.a(), self.c())

    # ----------------------------------------------------------------------------------------------

    def is_my_vertex(self, p):
        """
        Check if point is vertex.

        :param p: Point.
        :return:
            True - if point is my vertex,
            False - otherwise.
        """

        return self.a().is_near(p, eps) or self.b().is_near(p, eps) or self.c().is_near(p, eps)

    # ----------------------------------------------------------------------------------------------

    def is_my_segment(self, s):
        """
        Check if segment is edge of this triangle.

        :param s: Segment.
        :return:
            True - if segment is an edge of triangle,
            False - otherwise.
        """

        return s.is_eq(self.ab()) or s.is_eq(self.bc()) or s.is_eq(self.ac())

    # ----------------------------------------------------------------------------------------------

    def is_eq(self, t):
        """
        Check if eq to another triangle.

        :param t: Triangle.
        :return:
            True - if near to another triangle,
            False - otherwise.
        """

        a, b, c, ab, bc, ac = self.a(), self.b(), self.c(), self.ab(), self.bc(), self.ac()
        ta, tb, tc, tab, tbc, tac = t.a(), t.b(), t.c(), t.ab(), t.bc(), t.ac()

        return (ab.is_eq(tab) and c.is_eq(tc)) \
               or (ab.is_eq(tbc) and c.is_eq(ta)) \
               or (ab.is_eq(tac) and c.is_eq(tb))

    # ----------------------------------------------------------------------------------------------

    def normal(self):
        """
        Normal.

        :return: Normal.
        """

        a, b, c = self.a(), self.b(), self.c()

        return (b - a).cross_product(c - a).orth()

    # ----------------------------------------------------------------------------------------------

    def is_parallel_to_triangle(self, t):
        """
        Check if parallel to triangle.

        :param t: Triangle.
        :return:
            True - if parallel to triangle,
            False - otherwise.
        """

        n1 = self.normal()
        n2 = self.normal()

        return n1.angle_cos(n2) < eps

    # ----------------------------------------------------------------------------------------------

    def box_corner_ldb(self):
        """
        Box corner Left-Down-Back.

        :return: Box corner.
        """

        a, b, c, = self.a(), self.b(), self.c()

        return Vect(min(a.X, b.X, c.X), min(a.Y, b.Y, c.Y), min(a.Z, b.Z, c.Z))

    # ----------------------------------------------------------------------------------------------

    def box_corner_ruf(self):
        """
        Box corner Right-Up-Front.

        :return: Box corner.
        """

        a, b, c, = self.a(), self.b(), self.c()

        return Vect(max(a.X, b.X, c.X), max(a.Y, b.Y, c.Y), max(a.Z, b.Z, c.Z))

    # ----------------------------------------------------------------------------------------------

    def box_corners(self):
        """
        Box corners points.

        :return: Corner points.
        """

        return (self.box_corner_ldb(), self.box_corner_ruf())

    # ----------------------------------------------------------------------------------------------

    def move(self, v):
        """
        Move with vector.

        :param v: Vector.
        """

        a, b, c = self.a(), self.b(), self.c()
        self.Points = [a + v, b + v, c + v]

    # ----------------------------------------------------------------------------------------------

    def is_no_intersection_with_triangle_by_boxes(self, t):
        """
        Check if there is no intersection with another triangle if we use boxes.

        :param t: Triangle.
        :return:
            True - if self and t boxes does not intersect,
            False - if self and t boxes intersect.
        """

        (left1, right1) = self.box_corners()
        (left2, right2) = t.box_corners()

        return (left2.X > right1.X) or (right2.X < left1.X) \
               or (left2.Y > right1.Y) or (right2.Y < left1.Y) \
               or (left2.Z > right1.Z) or (right2.Z < left1.Z)

    # ----------------------------------------------------------------------------------------------

    def has_common_edge_with_triangle(self, t):
        """
        Check if has common edge with another trianngle.

        :param t: Triangle.
        :return:
            True - if has common edge with another triangle,
            False - otherwise.
        """

        ab, bc, ac = self.ab(), self.bc(), self.ac()
        tab, tbc, tac = t.ab(), t.bc(), t.ac()

        return ab.is_eq(tab) or ab.is_eq(tbc) or ab.is_eq(tac) \
               or bc.is_eq(tab) or bc.is_eq(tbc) or bc.is_eq(tac) \
               or ac.is_eq(tab) or ac.is_eq(tbc) or ac.is_eq(tac)

    # ----------------------------------------------------------------------------------------------

    def area(self):
        """
        Area.

        :return: Area.
        """

        ab = self.b() - self.a()
        ac = self.c() - self.a()

        return 0.5 * abs(Vect.cross_product(ab, ac).mod())

    # ----------------------------------------------------------------------------------------------

    def is_point_inside(self, p):
        """
        Triangle has point inside it.

        :param p: Point.
        :return:
            True - if point is inside triangle,
            False - otherwise.
        """

        a, b, c = self.a(), self.b(), self.c()
        ar1 = self.area()
        t1 = Triangle(a, b, p)
        t2 = Triangle(b, c, p)
        t3 = Triangle(a, c, p)
        ar2 = t1.area() + t2.area() + t3.area()

        return abs(ar1 - ar2) < eps

    # ----------------------------------------------------------------------------------------------

    def intersect_with_segment(self, s):
        """
        Intersection with segment.

        :param s: Segment.
        """

        # Triangle ABC equation.
        #   x = xa + (xb - xa) * alf + (xc - xa) * bet
        #   y = ya + (yb - ya) * alf + (yc - ya) * bet
        #   z = za + (zb - za) * alf + (zc - za) * bet
        #   alf >= 0.0
        #   bet >= 0.0
        #   alf + bet <= 1.0
        # Segment DE equation.
        #   x = xd + (xe - xd) * gam
        #   y = yd + (ye - yd) * gam
        #   z = zd + (ze - zd) * gam
        #   gam >= 0.0
        #   gam <= 1.0
        # Rewrite equations.
        #   xa + (xb - xa) * alf + (xc - xa) * bet = xd + (xe - xd) * gam
        #   ya + (yb - ya) * alf + (yc - ya) * bet = yd + (ye - yd) * gam
        #   za + (zb - za) * alf + (zc - za) * bet = zd + (ze - zd) * gam
        # ...
        #   (xb - xa) * alf + (xc - xa) * bet + (xd - xe) * gam = xd - xa
        #   (yb - ya) * alf + (yc - ya) * bet + (yd - ye) * gam = yd - ya
        #   (zb - za) * alf + (zc - za) * bet + (zd - ze) * gam = zd - za

        a, b, c = self.a(), self.b(), self.c()
        d, e = s.a(), s.b()

        matrix_a = np.array([[b.X - a.X, c.X - a.X, d.X - e.X],
                             [b.Y - a.Y, c.Y - a.Y, d.Y - e.Y],
                             [b.Z - a.Z, c.Z - a.Z, d.Z - e.Z]])
        vector_b = np.array([d.X - a.X, d.Y - a.Y, d.Z - a.Z])

        # matrix_a * [alf, bet, gam] = vector_b

        if abs(lalg.det(matrix_a)) > eps:
            r = lalg.solve(matrix_a, vector_b)
            alf, bet, gam = r[0], r[1], r[2]
            is_alf_bet = (alf >= 0.0) and (bet >= 0.0) and (alf + bet <= 1.0)
            is_gam = (gam >= 0.0) and (gam <= 1.0)
            if is_alf_bet and is_gam:
                # print('POINT', d + (e - d) * gam)
                return d + (e - d) * gam

        return None

# ==================================================================================================


class Face:
    """
    Face.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, t):
        """
        Constructor.

        :param t: Triangle.
        """

        self.T = t

        # Mark.
        self.M = 0

    # ----------------------------------------------------------------------------------------------

    def add_point_to_sp(self, p):
        """
        Add point to SP.

        :param p: Point.
        """

        # Do not add None point.
        if p is None:
            return

        # Do not add vertex point.
        if self.T.is_my_vertex(p):
            return

        # Double point check.
        for pi in self.SP:
            if pi.is_near(p, eps):
                return

        # Eventually add it.
        self.SP.append(p)

    # ----------------------------------------------------------------------------------------------

    def shred(self):
        """
        Shred.
        """

        fs = [self]

        for p in self.SP:
            nfs = []
            for f in fs:
                t = f.T
                if f.T.is_point_inside(p):
                    a, b, c = t.a(), t.b(), t.c()
                    t1 = Triangle(a, b, p)
                    t2 = Triangle(b, c, p)
                    t3 = Triangle(c, a, p)
                    if t1.area() > eps:
                        nfs.append(Face(t1))
                    if t2.area() > eps:
                        nfs.append(Face(t2))
                    if t3.area() > eps:
                        nfs.append(Face(t3))
                else:
                    nfs.append(f)
            fs = nfs

        for f in fs:
            f.M = 1.0

        return fs

    # ----------------------------------------------------------------------------------------------

    def edge_neighbours(self, s):
        """
        Neighbours over the segment.

        :param s: Segment.
        :return: List of neighbours.
        """

        return [l for l in self.L if l.T.is_my_segment(s)]

# ==================================================================================================


class Mesh:
    """
    Mesh.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.Faces = []

    # ----------------------------------------------------------------------------------------------

    def faces_count(self):
        """
        Get faces count.

        :return: Faces count.
        """

        return len(self.Faces)

    # ----------------------------------------------------------------------------------------------

    def nodes_count(self):
        """
        Get nodes count.

        :return: Nodes count.
        """

        return 3 * self.faces_count()

    # ----------------------------------------------------------------------------------------------

    def box_corner_ldb(self):
        """
        Box corner Left-Down-Back.

        :return: Box corner.
        """

        ps = [f.T.box_corner_ldb() for f in self.Faces]

        return Vect(min([p.X for p in ps]),
                    min([p.Y for p in ps]),
                    min([p.Z for p in ps]))

    # ----------------------------------------------------------------------------------------------

    def box_corner_ruf(self):
        """
        Box corner Right-Up-Front.

        :return: Box corner.
        """

        ps = [f.T.box_corner_ruf() for f in self.Faces]

        return Vect(max([p.X for p in ps]),
                    max([p.Y for p in ps]),
                    max([p.Z for p in ps]))

    # ----------------------------------------------------------------------------------------------

    def box_corners(self):
        """
        Box corners points.

        :return: Corner points.
        """

        return (self.box_corner_ldb(), self.box_corner_ruf())

    # ----------------------------------------------------------------------------------------------

    def add_from_gsu_grid(self, g):
        """
        Import from GSU grid.

        :param g: Grid.
        """

        for f in g.Faces:
            a = Vect(f.Nodes[0].P[0], f.Nodes[0].P[1], f.Nodes[0].P[2])
            b = Vect(f.Nodes[1].P[0], f.Nodes[1].P[1], f.Nodes[1].P[2])
            c = Vect(f.Nodes[2].P[0], f.Nodes[2].P[1], f.Nodes[2].P[2])
            self.Faces.append(Face(Triangle(a, b, c)))

    # ----------------------------------------------------------------------------------------------

    def add_vertical_square(self, x1, y1, x2, y2):
        """
        Add vertical square, divided by two triangles.

        :param x1: First X.
        :param y1: First Y.
        :param x2: Second X.
        :param y2: Second Y.
        """

        self.Faces.append(Face(Triangle(Vect(x1, y1, 0.0),
                                        Vect(x2, y2, 0.0),
                                        Vect(x2, y2, 1.0))))
        self.Faces.append(Face(Triangle(Vect(x1, y1, 0.0),
                                        Vect(x2, y2, 1.0),
                                        Vect(x1, y1, 1.0))))

    # ----------------------------------------------------------------------------------------------

    def move(self, v):
        """
        Move with vector.

        :param v: Vector.
        """

        for f in self.Faces:
            f.T.move(v)

    # ----------------------------------------------------------------------------------------------

    def filter(self, fun):
        """
        Filter faces.

        :param fun: Function.
        """

        self.Faces = [f for f in self.Faces if fun(f)]

    # ----------------------------------------------------------------------------------------------

    def load(self, filename):
        """
        Load mesh.

        :param filename: Name of file.
        """

        with open(filename, 'r') as f:

            # Zone head.
            hl = f.readline()
            hl = f.readline()
            hls = hl.split()
            faces_count = int(hls[3].split('=')[1])

            # Coordinates.
            xs = [float(xsi) for xsi in f.readline().split()]
            ys = [float(ysi) for ysi in f.readline().split()]
            zs = [float(zsi) for zsi in f.readline().split()]

            # Create faces.
            for i in range(faces_count):
                a = Vect(xs[3 * i], ys[3 * i], zs[3 * i])
                b = Vect(xs[3 * i + 1], ys[3 * i + 1], zs[3 * i + 1])
                c = Vect(xs[3 * i + 2], ys[3 * i + 2], zs[3 * i + 2])
                self.Faces.append(Face(Triangle(a, b, c)))

            # Load marks.
            ms = [float(msi) for msi in f.readline().split()]
            for i in range(faces_count):
                self.Faces[i].M = ms[i]

            # We do not need to load links between faces and nodes.

            f.close()

    # ----------------------------------------------------------------------------------------------

    def store(self, filename):
        """
        Store mesh.

        :param filename: Name of file.
        """

        with open(filename, 'w', newline='\n') as of:

            of.write('VARIABLES="X", "Y", "Z", "M"\n')
            of.write('ZONE T="SINGLE" NODES={0} ELEMENTS={1} '
                     'DATAPACKING=BLOCK ZONETYPE=FETRIANGLE '
                     'VARLOCATION=([4-4]=CELLCENTERED)\n'.format(self.nodes_count(),
                                                                 self.faces_count()))
            of.write(' '.join(['{0} {1} {2}'.format(f.T.a().X, f.T.b().X, f.T.c().X)
                                 for f in self.Faces]))
            of.write('\n')
            of.write(' '.join(['{0} {1} {2}'.format(f.T.a().Y, f.T.b().Y, f.T.c().Y)
                                 for f in self.Faces]))
            of.write('\n')
            of.write(' '.join(['{0} {1} {2}'.format(f.T.a().Z, f.T.b().Z, f.T.c().Z)
                                 for f in self.Faces]))
            of.write('\n')
            of.write(' '.join([str(f.M) for f in self.Faces]))
            of.write('\n')

            for i in range(self.faces_count()):
                of.write('{0} {1} {2}\n'.format(3 * i + 1, 3 * i + 2, 3 * i + 3))

            of.close()

    # ----------------------------------------------------------------------------------------------

    def shred(self):
        """
        Shred mesh.
        """

        fc = self.faces_count()

        for f in self.Faces:
            f.SP = []

        for i in range(fc):

            fi = self.Faces[i]
            print('processing {0} of {1}'.format(i, fc))

            for j in range(i + 1, fc):

                fj = self.Faces[j]

                # Check boxes.
                if fi.T.is_no_intersection_with_triangle_by_boxes(fj.T):
                    continue

                # Intersections.
                ps = [fi.T.intersect_with_segment(fj.T.ab()),
                      fi.T.intersect_with_segment(fj.T.bc()),
                      fi.T.intersect_with_segment(fj.T.ac()),
                      fj.T.intersect_with_segment(fi.T.ab()),
                      fj.T.intersect_with_segment(fi.T.bc()),
                      fj.T.intersect_with_segment(fi.T.ac())]

                for p in ps:
                    fi.add_point_to_sp(p)
                    fj.add_point_to_sp(p)

        # Array for new faces.
        new_faces = []

        # Ready to shred.
        for f in self.Faces:
            pc = len(f.SP)
            # print(f.T, pc)
            if pc == 0:
                pass
            else:
                f.M = -1
                new_faces = new_faces + f.shred()

        self.Faces = self.Faces + new_faces
        self.filter(lambda f: f.M > -1)

    # ----------------------------------------------------------------------------------------------

    def wrap(self):
        """
        Wrap.
        """

        # Prepare edges.
        for f in self.Faces:
            f.L = []
            f.M = 0

        fc = self.faces_count()

        # Imitate edges.
        for i in range(fc):
            print('process {0} of {1}'.format(i, fc))
            fi = self.Faces[i]
            for j in range(i + 1, fc):
                fj = self.Faces[j]
                if fi.T.has_common_edge_with_triangle(fj.T):
                    fi.L.append(fj)
                    fj.L.append(fi)
        for f in self.Faces:
            abl = f.edge_neighbours(f.T.ab())
            bcl = f.edge_neighbours(f.T.bc())
            acl = f.edge_neighbours(f.T.ac())
            if (len(abl) > 1) or (len(bcl) > 1) or (len(acl) > 1):
                f.M = 1

        # TODO: initial nodes.
        # stack = []
        for f in stack:
            f.M = 1
        while len(stack) > 0:
            f = stack.pop(0)
            for ff in f.L:
                if ff.M == 0:
                    ff.M = 1
                    stack.append(ff)

        self.filter(lambda f: f.M == 0)

    # ----------------------------------------------------------------------------------------------

    def select_neighbour(self, f, l):
        """
        Select from list.

        :param f: Face.
        :param l: List.
        """

        f.M = 1.0

# ==================================================================================================


if __name__ == '__main__':
    """
    Main function.
    """

    pass

# ==================================================================================================
