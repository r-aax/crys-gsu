"""
GSU.
Geometry aspects of Grid Surface Unstructured.
"""

import math
import gsu

# ==================================================================================================

eps = 1.0e-10

# ==================================================================================================


class Vect:
    """
    Vect.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, x=0.0, y=0.0, z=0.0):
        """
        Constructor.
        :param x: X coord.
        :param y: Y coord.
        :param z: Z coord.
        """

        self.X = x
        self.Y = y
        self.Z = z

    # ----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.
        :return: String.
        """

        return '({0}, {1}, {2})'.format(self.X, self.Y, self.Z)

    # ----------------------------------------------------------------------------------------------

    def __add__(self, v):
        """
        Add.
        :param v: Vector.
        :return: Result.
        """

        return Vect(self.X + v.X, self.Y + v.Y, self.Z + v.Z)

    # ----------------------------------------------------------------------------------------------

    def __sub__(self, v):
        """
        Sub.
        :param v: Vector.
        :return: Result.
        """

        return Vect(self.X - v.X, self.Y - v.Y, self.Z - v.Z)

    # ----------------------------------------------------------------------------------------------

    def __mul__(self, k):
        """
        Mul.
        :param k: Value.
        :return: Result.
        """

        return Vect(self.X * k, self.Y * k, self.Z * k)

    # ----------------------------------------------------------------------------------------------

    def __truediv__(self, k):
        """
        Div.
        :param k: Value.
        :return: Result.
        """

        return self * (1.0 / k)

    # ----------------------------------------------------------------------------------------------

    def norm2(self):
        """
        Square norm.
        :return: Square norm.
        """

        return self.X * self.X + self.Y * self.Y + self.Z * self.Z

    # ----------------------------------------------------------------------------------------------

    def norm(self):
        """
        Norm.
        :return: Norm.
        """

        return math.sqrt(self.norm2())

    # ----------------------------------------------------------------------------------------------

    def orth(self):
        """
        Orth of the vector.
        :return: Orth.
        """

        n = self.norm()

        return Vect(self.X / n, self.Y / n, self.Z / n)

    # ----------------------------------------------------------------------------------------------

    def dist_to(self, v):
        """
        Distance to.
        :param v: Vector.
        :return: Distance.
        """

        return (self - v).norm()

    # ----------------------------------------------------------------------------------------------

    def is_eq(self, v):
        """
        Check if eq to another vector.
        :param v: Vector.
        :return: True - if near to another vector,
                 False - otherwise.
        """

        return abs(self.X - v.X) + abs(self.Y - v.Y) + abs(self.Z - v.Z) < eps

    # ----------------------------------------------------------------------------------------------

    def dot_product(self, v):
        """
        Dot product.
        :param v: Vector.
        :return: Dot product.
        """

        return self.X * v.X + self.Y * v.Y + self.Z * v.Z

    # ----------------------------------------------------------------------------------------------

    def angle_cos(self, v):
        """
        Cosine of angle between vectors.
        :param v: Vector.
        :return: Angle cosine.
        """

        return self.dot_product(v) / (self.norm() * v.norm())

    # ----------------------------------------------------------------------------------------------

    def cross_product(self, v):
        """
        Cross product.
        :param v: Vector.
        :return: Cross product.
        """

        return Vect(self.Y * v.Z - self.Z * v.Y,
                    self.Z * v.X - self.X * v.Z,
                    self.X * v.Y - self.Y * v.X)

# ==================================================================================================


class Segment:
    """
    Segment.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, a, b):
        """
        Constructor.
        :param a: A point.
        :param b: B point.
        """

        self.Points = [a, b]

    # ----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.
        :return: String.
        """

        return 'S: {0} - {1}'.format(self.a(), self.b())

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

    def is_eq(self, s):
        """
        Check if eq to another segment.
        :param s: Segment.
        :return: True - if near to another segment,
                 False - otherwise.
        """

        a, b = self.a(), self.b()
        sa, sb = s.a(), s.b()

        return (a.is_eq(sa) and b.is_eq(sb)) \
               or (a.is_eq(sb) and b.is_eq(sa))

    # ----------------------------------------------------------------------------------------------

    def len(self):
        """
        Length.
        :return: Length.
        """

        return (self.a() - self.b()).norm()

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

    def is_eq(self, t):
        """
        Check if eq to another triangle.
        :param t: Triangle.
        :return: True - if near to another triangle,
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
        :return: True - if parallel to triangle,
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
        :return: True - if self and t boxes does not intersect,
                 False - if self and t boxes intersect.
        """

        (left1, right1) = self.box_corners()
        (left2, right2) = t.box_corners()

        return (left2.X > right1.X) or (right2.X < left1.X) \
               or (left2.Y > right1.Y) or (right2.Y < left1.Y) \
               or (left2.Z > right1.Z) or (right2.Z < left1.Z)

    # ----------------------------------------------------------------------------------------------

    def is_has_common_edge_with_triangle(self, t):
        """
        Check if has common edge with another trianngle.
        :param t: Triangle.
        :return: True - if has common edge with another triangle,
                 False - otherwise.
        """

        ab, bc, ac = self.ab(), self.bc(), self.ac()
        tab, tbc, tac = t.ab(), t.bc(), t.ac()

        return ab.is_eq(tab) or ab.is_eq(tbc) or ab.is_eq(tac) \
               or bc.is_eq(tab) or bc.is_eq(tbc) or bc.is_eq(tac) \
               or ac.is_eq(tab) or ac.is_eq(tbc) or ac.is_eq(tac)

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

# ==================================================================================================


if __name__ == '__main__':
    """
    Main function.
    """

    # No intersection by boxes.
    t1 = Triangle(Vect(0.0, 0.0, 0.0), Vect(1.0, 0.0, 0.0), Vect(0.0, 1.0, 0.0))
    t2 = Triangle(Vect(6.0, 0.0, 0.0), Vect(5.0, 1.0, 0.0), Vect(5.0, 0.0, 0.0))
    assert t1.is_no_intersection_with_triangle_by_boxes(t2)

    # Exact equal.
    t1 = Triangle(Vect(0.0, 0.0, 0.0), Vect(1.0, 0.0, 0.0), Vect(0.0, 1.0, 0.0))
    t2 = Triangle(Vect(1.0, 0.0, 0.0), Vect(0.0, 1.0, 0.0), Vect(0.0, 0.0, 0.0))
    assert t1.is_eq(t2)

    # Common edge.
    t1 = Triangle(Vect(0.0, 0.0, 0.0), Vect(1.0, 0.0, 0.0), Vect(0.0, 1.0, 0.0))
    t2 = Triangle(Vect(1.0, 0.0, 0.0), Vect(0.0, 1.0, 0.0), Vect(1.0, 1.0, 0.0))
    assert t1.is_has_common_edge_with_triangle(t2)

# ==================================================================================================
