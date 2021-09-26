"""
GSU.
Geometry aspects of Grid Surface Unstructured.
"""

import math
import gsu

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

    def __sub__(self, v):
        """
        Sub.
        :param v: Vector.
        :return: Result.
        """

        return Vect(self.X - v.X, self.Y - v.Y, self.Z - v.Z)

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

        return sqrt(self.norm2())

    # ----------------------------------------------------------------------------------------------

    def dist_to(self, v):
        """
        Distance to.
        :param v: Vector.
        :return: Distance.
        """

        return (self - v).norm()

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

    def is_no_intersection_with_triangle_by_boxes(self, t):
        """
        Check if there is no intersection with another triangle if we use boxes.
        :param t: Triangle.
        :return: True - if self and t boxes does not intersect,
                 False - if self and t boxes intersect.
        """

        (left1, right1) = self.box_corners()
        (left2, right2) = t.box_corners()

        print(left1, right1, left2, right2)

        return (left2.X > right1.X) or (right2.X < left1.X) \
               or (left2.Y > right1.Y) or (right2.Y < left1.Y) \
               or (left2.Z > right1.Z) or (right2.Z < left1.Z)

    # ----------------------------------------------------------------------------------------------

    def intersection_with_triangle(self, t):
        """
        Calculate intersection with another triangle.
        :param t: Triangle.
        :return: Intersection information.
        """

        # Intersection between two triangles is described as array of points.
        #   - empty array - no intersection,
        #   - 1 point - one triangle touch another triangle with its node,
        #   - 2 points - intersection by segment,
        #   - 3 points - triangles lay in one plane,
        #   - 4 points - triangles lay in one plane.

        # First of all check intersection by boxes.
        if self.is_no_intersection_with_triangle_by_boxes(t):
            return []

        return []

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

    def import_from_gsu_grid(self, g):
        """
        Import from GSU grid.
        :param g: Grid.
        """

        self.Faces = []

        for f in g.Faces:
            a = Vect(f.Nodes[0].P[0], f.Nodes[0].P[1], f.Nodes[0].P[2])
            b = Vect(f.Nodes[1].P[0], f.Nodes[1].P[1], f.Nodes[1].P[2])
            c = Vect(f.Nodes[2].P[0], f.Nodes[2].P[1], f.Nodes[2].P[2])
            self.Faces.append(Face(Triangle(a, b, c)))

    # ----------------------------------------------------------------------------------------------

    def store(self, filename):
        """
        Store mesh.
        :param filename: Name of file.
        """

        with open(filename, 'w', newline='\n') as of:

            of.write('VARIABLES="X", "Y", "Z"\n')
            of.write('ZONE T="SINGLE" NODES={0} ELEMENTS={1} '
                     'DATAPACKING=BLOCK ZONETYPE=FETRIANGLE '
                     'VARLOCATION=([4-3]=CELLCENTERED)\n'.format(self.nodes_count(),
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

            for i in range(self.faces_count()):
                of.write('{0} {1} {2}\n'.format(3 * i + 1, 3 * i + 2, 3 * i + 3))

            of.close()

# ==================================================================================================


if __name__ == '__main__':
    """
    Main function.
    """

    pass

# ==================================================================================================
