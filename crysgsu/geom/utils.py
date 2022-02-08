"""
Utils functions.
"""

from geom import Vect

# ==================================================================================================


def tetrahedron_volume(a: Vect,
                       b: Vect,
                       c: Vect,
                       d: Vect):
    """
    Tetrahedron volume.

    Parameters
    ----------
        a : Vector
            First vector
        b : Vector
            Second vector
        c : Vector
            Third vector
        d : Vector
            4-th vector

    Return
    ------
        Volume
    """

    return abs((a - d) * Vect.cross_product(b - d, c - d)) / 6.0

# --------------------------------------------------------------------------------------------------


def displaced_triangle_volume(a: Vect,
                              b: Vect,
                              c: Vect,
                              na: Vect,
                              nb: Vect,
                              nc: Vect):
    """
    Displaced triangle volume from 6 points.

    Parameters
    ----------
        a : Vector
            First vector
        b : Vector
            Second vector
        c : Vector
            Third vector
        na : Vector
            New position of the first vector
        nb : Vector
            New position of the second vector
        nc : Vector
            New position of the third vector

    Return
    ------
        Volume
    """

    return (tetrahedron_volume(a, b, c, nc) +
            tetrahedron_volume(a, b, nb, nc) +
            tetrahedron_volume(a, na, nb, nc))

# --------------------------------------------------------------------------------------------------


def prizmatoid_volume_coefs(v1: Vect,
                            v2: Vect,
                            v3: Vect,
                            n1: Vect,
                            n2: Vect,
                            n3: Vect,
                            n: Vect):
    """
    Prizmatoid volume coefficients (from Tong's article).

    Parameters
    ----------
        p1 : Vector
            First vector
        p2 : Vector
            Second vector
        p3 : Vector
            Third vector
        n1 : Vector
            Nodal normal of the firts vector
        n2 : Vector
            Nodal normal of the second vector
        n3 : Vector
            Nodal normal of the third vector
        n : Vector
            Face normal

    Return
    ------
        a, b, c coefficients
    """

    # p_ij, u_ij
    v21 = v2 - v1
    v31 = v3 - v1
    u1 = n1 / Vect.dot_product(n, n1)
    u2 = n2 / Vect.dot_product(n, n2)
    u3 = n3 / Vect.dot_product(n, n3)
    u21 = u2 - u1
    u31 = u3 - u1

    # a, b, c coefficients
    a = 0.5 * Vect.cross_product(v21, v31).mod()
    b = 0.25 * ((Vect.cross_product(v21, u31) + Vect.cross_product(u21, v31)) * n)
    c = (1.0 / 6.0) * (Vect.cross_product(u21, u31) * n)

    return a, b, c

# ==================================================================================================
