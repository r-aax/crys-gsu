"""
Utils functions.
"""

import re
import pathlib
import math
import numpy as np

# --------------------------------------------------------------------------------------------------


def norm2(v):
    """
    Norm^2 for vector.
    :param v: vector
    :return: norm^2
    """

    (vx, vy, vz) = v

    return vx * vx + vy * vy + vz * vz

# --------------------------------------------------------------------------------------------------


def norm(v):
    """
    Norm.
    :param v: vector
    :return: norm
    """

    return math.sqrt(norm2(v))

# --------------------------------------------------------------------------------------------------


def dist2(a, b):
    """
    Distance^2 between points.
    :param a: point
    :param b: point
    :return: distance^2
    """

    (ax, ay, az) = a
    (bx, by, bz) = b

    return norm2((ax - bx, ay - by, az - bz))

# --------------------------------------------------------------------------------------------------


def dist(a, b):
    """
    Distance between points.
    :param a: point
    :param b: point
    :return: distance
    """

    return math.sqrt(dist2(a, b))

# --------------------------------------------------------------------------------------------------


def a_kb(a, k, b):
    """
    Calculate a + k * b
    :param a: point
    :param k: coeff
    :param b: point
    :return: new point
    """

    (ax, ay, az) = a
    (bx, by, bz) = b

    return (ax + k * bx, ay + k * by, az + k * bz)

# --------------------------------------------------------------------------------------------------


def cross_product(a, b):
    """
    Cross product of two vectors.
    :param a: vector
    :param b: vector
    :return: result of cross product
    """

    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

# --------------------------------------------------------------------------------------------------


def normalized(v):
    """
    Calculate normalized vector.
    :param v: vector
    :return: normalized vecror
    """

    (vx, vy, vz) = v
    n = norm(v)

    return (vx / n, vy / n, vz / n)

# --------------------------------------------------------------------------------------------------


def is_triangle_and_segment_intersect(a, b, c, p, q):
    """
    Check if triangle and segment inntersect.
    :param a: point of triangle
    :param b: point of triangle
    :param c: point of triangle
    :param p: point of segment
    :param q: point of segment
    :return: True - if there is intersection, False - otherwise.
    """

    #
    # Point (x, y, z) inside triangle can be represented as
    # x = x_a + (x_b - x_a) * alfa + (x_c - x_a) * beta
    # y = y_a + (y_b - y_a) * alfa + (y_c - y_a) * beta
    # z = z_a + (z_b - z_a) * alfa + (z_c - z_a) * beta
    # ...
    # x = x_a + x_ba * alfa + x_ca * beta
    # y = y_a + y_ba * alfa + y_ca * beta
    # z = z_a + z_ba * alfa + z_ca * beta
    #
    # Point (x, y, z) on segment can be represented as
    # x = x_p + (x_q - x_p) * phi
    # y = y_p + (y_q - y_p) * phi
    # x = z_p + (z_q - z_p) * phi
    # ...
    # x = x_p + x_qp * phi
    # y = y_p + y_qp * phi
    # x = z_p + z_qp * phi
    #
    # So to find intersection we have to solve system
    # x_a + x_ba * alfa + x_ca * beta = x_p + x_qp * phi
    # y_a + y_ba * alfa + y_ca * beta = y_p + y_qp * phi
    # z_a + z_ba * alfa + z_ca * beta = z_p + z_qp * phi
    # ...
    # x_ba * alfa + x_ca * beta + (-x_qp) * phi = x_p - x_a
    # y_ba * alfa + y_ca * beta + (-y_qp) * phi = y_p - y_a
    # z_ba * alfa + z_ca * beta + (-z_qp) * phi = z_p - z_a
    # ...
    # x_ba * alfa + x_ca * beta + x_pq * phi = x_pa
    # y_ba * alfa + y_ca * beta + y_pq * phi = y_pa
    # z_ba * alfa + z_ca * beta + z_pq * phi = z_pa
    #
    # Matrix of this system can be written in the following view
    # [x_ba x_ca x_pq
    #  y_ba y_ca y_pq
    #  z_ba z_ca z_pq]

    (xa, ya, za) = a
    (xb, yb, zb) = b
    (xc, yc, zc) = c
    (xp, yp, zp) = p
    (xq, yq, zq) = q

    m = np.array([[xb - xa, xc - xa, xp - xq],
                  [yb - ya, yc - ya, yp - yq],
                  [zb - za, zc - za, zp - zq]])
    d = np.linalg.det(m)

    if abs(d) < 1e-6:
        return False

    im = np.linalg.inv(m)
    r = im.dot(np.array([xp - xa, yp - ya, zp - za]))
    alfa, beta, phi = r[0], r[1], r[2]

    return (alfa >= 0.0) and (beta >= 0.0) and (alfa + beta <= 1.0)\
           and (phi >= 0.0) and (phi <= 1.0)

# --------------------------------------------------------------------------------------------------


def flatten(ar):
    """
    Flatten array.
    :param ar: array
    :return: flat array
    """

    r = []

    for e in ar:
        if isinstance(e, list):
            r += flatten(e)
        else:
            r.append(e)

    return r

# --------------------------------------------------------------------------------------------------


def has_filename_timestamp(s):
    """
    Check if filename has timestamp.
    :param s: name of file without extension
    :return: True - if filename has timestamp, False - otherwise
    """

    return re.search(r'_\d\d\d\d\d\d\d\d\d\d\d\d$', s) is not None

# --------------------------------------------------------------------------------------------------


def get_filename_and_timestamp_pair(s):
    """
    Get tuple of base filename and timestamp.
    :param s: filename string
    :return: tuple
    """

    tsl = 12

    if has_filename_timestamp(s):
        return s[:-13], s[-12:]
    else:
        return s, '000000000000'

# --------------------------------------------------------------------------------------------------


def dummy_split_filename(s):
    """
    Dummy function of splitting filename into 4 pieces:
        - stem
        - mpi id (if it is digit) or None
        - timestamp (if it is digit) or None
        - extension
    :param s: filename
    :return: tuple if four elements (str, int|None, int|None, str)
    """

    s, m, t, e = s[:-23], s[-22:-17], s[-16:-4], s[-3:]

    if m.isdigit():
        md = int(m)
    else:
        md = None

    if t.isdigit():
        td = int(t)
    else:
        td = None

    return s, md, td, e

# --------------------------------------------------------------------------------------------------


def is_filename_correct_crys_cry_file(fn, stem):
    """
    Check filename for correct crys cry.
    :param fn: filename
    :param stem: stem for check
    :return: True - if file is correct crys cry filename, False - otherwise
    """

    s, m, t, e = dummy_split_filename(fn)

    return (s == stem) and (m is not None) and (t is not None) and (e == 'cry')

# --------------------------------------------------------------------------------------------------


def is_filename_correct_crys_txt_file(fn, stem):
    """
    Check filename for correct crys txt.
    :param fn: filename
    :param stem: stem for check
    :return: True - if file is correct crys txt filename, False - otherwise
    """

    s, m, t, e = dummy_split_filename(fn)

    return (s == stem) and (m is not None) and (t is not None) and (e == 'txt')

# --------------------------------------------------------------------------------------------------


def group_txt_files_by_timestamps(fs):
    """
    Group list of files [f1, f2, ..., fn] into set { timestamp : list of files }.
    :param fs: files
    :return: set of timestamp : list of files.
    """

    d = dict()

    for f in fs:
        tm = f[-16:-4]
        if tm in d:
            v = d.get(tm)
            v.append(f)
            d.update([(tm, v)])
        else:
            d.update([(tm, [f])])

    return d

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    # is_triangle_and_segment_intersect
    assert(is_triangle_and_segment_intersect((0.0, 0.0, 0.0),
                                             (1.0, 0.0, 0.0),
                                             (0.0, 0.0, 1.0),
                                             (0.25, 1.0, 0.25),
                                             (0.25, -1.0, 0.25)))
    assert(not is_triangle_and_segment_intersect((0.0, 0.0, 0.0),
                                                 (1.0, 0.0, 0.0),
                                                 (0.0, 0.0, 1.0),
                                                 (1.0, 1.0, 1.0),
                                                 (1.0, -1.0, 1.0)))

    # flatten
    assert(flatten([]) == [])
    assert(flatten([1]) == [1])
    assert(flatten([1, 2]) == [1, 2])
    assert(flatten([1, [2, 3]]) == [1, 2, 3])
    assert(flatten([1, [2, [3], 4], 5]) == [1, 2, 3, 4, 5])

    # has_filename_timestamp
    assert(not has_filename_timestamp(''))
    assert(not has_filename_timestamp('short'))
    assert(not has_filename_timestamp('name_000000aaa000'))
    assert(not has_filename_timestamp('nameX000111222333'))
    assert(has_filename_timestamp('name_000000000000'))
    assert(has_filename_timestamp('name_111222333444'))
    assert(not has_filename_timestamp('name_111222333444_r'))

    # get_filename_and_timestamp_pair
    assert(get_filename_and_timestamp_pair('name') == ('name', '000000000000'))
    assert(get_filename_and_timestamp_pair('name_000111222333') == ('name', '000111222333'))

    # dummy_split_filename
    assert(dummy_split_filename('bunny_00001_111222333444.cry') == \
           ('bunny', 1, 111222333444, 'cry'))

    # is_filename_correct_crys_{cry/txt}_file
    assert(is_filename_correct_crys_cry_file('bunny_00004_444333222111.cry', 'bunny'))
    assert(is_filename_correct_crys_txt_file('bunny_00004_444333222111.txt', 'bunny'))

# --------------------------------------------------------------------------------------------------
