"""
Utils functions.
"""

import re
import pathlib


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
