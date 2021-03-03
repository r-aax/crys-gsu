"""
Utils functions.
"""

import re

# ----------------------------------------------------------------------------------------------------------------------


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

# ----------------------------------------------------------------------------------------------------------------------


def has_filename_timestamp(s):
    """
    Check is filename has timestamp.
    :param s: name of file without extension
    :return: True - if filename has timestamp, False - otherwise
    """

    return re.search(r'_\d\d\d\d\d\d\d\d\d\d\d\d$', s) is not None

# ----------------------------------------------------------------------------------------------------------------------


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

# ----------------------------------------------------------------------------------------------------------------------


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

# ----------------------------------------------------------------------------------------------------------------------
