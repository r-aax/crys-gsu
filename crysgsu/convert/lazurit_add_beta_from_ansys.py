#!/usr/bin/env python3
"""
Добавление поля beta из сетки ansys.
Во время этого действия чтение и запись выполняеся с помощью классического gsu,
то есть формат приводится к актуальному виду.
"""

import sys
import os
import re
from functools import reduce
import itertools
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from gsu.gsu import Grid
from geom.vect import Vect
import numpy as np

# --------------------------------------------------------------------------------------------------

# Признак режима тишины (сообщения о ходе работы скрипта не выводятся).
is_silent = False

# --------------------------------------------------------------------------------------------------


def say(obj):
    """
    Вывод сообщения.

    Arguments:
        obj -- данные, которые модут быть выведены функцией print.
    """

    if not is_silent:
        print(obj)

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='lazurit_add_beta_from_ansys',
                                     description='Add Beta field from Ansys to Lazurit data.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('lazurit', help='Lazurit grid')
    parser.add_argument('ansys', help='Ansys grid')
    parser.add_argument('out', help='out file')
    args = parser.parse_args()

    say('Запуск lazurit_add_beta_from_ansys.')

    say('Входные данные :')
    say('  lazurit = %s,' % args.lazurit)
    say('  ansys   = %s,' % args.ansys)
    say('  out     = %s'  % args.out)

    # Проверка существования файла.
    if not os.path.exists(args.lazurit):
        raise Exception('файл %s не существует' % args.lazurit)
    elif not os.path.exists(args.ansys):
        raise Exception('файл %s не существует' % args.ansys)
    else:
        say('Файлы найдены.')

    lgrid = Grid()
    lgrid.load(args.lazurit)
    agrid = Grid()
    agrid.load(args.ansys)

    #
    adata = [(Vect.from_iterable(f.get_center()), f['Beta']) for f in agrid.Faces]
    for f in lgrid.Faces:
        c = Vect.from_iterable(f.get_center())
        a = np.array([(c - ac).mod2() for (ac, _) in adata])
        f['Beta'] = adata[a.argmin()][1]
    #

    lgrid.store(args.out)
    say('Файл %s сохранен.' % args.out)

# --------------------------------------------------------------------------------------------------
