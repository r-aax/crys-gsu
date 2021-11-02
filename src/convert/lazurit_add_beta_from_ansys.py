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
import gsu
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

    if (len(sys.argv)) < 3:
        print('Использование скрипта:')
        print('  ./lazurit_add_beta_from_ansys lazurit=<сетка lazurit>')
        print('                                ansys=<сетка ansys>')
        print('                                out=<выходная сетка>')
        exit(0)

    say('Запуск lazurit_add_beta_from_ansys.')

    for arg in sys.argv[1:]:
        [par, val] = arg.split('=')

        if par == 'lazurit':
            lazurit_grid = val
        elif par == 'ansys':
            ansys_grid = val
        elif par == 'out':
            out_grid = val
        else:
            raise Exception('неизвестный параметр')

    say('Входные данные :')
    say('  lazurit = %s,' % lazurit_grid)
    say('  ansys   = %s,' % ansys_grid)
    say('  out     = %s'  % out_grid)

    # Проверка существования файла.
    if not os.path.exists(lazurit_grid):
        raise Exception('файл %s не существует' % lazurit_grid)
    elif not os.path.exists(ansys_grid):
        raise Exception('файл %s не существует' % ansys_grid)
    else:
        say('Файлы найдены.')

    lgrid = gsu.Grid()
    lgrid.load(lazurit_grid)
    agrid = gsu.Grid()
    agrid.load(ansys_grid)

    #
    adata = [(Vect.from_iterable(f.get_center()), f.get_beta()) for f in agrid.Faces]
    for f in lgrid.Faces:
        c = Vect.from_iterable(f.get_center())
        a = np.array([(c - ac).mod2() for (ac, _) in adata])
        f.set_beta(adata[a.argmin()][1])
    #

    lgrid.store(out_grid)
    say('Файл %s сохранен.' % out_grid)

# --------------------------------------------------------------------------------------------------
