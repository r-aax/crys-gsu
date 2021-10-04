#!/usr/bin/env python3
"""
Скрипт для конвертации входных данных в формате lazurit во входные данные crystal.
Входные данные в формате lazurit представляет собой единый файл:
"""

import sys
import os
import re
import gsu_converter
from functools import reduce
import itertools

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

def leap_lines(lines):
    """
    Соединение строк.
    :param lines: Строки.
    :return: Строки элементов.
    """

    nlines = [[s for s in line.split()] for line in lines]

    return list(itertools.chain(*nlines))

# --------------------------------------------------------------------------------------------------

def get_block_name(s):
    """
    Get name of block.
    :param s: Name of block.
    :return: Name of block.
    """

    for i in range(len(s) - 1, -1, -1):
        if s[i] == 'B':
            return s[i:-4]

# --------------------------------------------------------------------------------------------------

def double_list(li):
    """
    Double elements in the list.
    :param li: List.
    :return: Doubled list.
    """

    dli = [[e, e] for e in li]

    return list(itertools.chain(*dli))


# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    if (len(sys.argv)) < 3:
        print('Использование скрипта:')
        print('  ./lazurit_to_crystal in_dir=<имя входной директории>')
        print('                       out_file=<имя выходного файла>')
        exit(0)

    say('Запуск lazurit_to_crystal.')

    for arg in sys.argv[1:]:
        [par, val] = arg.split('=')

        if par == 'in_dir':
            in_dir = val
        elif par == 'out_file':
            out_file = val
        else:
            raise Exception('неизвестный параметр')

    say('Входные данные :')
    say('  in_dir = %s,' % in_dir)
    say('  out_file = %s' % out_file)

    # Проверка существования файла.
    if not os.path.exists(in_dir):
        raise Exception('директория %s не существует' % in_dir)
    else:
        say('Входная директория найдена.')

    # В списке всех файлов сначала идут сетки, потом решения.
    all_files = [fn for fn in os.listdir(in_dir) if '.TEC' in fn]
    n = len(all_files)

    with open(out_file, 'w') as of:
        of.write('# EXPORT_MODE=CHECK_POINT\n')
        of.write('TITLE="FE Surface Data ASCII"\n')
        of.write('VARIABLES="X", "Y", "Z", "T", "Hw", "Hi", "HTC", "Beta", "TauX", "TauY", "TauZ"\n')

        # Проходим по парам всех файлов.
        for fi in range(n // 2):
            gfn, sfn = all_files[fi], all_files[fi + n // 2]
            say('Обработка файлов {0}, {1}'.format(gfn, sfn))

            # Открываем пару файлов.
            gf = open('{0}/{1}'.format(in_dir, gfn), 'r')
            sf = open('{0}/{1}'.format(in_dir, sfn), 'r')

            # Читаем все, что нужно.
            for _ in range(4):
                gl = gf.readline()
            ss = gl.split()
            zone_i, zone_j, zone_k = int(ss[1]), int(ss[3]), int(ss[5])
            nc = zone_i * zone_j * zone_k
            hec = (zone_i - 1) * (zone_j - 1)
            assert zone_k == 1
            for _ in range(2):
                gl = gf.readline()
            for _ in range(6):
                sl = sf.readline()

            glines = gf.readlines()
            slines = sf.readlines()
            gelems = leap_lines(glines)
            selems = leap_lines(slines)

            of.write('ZONE T="{0}"\n'.format(get_block_name(gfn)))
            of.write('NODES={0}\n'.format(nc))
            of.write('ELEMENTS={0}\n'.format(2 * hec))
            of.write('DATAPACKING=BLOCK\n')
            of.write('ZONETYPE=FETRIANGLE\n')
            of.write('VARLOCATION=([4-11]=CELLCENTERED)\n')

            # Печать координат.
            of.write(' '.join(gelems[0: nc]) + '\n')
            of.write(' '.join(gelems[nc: 2 * nc]) + '\n')
            of.write(' '.join(gelems[2 * nc: 3 * nc]) + '\n')

            # Печать данных.
            of.write(' '.join(double_list(selems[0: hec])) + '\n')
            of.write(' '.join(double_list(selems[hec: 2 * hec])) + '\n')
            of.write(' '.join(double_list(selems[2 * hec: 3 * hec])) + '\n')
            of.write(' '.join(double_list(selems[3 * hec: 4 * hec])) + '\n')
            of.write(' '.join(double_list(selems[4 * hec: 5 * hec])) + '\n')
            of.write(' '.join(double_list(selems[5 * hec: 6 * hec])) + '\n')
            of.write(' '.join(double_list(selems[6 * hec: 7 * hec])) + '\n')
            of.write(' '.join(double_list(selems[7 * hec: 8 * hec])) + '\n')

            # Печать линков.
            for i in range(zone_i - 1):
                for j in range(zone_j - 1):
                    ni00 = (i)     + (j)     * zone_i
                    ni01 = (i)     + (j + 1) * zone_i
                    ni10 = (i + 1) + (j)     * zone_i
                    ni11 = (i + 1) + (j + 1) * zone_i
                    of.write('{0} {1} {2}\n'.format(ni00 + 1, ni01 + 1, ni11 + 1))
                    of.write('{0} {1} {2}\n'.format(ni00 + 1, ni11 + 1, ni10 + 1))

            # Закрываем пару файлов.
            gf.close()
            sf.close()

        of.close()

# --------------------------------------------------------------------------------------------------
