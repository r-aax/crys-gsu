#!/usr/bin/env python3
"""
Скрипт для конвертации входных данных в формате lazurit во входные данные crystal.
Входные данные в формате lazurit представляет собой единый файл:
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
import re
from functools import reduce
import itertools
import gsu

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


def double_and_neg_list(li):
    """
    Double elements in the list and negate all of them.
    :param li: List.
    :return: Doubled and negated list.
    """

    dli = [[str(-float(e)), str(-float(e))] for e in li]

    return list(itertools.chain(*dli))

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='lazurit_to_crystal',
                                     description='Convert data from Lazurit to Crystal.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_dir', help='input dir')
    parser.add_argument('out_file', help='output file')
    args = parser.parse_args()

    say('Запуск lazurit_to_crystal.')

    say('Входные данные :')
    say('  in_dir = %s,' % args.in_dir)
    say('  out_file = %s' % args.out_file)

    # Проверка существования файла.
    if not os.path.exists(args.in_dir):
        raise Exception('директория %s не существует' % args.in_dir)
    else:
        say('Входная директория найдена.')

    # В списке всех файлов сначала должны идти все сетки, потом все решения.
    # Так как функция listdir этого не гарантирует то сортируем список вручную.
    # Ниже есть проверки корректности списка.
    all_files = sorted([fn for fn in os.listdir(args.in_dir) if '.TEC' in fn])
    n = len(all_files)

    # Если все верно, то количество файлов должно быть четным.
    if n % 2 == 1:
        raise Exception('количество файлов *.TEC во входной директории должно быть четным')

    with open(args.out_file, 'w') as of:
        of.write('# EXPORT_MODE=CHECK_POINT\n')
        of.write('TITLE="FE Surface Data ASCII"\n')
        of.write('VARIABLES="X", "Y", "Z", "T", "Hw", "Hi", "HTC", "Beta", "TauX", "TauY", "TauZ"\n')

        # Проходим по парам всех файлов.
        for fi in range(n // 2):
            gfn, sfn = all_files[fi], all_files[fi + n // 2]
            say('Обработка файлов {0}, {1}'.format(gfn, sfn))

            # Первый файл должен быть файлом сетки.
            if gfn[:14] != 'TEC_WALL_Grid_':
                raise Exception('файл {0} не является файлом сетки'.format(gfn))

            # Второй файл должен быть файлом решения.
            if sfn[:18] != 'TEC_WALL_Solution_':
                raise Exception('файл {0} не является файлом решения'.format(sfn))

            # Файл сетки и файл решения должны принадлежать одному и тому же блоку.
            if gfn[14:] != sfn[18:]:
                raise Exception('файлы {0} и {1} не соответствуют одному блоку'.format(gfn, sfn))

            # Открываем пару файлов.
            gf = open('{0}/{1}'.format(args.in_dir, gfn), 'r')
            sf = open('{0}/{1}'.format(args.in_dir, sfn), 'r')

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
            # Данные нужно удваивать, так как на вход подаются 4-угольники,
            # а нам нужны треугольники (каждый четырехугольник делится пополам).
            of.write(' '.join(double_list(selems[0: hec])) + '\n')
            of.write(' '.join(double_list(selems[hec: 2 * hec])) + '\n')
            of.write(' '.join(double_list(selems[2 * hec: 3 * hec])) + '\n')
            of.write(' '.join(double_list(selems[3 * hec: 4 * hec])) + '\n')
            of.write(' '.join(double_list(selems[4 * hec: 5 * hec])) + '\n')
            of.write(' '.join(double_and_neg_list(selems[5 * hec: 6 * hec])) + '\n')
            of.write(' '.join(double_and_neg_list(selems[6 * hec: 7 * hec])) + '\n')
            of.write(' '.join(double_and_neg_list(selems[7 * hec: 8 * hec])) + '\n')

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

    # Читаем сетку с помощью gsu и сохраняем обратно -> появятся нужные поля.
    g2 = gsu.Grid()
    g2.load(args.out_file)
    g2.store(args.out_file)

# --------------------------------------------------------------------------------------------------
