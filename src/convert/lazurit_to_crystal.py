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

def extract_numbers_from_lines(lines):
    """
    Извлечение чисел из строк.
    :param lines: Строки.
    :return: Числа.
    """

    nlines = [[float(s) for s in line.split()] for line in lines]

    return list(itertools.chain(*nlines))

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
            i, j, k = int(ss[1]), int(ss[3]), int(ss[5])
            assert k == 1
            for _ in range(2):
                gl = gf.readline()
            for _ in range(6):
                sl = sf.readline()

            glines = gf.readlines()
            slines = sf.readlines()
            gnumbers = extract_numbers_from_lines(glines)
            snumbers = extract_numbers_from_lines(slines)

            print(len(gnumbers))
            print(len(snumbers))

            # Закрываем пару файлов.
            gf.close()
            sf.close()

        of.close()

# --------------------------------------------------------------------------------------------------
