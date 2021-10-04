#!/usr/bin/env python3
"""
Скрипт для конвертации входных данных в формате ANSYS во входные данные crystal.
Входные данные в формате ANSYS представляют собой пару файлов:
    <name>.fensap.dat
    <naem>.drop3d.dat
"""

import sys
import os
import re
import pathlib
import gsu_converter
from functools import reduce

# --------------------------------------------------------------------------------------------------

# Признак режима тишины (сообщения о ходе работы скрипта не выводятся).
is_silent = False

# Имена зон, которые надо извлечь.
zones_to_extract = []

# Названия переменных fensap (не обязательно полные), которые надо извлечь.
fensap_variables_to_extract = ['X', 'Y', 'Z', 'Static temperature', 'SF1', 'SF2', 'SF3', 'Classical heat flux']

# Названия переменных drop3d (не обязательно полные), которые надо извлечь.
drop3d_variables_to_extract = ['X', 'Y', 'Z', 'Collection efficiency-Droplet']

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


def is_str_title(s):
    """
    Проверка, что строка является титульной.

    Arguments:
        s -- строка.

    Result:
        True -- если строка является титульной,
        False -- в противном случае.
    """

    return 'TITLE' in s

# --------------------------------------------------------------------------------------------------


def is_str_variables(s):
    """
    Проверка, что строка является описанием переменных.

    Arguments:
        s -- строка.

    Result:
        True -- если строка является описанием переменных,
        False -- в противном случае.
    """

    return 'VARIABLES' in s

# --------------------------------------------------------------------------------------------------


def is_str_zone(s):
    """
    Проверка, что строка является описанием зоны.

    Arguments:
        s -- строка.

    Result:
        True -- если строка является описанием зоны,
        False -- в противном случае.
    """

    return 'ZONE' in s

# --------------------------------------------------------------------------------------------------


def is_str_meta(s):
    """
    Проверка, является ли строка строкой с метаданными.

    Arguments:
        s -- строка.

    Result:
        True -- если это строка с метаданными,
        False -- в противном случае.
    """

    return is_str_title(s) or is_str_variables(s) or is_str_zone(s)

# --------------------------------------------------------------------------------------------------


def is_str_zone_to_extract(s):
    """
    Проверка, описывает ли строка зону, которую необходимо извлечь.

    Arguments:
        s -- строка.

    Result:
        True -- если строка описывает зону, которую необходимо извлечь,
        False -- в противном случае.
    """

    if not is_str_zone(s):
        return False

    for zone_to_extract in zones_to_extract:
        str = 'T="%s"' % zone_to_extract
        if str in s:
            return True

    return False

# --------------------------------------------------------------------------------------------------


def get_zone_nodes_and_elements_count(s):
    """
    Получение из описания зоны количества узлов и элементов.

    Arguments:
        s -- описание зоны.

    Result:
        Кортеж, содержащий количество узлов и количество элементов.
    """

    # Отрезаем часть строки после 'N=' и ищем в ней все целые числа.
    nums = re.findall('(\d+)', s.split('N=')[1])

    return int(nums[0]), int(nums[1])

# --------------------------------------------------------------------------------------------------


def is_variable_names_match(fn, sn):
    """
    Проверка, что искомая переменная совпадает с переменной из файла.

    Arguments:
        fn -- переменная из файла,
        sn -- искомая переменная.

    Result:
        True -- если искомая переменная совпадает с переменной из файла,
        False -- в противном случае.
    """

    if len(sn) == 1:
        return sn == fn
    else:
        return sn in fn

# --------------------------------------------------------------------------------------------------


def variables_names_and_mask(fnames, snames):
    """
    Поиск маски переменных из файла (fnames), которые были найдены в списке snames.
    Переименование переменных в списке имен переменных из файла.

    Arguments:
        fnames -- имена переменных из файла,
        snames -- искомые имена.

    Result:
        Кортеж.
        Первый элемент кортежа - список переменных с переименовынными переменными.
        Второй элемент кортежа - маска переменных.
    """

    mask = [0] * len(fnames)

    for (fi, fn) in enumerate(fnames):
        for sn in snames:
            if is_variable_names_match(fn, sn):
                fnames[fi] = sn
                mask[fi] = 1

    return fnames, mask

# --------------------------------------------------------------------------------------------------


def filter_with_mask(a, m):
    """
    Фильтрация массива с использованием маски.

    Arguments:
        a -- массив,
        m -- маска.

    Result:
        Отфильтрованный массив.
    """

    z = list(zip(a, m))
    f = list(filter(lambda x: x[1] != 0, z))

    return [xa for (xa, xm) in f]

# --------------------------------------------------------------------------------------------------


def filter_zones(filename, ofilename):
    """
    Фильтрация зон.

    Arguments:
        filename -- имя входного файла,
        ofilename -- имя выходного файла.
    """

    say('  filter_zones : %s -> %s' % (filename, ofilename))

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Цикл по всем строкам входного файла.
            nodes_to_read, elements_to_read = 0, 0
            l = f.readline()
            while l:

                if is_str_meta(l):
                    if (nodes_to_read > 0) or (elements_to_read > 0):
                        raise Exception('внутренная ошибка при чтении сетки')
                    if is_str_title(l):
                        of.write(l)
                    elif is_str_variables(l):
                        of.write(l)
                    elif is_str_zone(l):
                        (nodes_count, elements_count) = get_zone_nodes_and_elements_count(l)
                        if is_str_zone_to_extract(l):
                            of.write(l)
                            nodes_to_read, elements_to_read = nodes_count, elements_count
                else:
                    if nodes_to_read > 0:
                        of.write(l)
                        nodes_to_read -= 1
                    elif elements_to_read > 0:
                        of.write(l)
                        elements_to_read -= 1

                l = f.readline()

    f.close()
    of.close()

# --------------------------------------------------------------------------------------------------


def filter_variables(filename, ofilename, variables_names):
    """
    Фильтрация переменных.

    Arguments:
        filename -- имя входного файла,
        ofilename -- имя выходного файла,
        variables_names -- список имен переменных.
    """

    say('  filter_variables : %s -> %s\n    переменные %s' % (filename, ofilename, variables_names))

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Цикл по всем строкам входного файла.
            nodes_to_read, elements_to_read = 0, 0
            l = f.readline()
            while l:

                if is_str_meta(l):
                    if (nodes_to_read > 0) or (elements_to_read > 0):
                        raise Exception('внутренная ошибка при чтении сетки')
                    if is_str_title(l):
                        of.write(l)
                    elif is_str_variables(l):
                        fnames = l.split('=')[1].replace('"', '').split(',')
                        (fnames, mask) = variables_names_and_mask(fnames, variables_names)
                        filtered_names = filter_with_mask(fnames, mask)
                        joined_names = ', '.join(['"' + name + '"' for name in filtered_names])
                        of.write('VARIABLES=' + joined_names + '\n')
                    elif is_str_zone(l):
                        (nodes_to_read, elements_to_read) = get_zone_nodes_and_elements_count(l)
                        of.write(l)
                else:
                    if nodes_to_read > 0:
                        floats = l.split()
                        filtered_floats = filter_with_mask(floats, mask)
                        joined_floats = ' '.join(filtered_floats)
                        of.write(joined_floats + '\n')
                        nodes_to_read -= 1
                    elif elements_to_read > 0:
                        of.write(l)
                        elements_to_read -= 1

                l = f.readline()

    f.close()
    of.close()

# --------------------------------------------------------------------------------------------------


def merge_fensap_drop3d(ffilename, dfilename, ofilename):
    """
    Слияние данный в один файл.

    Arguments:
        ffilename -- файл fensap,
        dfiliname -- файл drop3d,
        ofilename -- выходной файл.
    """

    say('  merge_fensap_drop3d : %s, %s -> %s' % (ffilename, dfilename, ofilename))

    with open(ffilename, 'r') as ff:
        with open(dfilename, 'r') as df:
            with open(ofilename, 'w') as of:

                # Цикл по всем строкам входного файла.
                nodes_to_read, elements_to_read = 0, 0
                fl, dl = ff.readline(), df.readline()
                while fl:

                    if is_str_meta(fl):
                        if (nodes_to_read > 0) or (elements_to_read > 0):
                            raise Exception('внутренная ошибка при чтении сетки')
                        if is_str_title(fl):
                            of.write(fl)
                        elif is_str_variables(fl):
                            of.write('VARIABLES = "X", "Y", "Z", "Flux", "Beta", ' \
                                     + '"TauX", "TauY", "TauZ", "Ts"\n')
                        elif is_str_zone(fl):
                            (nodes_to_read, elements_to_read) = get_zone_nodes_and_elements_count(fl)
                            of.write(fl)
                    else:
                        if nodes_to_read > 0:
                            ffields = fl.split()
                            dfields = dl.split()
                            fields = [ffields[0], ffields[1], ffields[2], ffields[7], dfields[3], \
                                      ffields[4], ffields[5], ffields[6], ffields[3]]
                            of.write(' '.join(fields) + '\n')
                            nodes_to_read -= 1
                        elif elements_to_read > 0:
                            of.write(fl)
                            elements_to_read -= 1

                    fl, dl = ff.readline(), df.readline()

    ff.close()
    df.close()
    of.close()

# --------------------------------------------------------------------------------------------------


def calculate_htc(filename, ofilename, ta):
    """
    Вычисление коэффициента теплоотдачи.

    Arguments:
        filename -- имя входного файла,
        ofilename -- имя выходного файла,
        ta -- температура свободного потока.
    """

    say('  calculate_htc : %s -> %s' % (filename, ofilename))

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Цикл по всем строкам входного файла.
            nodes_to_read, elements_to_read = 0, 0
            l = f.readline()
            while l:

                if is_str_meta(l):
                    if (nodes_to_read > 0) or (elements_to_read > 0):
                        raise Exception('внутренная ошибка при чтении сетки')
                    if is_str_title(l):
                        of.write(l)
                    elif is_str_variables(l):
                        of.write('VARIABLES = "X", "Y", "Z", "HTC", "Beta", ' \
                                 + '"TauX", "TauY", "TauZ"\n')
                    elif is_str_zone(l):
                        (nodes_to_read, elements_to_read) = get_zone_nodes_and_elements_count(l)
                        of.write(l)
                else:
                    if nodes_to_read > 0:

                        fields = l.split()
                        flux = float(fields[3])
                        static_temperature = float(fields[8])

                        # Вычисление HTC.
                        htc = flux / (ta - static_temperature)
                        fields[3] = str(htc)

                        # Отцепляем static temperature.
                        fields = fields[:-1]

                        of.write(' '.join(fields) + '\n')
                        nodes_to_read -= 1
                    elif elements_to_read > 0:
                        of.write(l)
                        elements_to_read -= 1

                l = f.readline()

    f.close()
    of.close()

# --------------------------------------------------------------------------------------------------


def data_from_nodes_to_faces(filename, ofilename, is_reversed_normals):
    """
    Перевод данных из узлов в ячейки.

    Arguments:
        filename -- имя входного файла,
        ofilename -- имя выходного файла,
        is_reversed_normals -- признак развернутых нормалей.
    """

    say('  data_from_nodes_to_faces : %s -> %s' % (filename, ofilename))
    g = gsu_converter.Grid()
    g.LoadNodesData(filename, is_reversed_normals)
    g.ConvertDataFromNodesToFaces()
    g.Export(ofilename)

# --------------------------------------------------------------------------------------------------


def print_help():
    """
    Печать help.
    """

    print('Использование скрипта:')
    print('  ./ansys_to_crystal basename=<базовое имя сетки>')
    print('                     zones=<список зон для фильтрации>')
    print('                     ta=<температура свободного потока>')
    print('                     normals=origin|reversed')

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    if len(sys.argv) < 5:
        print_help()
        exit(0)

    for arg in sys.argv[1:]:
        [par, val] = arg.split('=')

        if par == 'basename':
            pp = pathlib.PurePath(val)
            dirname, basename = pp.parent, pp.name
        elif par == 'zones':
            zones_to_extract = val.split(',')
        elif par == 'ta':
            free_stream_temperature = float(val)
        elif par == 'normals':
            if val == 'origin':
                reversed_normals = False
            elif val == 'reversed':
                reversed_normals = True
            else:
                raise Exception('неизвестный параметр определения нормалей')
        else:
            raise Exception('неизвестный параметр')

    say('Входные данные :')
    say('  dir = %s, base = %s' % (dirname, basename))
    say('  zones = %s' % zones_to_extract)
    say('  ta = %f' % free_stream_temperature)
    say('  normals = %s' % 'reversed' if reversed_normals else 'origin')

    # Функции сборки имени файла.
    def cc(d, ph, b, s):
        if ph == '':
            ph_str = ''
        else:
            ph_str = '{0}.'.format(ph)
        if s == '':
            s_str = ''
        else:
            s_str = '.{0}'.format(s)
        return '{0}/{1}{2}{3}.dat'.format(d, ph_str, b, s_str)
    def cc2(d, ph, b):
        return cc(d, ph, b, 'fensap'), cc(d, ph, b, 'drop3d')

    # Начальные имена.
    fensap_fn, drop3d_fn = cc2(dirname, '', basename)
    print(fensap_fn, drop3d_fn)

    # Проверка существования файлов.
    if not os.path.exists(fensap_fn):
        raise Exception('файл %s не существует' % fensap_fn)
    elif not os.path.exists(drop3d_fn):
        raise Exception('файл %s не существует' % drop3d_fn)
    else:
        say('Входные файлы найдены.')

    # Фильтрация зон.
    fensap_fn_ph1, drop3d_fn_ph1 = cc2(dirname, 'phase1', basename)
    say('phase1 : Фильтрация зон : %s' % str(zones_to_extract))
    filter_zones(fensap_fn, fensap_fn_ph1)
    filter_zones(drop3d_fn, drop3d_fn_ph1)

    # Фильтрация переменных.
    fensap_fn_ph2, drop3d_fn_ph2 = cc2(dirname, 'phase2', basename)
    say('phase2 : Фильтрация переменных\n  fensap %s\n  drop3d %s'
        % (str(fensap_variables_to_extract), str(drop3d_variables_to_extract)))
    filter_variables(fensap_fn_ph1, fensap_fn_ph2, fensap_variables_to_extract)
    filter_variables(drop3d_fn_ph1, drop3d_fn_ph2, drop3d_variables_to_extract)

    # Слияние файлов данных fensap и drop3d.
    fn_ph3 = cc(dirname, 'phase3', basename, '')
    say('phase3 : Слияние в единый файл')
    merge_fensap_drop3d(fensap_fn_ph2, drop3d_fn_ph2, fn_ph3)

    # Пересчет HTC.
    fn_ph4 = cc(dirname, 'phase4', basename, '')
    say('phase4 : Вычисление HTC')
    calculate_htc(fn_ph3, fn_ph4, free_stream_temperature)

    # Перевод данных с узлов в ячейки.
    fn_ph5 = cc(dirname, '', basename, '')
    say('phase5 : Перевод данных из узлов на грани')
    data_from_nodes_to_faces(fn_ph4, fn_ph5, reversed_normals)

# --------------------------------------------------------------------------------------------------
