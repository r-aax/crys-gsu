#!/usr/bin/env python3
"""
Скрипт для конвертации входных данных в формате ANSYS во входные данные crystal.
Входные данные в формате ANSYS представляют собой пару файлов:
    <name>.fensap.dat
    <naem>.drop3d.dat
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
import re
import pathlib
import gsu_converter
from functools import reduce
from gsu.gsu import Grid

# --------------------------------------------------------------------------------------------------

# Признак режима тишины (сообщения о ходе работы скрипта не выводятся).
is_silent = False

# Имена зон, которые надо извлечь.
zones_to_extract = []

# Названия переменных fensap (не обязательно полные), которые надо извлечь.
fensap_variables_to_extract = ['X', 'Y', 'Z', 'Static temperature', 'SF1', 'SF2', 'SF3', 'Classical heat flux']

# Названия переменных drop3d (не обязательно полные), которые надо извлечь.
drop3d_variables_to_extract = ['X', 'Y', 'Z', 'Collection efficiency-Droplet']

# Названия переменных adiabatic (не обязательно полные), которые нужно извлечь.
adiabatic_variables_to_extract = ['X', 'Y', 'Z', 'Static temperature']

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

    say('    filter_zones : %s -> %s' % (filename, ofilename))

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

    say('    filter_variables : %s -> %s\n    переменные %s' % (filename, ofilename, variables_names))

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


def merge_fensap_drop3d_adiabatic(ffilename, dfilename, afilename, ofilename):
    """
    Слияние данный в один файл.

    Arguments:
        ffilename -- файл fensap,
        dfiliname -- файл drop3d,
        afilename -- файл adiabatic,
        ofilename -- выходной файл.
    """

    say('    merge_fensap_drop3d_adiabatic : \n    %s, \n    %s, \n    %s \n    -> %s' % (ffilename, dfilename, afilename, ofilename))

    with open(ffilename, 'r') as ff:
        with open(dfilename, 'r') as df:
            with open(afilename, 'r') as af:
                with open(ofilename, 'w') as of:

                    # Цикл по всем строкам входного файла.
                    nodes_to_read, elements_to_read = 0, 0
                    fl, dl, al = ff.readline(), df.readline(), af.readline()
                    while fl:

                        if is_str_meta(fl):
                            if (nodes_to_read > 0) or (elements_to_read > 0):
                                raise Exception('внутренная ошибка при чтении сетки')
                            if is_str_title(fl):
                                of.write(fl)
                            elif is_str_variables(fl):
                                of.write('VARIABLES = "X", "Y", "Z", "Flux", "Beta", ' \
                                         + '"TauX", "TauY", "TauZ", "Ts", "StatT"\n')
                            elif is_str_zone(fl):
                                (nodes_to_read, elements_to_read) = get_zone_nodes_and_elements_count(fl)
                                of.write(fl)
                        else:
                            if nodes_to_read > 0:
                                ffields = fl.split()
                                dfields = dl.split()
                                afields = al.split()
                                fields = [ffields[0], ffields[1], ffields[2], ffields[7], dfields[3], \
                                          ffields[4], ffields[5], ffields[6], ffields[3], afields[3]]
                                of.write(' '.join(fields) + '\n')
                                nodes_to_read -= 1
                            elif elements_to_read > 0:
                                of.write(fl)
                                elements_to_read -= 1

                        fl, dl, al = ff.readline(), df.readline(), af.readline()

        ff.close()
        df.close()
        af.close()
        of.close()

# --------------------------------------------------------------------------------------------------


def calculate_htc_recovery_factor(filename, ofilename, air_temperature, air_velocity, heat_capacity):
    """
    Вычисление коэффициента теплоотдачи.

    Arguments:
        filename -- имя входного файла,
        ofilename -- имя выходного файла,
        air_temperature -- температура свободного потока,
        air_velocity -- скорость свободного потока,
        heat_capacity -- темлоемкость воздуха.
    """

    say('    calculate_htc : %s -> %s' % (filename, ofilename))

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
                                 + '"TauX", "TauY", "TauZ", "RecoveryFactor"\n')
                    elif is_str_zone(l):
                        (nodes_to_read, elements_to_read) = get_zone_nodes_and_elements_count(l)
                        of.write(l)
                else:
                    if nodes_to_read > 0:

                        fields = l.split()
                        flux = -float(fields[3])
                        t_surf = float(fields[8])
                        static_temperature = float(fields[9])

                        # Вычисление RecoveryFactor.
                        k = (2.0 * heat_capacity) / (air_velocity * air_velocity)
                        recovery_factor = (static_temperature - air_temperature) * k

                        # Вычисление HTC.
                        htc = flux / (t_surf - (air_temperature + recovery_factor / k))
                        fields[3] = str(htc)

                        # Отцепляем последнее поле.
                        fields = fields[:-1]

                        # И записываем туда RecoveryFactor.
                        fields[-1] = str(recovery_factor)

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

    g = gsu_converter.Grid()
    g.LoadNodesData(filename, is_reversed_normals)
    g.ConvertDataFromNodesToFaces()
    g.Export(ofilename)

    # Читаем сетку с помощью gsu и сохраняем обратно -> появятся нужные поля.
    g2 = Grid()
    g2.load(ofilename)
    g2.store(ofilename)

    say('    data_from_nodes_to_faces : %s -> %s' % (filename, ofilename))

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='ansys_to_crystal',
                                     description='Convert data from Ansys to Crystal.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('basename',
                        help='basename for *.fensap.dat, *.drop3d.dat, *.adiabatic.dat files')
    parser.add_argument('zones', help='zones to extract')
    parser.add_argument('-t', '--air_termerature', dest='air_temperature',
                        type=float, help='temperature of free stream (K)')
    parser.add_argument('-v', '--air_velocity', dest='air_velocity',
                        type=float, help='velocity of free stream (m / s)')
    parser.add_argument('-cp', '--heat_capacity', dest='heat_capacity',
                        type=float, help='heat capacity of free stream (J / (kg * K))')
    parser.add_argument('-n', '--normals', dest='normals', choices=['origin', 'reversed'], default='reversed',
                        help='normals orientation')
    args = parser.parse_args()

    pp = pathlib.PurePath(args.basename)
    dirn, basen = pp.parent, pp.name
    zones_to_extract = args.zones.split(',')
    reversed_normals = (args.normals == 'reversed')

    say('Входные данные :')
    say('    dir / base      = %s / %s' % (dirn, basen))
    say('    zones           = %s' % zones_to_extract)
    say('    air_temperature = %f' % args.air_temperature)
    say('    air_velocity    = %f' % args.air_velocity)
    say('    heat_capacity   = %f' % args.heat_capacity)
    say('    normals         = %s' % args.normals)

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
    def cc3(d, ph, b):
        return cc(d, ph, b, 'fensap'), cc(d, ph, b, 'drop3d'), cc(d, ph, b, 'adiabatic')

    # Начальные имена.
    fensap_fn, drop3d_fn, adiabatic_fn = cc3(dirn, '', basen)
    print(fensap_fn, drop3d_fn, adiabatic_fn)

    # Проверка существования файлов.
    if not os.path.exists(fensap_fn):
        raise Exception('файл %s не существует' % fensap_fn)
    elif not os.path.exists(drop3d_fn):
        raise Exception('файл %s не существует' % drop3d_fn)
    elif not os.path.exists(adiabatic_fn):
        raise Exception('файл %s не существует' % adiabatic_fn)
    else:
        say('Входные файлы найдены.')

    # Фильтрация зон.
    fensap_fn_ph1, drop3d_fn_ph1, adiabatic_fn_ph1 = cc3(dirn, 'phase1', basen)
    say('phase1 : Фильтрация зон : %s' % str(zones_to_extract))
    filter_zones(fensap_fn, fensap_fn_ph1)
    filter_zones(drop3d_fn, drop3d_fn_ph1)
    filter_zones(adiabatic_fn, adiabatic_fn_ph1)

    # Фильтрация переменных.
    fensap_fn_ph2, drop3d_fn_ph2, adiabatic_fn_ph2 = cc3(dirn, 'phase2', basen)
    say('phase2 : Фильтрация переменных\n    fensap %s\n    drop3d %s\n    adiabatic %s'
        % (str(fensap_variables_to_extract), str(drop3d_variables_to_extract), str(adiabatic_variables_to_extract)))
    filter_variables(fensap_fn_ph1, fensap_fn_ph2, fensap_variables_to_extract)
    filter_variables(drop3d_fn_ph1, drop3d_fn_ph2, drop3d_variables_to_extract)
    filter_variables(adiabatic_fn_ph1, adiabatic_fn_ph2, adiabatic_variables_to_extract)

    # Слияние файлов данных fensap, drop3d и adiabatic.
    fn_ph3 = cc(dirn, 'phase3', basen, '')
    say('phase3 : Слияние в единый файл')
    merge_fensap_drop3d_adiabatic(fensap_fn_ph2, drop3d_fn_ph2, adiabatic_fn_ph2, fn_ph3)

    # Пересчет HTC и RecoveryFactor.
    fn_ph4 = cc(dirn, 'phase4', basen, '')
    say('phase4 : Вычисление HTC')
    calculate_htc_recovery_factor(fn_ph3, fn_ph4, args.air_temperature, args.air_velocity, args.heat_capacity)

    # Перевод данных с узлов в ячейки.
    fn_ph5 = cc(dirn, '', basen, '')
    say('phase5 : Перевод данных из узлов на грани')
    data_from_nodes_to_faces(fn_ph4, fn_ph5, reversed_normals)

# --------------------------------------------------------------------------------------------------
