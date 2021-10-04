#!/usr/bin/env python3
"""
Скрипт для конвертации входных данных в формате lazurit во входные данные crystal.
Входные данные в формате lazurit представляет собой единый файл:
"""

import sys
import os
import re
import gsu
from functools import reduce

#---------------------------------------------------------------------------------------------------
# Глобальные переменные.
#---------------------------------------------------------------------------------------------------

# Признак режима тишины (сообщения о ходе работы скрипта не выводятся).
is_silent = False

#---------------------------------------------------------------------------------------------------
# Общие функции.
#---------------------------------------------------------------------------------------------------

def say(obj):
    """
    Вывод сообщения.

    Arguments:
        obj -- данные, которые модут быть выведены функцией print.
    """

    if not is_silent:
        print(obj)

#---------------------------------------------------------------------------------------------------
# Запуск скрипта.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    say('Запуск lazurit_to_crystal.')

    if (len(sys.argv)) < 3:
        print('Использование скрипта:')
        print('  ./lazurit_to_crystal in=<имя входного файла>')
        print('                       out=<имя выходного файла>')
        raise Exception('переданы не все входные данные')
    else:
        for arg in sys.argv[1:]:
            splitarg = arg.split('=')
            par = splitarg[0]
            val = splitarg[1]

            if par == 'in':
                in_filename = val
            elif par == 'out':
                out_filename = val
            else:
                raise Exception('неизвестный параметр')

        say('Входные данные :')
        say('  in = %s,' % in_filename)
        say('  out = %s.' % out_filename)

    # Проверка существования файла.
    if not os.path.exists(in_filename):
        raise Exception('файл %s не существует' % in_filename)
    else:
        say('Входной файл найден.')

    # Чтение файла.
    g = gsu.Grid()
    with open(in_filename, 'r') as f:

        l = f.readline()
        while l:

            if 'VARIABLES' in l:
                pass
            elif 'ZONE' in l:
                ss = l.split()
                zone_title = ss[1][3:-1]
                zone_i = int(ss[3])
                zone_j = int(ss[5])
                zone_k = int(ss[7])
                if zone_k != 1:
                    raise Exception('размер зоны K должен быть равен единице')
                nodes_count = zone_i * zone_j * zone_k
                # Добавление новой зоны.
                zone = gsu.Zone(zone_title)
                g.Zones.append(zone)
                # Зачитываем все узлы.
                for i in range(nodes_count):
                    l = f.readline()
                    ss = l.split()
                    zone.Nodes.append(gsu.Node(len(zone.Nodes),
                                               float(ss[0]), float(ss[1]), float(ss[2]),
                                               float(ss[13]), float(ss[14]),
                                               float(ss[10]), float(ss[11]), float(ss[12])))
                # Теперь добавляем все зоны (четырехугольные).
                for i in range(zone_i - 1):
                    for j in range(zone_j - 1):
                        ni00 = (i)     + (j)     * zone_i
                        ni01 = (i)     + (j + 1) * zone_i
                        ni10 = (i + 1) + (j)     * zone_i
                        ni11 = (i + 1) + (j + 1) * zone_i
                        n00 = zone.Nodes[ni00]
                        n01 = zone.Nodes[ni01]
                        n10 = zone.Nodes[ni10]
                        n11 = zone.Nodes[ni11]
                        zone.Faces.append(gsu.Face(n00, n01, n11))
                        zone.Faces.append(gsu.Face(n00, n11, n10))
            else:
                raise Exception('unexpected line : %s' % l)

            l = f.readline()

        # Сохраняем.
        g.ConvertDataFromNodesToFaces()
        g.Export(out_filename)

#---------------------------------------------------------------------------------------------------
