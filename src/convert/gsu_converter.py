#!/usr/bin/env python3
"""
Функционал для работы с сеткой.
"""

import re
from functools import reduce

# --------------------------------------------------------------------------------------------------

class Vector:

    # ----------------------------------------------------------------------------------------------

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        """
        Конструктор.

        Arguments:
            x -- первая координата,
            y -- вторая координата,
            z -- третья координата.
        """

        self.X = x
        self.Y = y
        self.Z = z

# --------------------------------------------------------------------------------------------------

class Node:

    # ----------------------------------------------------------------------------------------------

    def __init__(self, id = 0,
                 x = 0.0, y = 0.0, z = 0.0, htc = 0.0, beta = 0.0,
                 taux = 0.0, tauy = 0.0, tauz = 0.0):
        """
        Конструктор.

        Arguments:
            id -- идентификатор узла,
            x -- первая координата,
            y -- вторая координата,
            z -- третья координата,
            htc -- коэффициент теплоотдачи,
            beta -- коэффициент улавливания влаги,
            taux -- первая координата вектора напряжения сдвига,
            tauy -- вторая координата вектора напряжения сдвига,
            tauz -- третья координата вектора напряжения сдвига.
        """

        self.Id = id
        self.Point = Vector(x, y, z)
        self.T = 0.0
        self.Hw = 0.0
        self.Hi = 0.0
        self.HTC = htc
        self.Beta = beta
        self.Tau = Vector(taux, tauy, tauz)

    # ----------------------------------------------------------------------------------------------

    def Copy(self):
        """
        Копия.

        Result:
            Копия.
        """

        return Node(id = self.Id,
                    x = self.Point.X,
                    y = self.Point.Y,
                    z = self.Point.Z,
                    htc = self.HTC,
                    beta = self.Beta,
                    taux = self.Tau.X,
                    tauy = self.Tau.Y,
                    tauz = self.Tau.Z)

# --------------------------------------------------------------------------------------------------

class Face:

    # ----------------------------------------------------------------------------------------------

    def __init__(self, n1, n2, n3,
                 htc = 0.0, beta = 0.0, taux = 0.0, tauy = 0.0, tauz = 0.0):
        """
        Конструктор.

        Arguments:
            n1 -- первый узел,
            n2 -- второй узел,
            n3 -- третий узел,
            htc -- коэффициент теплоотдачи,
            beta -- коэффициент улавливания влаги,
            taux -- первая координата вектора напряжения сдвига,
            tauy -- вторая координата вектора напряжения сдвига,
            tauz -- третья координата вектора напряжения сдвига.
        """

        self.Nodes = [n1, n2, n3]
        self.T = 0.0
        self.Hw = 0.0
        self.Hi = 0.0
        self.HTC = htc
        self.Beta = beta
        self.Tau = Vector(taux, tauy, tauz)

# --------------------------------------------------------------------------------------------------

class Zone:

    # ----------------------------------------------------------------------------------------------

    def __init__(self, title):
        """
        Конструктор.

        Arguments:
            title -- имя зоны.
        """

        self.Title = title
        self.Nodes = []
        self.Faces = []

    # ----------------------------------------------------------------------------------------------

    def ReplaceNode(self, on, nn):
        """
        Замена узла.

        Arguments:
            on -- старый узел,
            nn -- новый узел.
        """

        for (i, n) in enumerate(self.Nodes):
            if n == on:
                self.Nodes[i] = nn
        for f in self.Faces:
            for (i, n) in enumerate(f.Nodes):
                if n == on:
                    f.Nodes[i] = nn

    # ----------------------------------------------------------------------------------------------

    def RestoreNodesIndices(self):
        """
        Восстановление индексов узлов.
        """

        for (i, n) in enumerate(self.Nodes):
            n.Id = i

    # ----------------------------------------------------------------------------------------------

    def DeleteExtraNodes(self):
        """
        Удаление лишних узлов.
        """

        ns = []
        for f in self.Faces:
            ns = ns + f.Nodes
        self.Nodes = [n for n in self.Nodes if n in ns]

# --------------------------------------------------------------------------------------------------

class Grid:

    # ----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Конструктор.
        """

        self.Zones = []

    # ----------------------------------------------------------------------------------------------

    @staticmethod
    def parse_coordinate(s : str) -> float:
        try:
            x = float(s)
        except ValueError:
            x = 0.0
            print('Не удалось распознать {}\nПринудительно установлено значение 0.0'.format(s))
        return x

    # ----------------------------------------------------------------------------------------------

    def LoadNodesData(self, filename, is_reversed_normals):
        """
        Загрузка данных в узлах.

        Arguments:
            filename -- имя файла,
            is_reversed_normals -- признак развернутых нормалей.
        """

        with open(filename, 'r') as f:

            # Цикл по всем строкам входного файла.
            nodes_to_read, elements_to_read = 0, 0
            l = f.readline()
            while l:

                if 'TITLE' in l:
                    pass
                elif 'VARIABLES' in l:
                    pass
                elif 'ZONE' in l:
                    nums = re.findall('(\d+)', l.split('N=')[1])
                    nodes_to_read, elements_to_read = int(nums[0]), int(nums[1])
                    self.Zones.append(Zone(l.split('T=')[1].split('"')[1]))
                else:
                    if nodes_to_read > 0:
                        [x, y, z, htc, beta, taux, tauy, tauz] = l.split()

                        # Считываем координаты и если находим значение, которое
                        # не удается распарсить, то устанавливаем его в 0.0 и выводим предупрежд.
                        x = self.parse_coordinate(x)
                        y = self.parse_coordinate(y)
                        z = self.parse_coordinate(z)

                        zone = self.Zones[-1]
                        zone.Nodes.append(Node(len(zone.Nodes),
                                               float(x), float(y), float(z),
                                               float(htc), float(beta),
                                               float(taux), float(tauy), float(tauz)))
                        nodes_to_read -= 1
                    elif elements_to_read > 0:
                        [s1, s2, s3, s4] = l.split()
                        zone = self.Zones[-1]

                        # Считаем, что зоны только четырехугольные.
                        n1 = zone.Nodes[int(s1) - 1]
                        n2 = zone.Nodes[int(s2) - 1]
                        n3 = zone.Nodes[int(s3) - 1]
                        n4 = zone.Nodes[int(s4) - 1]

                        faces = zone.Faces

                        # #bug 388
                        # Нормали могу использоваться в исходном виде,
                        # а могут в развернутом.
                        if is_reversed_normals:
                            n1, n2, n3, n4 = n4, n3, n2, n1

                        # #bug 384
                        # Четыре узла могут описывать треугольник,
                        # в этом случае какие-то две вершины должны совпадать,
                        # причем совпадать могут любые две вершины.
                        # Проверяем все варианты.
                        if n1 == n2:
                            faces.append(Face(n1, n3, n4))
                        elif n1 == n3:
                            faces.append(Face(n1, n2, n4))
                        elif n1 == n4:
                            faces.append(Face(n1, n2, n3))
                        elif n2 == n3:
                            faces.append(Face(n1, n2, n4))
                        elif n2 == n4:
                            faces.append(Face(n1, n2, n3))
                        elif n3 == n4:
                            faces.append(Face(n1, n2, n3))
                        else:
                            faces.append(Face(n1, n2, n3))
                            faces.append(Face(n1, n3, n4))

                        elements_to_read -= 1

                l = f.readline()

        f.close()

    # ----------------------------------------------------------------------------------------------

    def PrintInfo(self):
        """
        Печать информации.
        """

        for zone in self.Zones:
            print("ZONE %s (%d nodes, %d faces)" % (zone.Title, len(zone.Nodes), len(zone.Faces)))

    # ----------------------------------------------------------------------------------------------

    def ConvertDataFromNodesToFaces(self):
        """
        Конвертация данных из узлов в грани.
        """

        for z in self.Zones:
            for f in z.Faces:
                [n1, n2, n3] = f.Nodes
                f.T = (n1.T + n2.T + n3.T) / 3.0
                f.Hw = (n1.Hw + n2.Hw + n3.Hw) / 3.0
                f.Hi = (n1.Hi + n2.Hi + n3.Hi) / 3.0
                f.HTC = (n1.HTC + n2.HTC + n3.HTC) / 3.0
                f.Beta = (n1.Beta + n2.Beta + n3.Beta) / 3.0
                f.Tau = Vector((n1.Tau.X + n2.Tau.X + n3.Tau.X) / 3.0,
                               (n1.Tau.Y + n2.Tau.Y + n3.Tau.Y) / 3.0,
                               (n1.Tau.Z + n2.Tau.Z + n3.Tau.Z) / 3.0)

    # ----------------------------------------------------------------------------------------------

    def Import(self, filename):
        """
        Импорт сетки.

        Arguments:
            filename -- имя файла.
        """

        with open(filename, 'r') as f:
            l = f.readline()
            while l:

                # Обработка строки.
                if 'EXPORT_MODE' in l:
                    # Пропускаем шапку.
                    pass
                elif 'TITLE' in l:
                    # Пропускаем заголовок.
                    pass
                elif 'VARIABLES' in l:
                    # Пропускаем список переменных.
                    pass
                elif 'ZONE T=' in l:
                    # Зачитываем зону.
                    z = Zone(l.split('=')[1][1:-2])
                    self.Zones.append(z)
                    l = f.readline()
                    nodes_count = int(l.split('=')[1])
                    l = f.readline()
                    elements_count = int(l.split('=')[1])
                    l = f.readline()
                    l = f.readline()
                    l = f.readline()
                    # Пошли координаты.
                    l = f.readline()
                    xs = l.split()
                    l = f.readline()
                    ys = l.split()
                    l = f.readline()
                    zs = l.split()
                    assert len(xs) == nodes_count
                    assert len(ys) == nodes_count
                    assert len(zs) == nodes_count
                    for i in range(nodes_count):
                        z.Nodes.append(Node(id = i,
                                            x = float(xs[i]),
                                            y = float(ys[i]),
                                            z = float(zs[i])))
                    # Пошли данные, но три строки пропускаем.
                    l = f.readline()
                    l = f.readline()
                    l = f.readline()
                    l = f.readline()
                    htcs = l.split()
                    l = f.readline()
                    betas = l.split()
                    l = f.readline()
                    tauxs = l.split()
                    l = f.readline()
                    tauys = l.split()
                    l = f.readline()
                    tauzs = l.split()
                    assert len(betas) == elements_count
                    assert len(htcs) == elements_count
                    assert len(tauxs) == elements_count
                    assert len(tauys) == elements_count
                    assert len(tauzs) == elements_count
                    for i in range(elements_count):
                        l = f.readline()
                        ns = l.split()
                        z.Faces.append(Face(z.Nodes[int(ns[0]) - 1],
                                            z.Nodes[int(ns[1]) - 1],
                                            z.Nodes[int(ns[2]) - 1],
                                            htc = float(htcs[i]),
                                            beta = float(betas[i]),
                                            taux = float(tauxs[i]),
                                            tauy = float(tauys[i]),
                                            tauz = float(tauzs[i])))
                else:
                    pass

                l = f.readline()

        f.close()

    # ----------------------------------------------------------------------------------------------

    def Export(self, filename):
        """
        Стандартный экспорт.

        Arguments:
            filename -- имя файла.
        """

        fl = open(filename, 'w')
        fl.write('# EXPORT_MODE=CHECK_POINT\n')
        fl.write('TITLE="FE Surface Data ASCII"\n')
        fl.write('VARIABLES="X", "Y", "Z", "T", "Hw", "Hi", "HTC", "Beta", "TauX", "TauY", "TauZ"\n')

        for z in self.Zones:
            fl.write('ZONE T="%s"\n' % z.Title)
            fl.write('NODES=%s\n' % len(z.Nodes))
            fl.write('ELEMENTS=%s\n' % len(z.Faces))
            fl.write('DATAPACKING=BLOCK\n')
            fl.write('ZONETYPE=FETRIANGLE\n')
            fl.write('VARLOCATION=([4-11]=CELLCENTERED)\n')
            for n in z.Nodes:
                fl.write('%f ' % n.Point.X)
            fl.write('\n')
            for n in z.Nodes:
                fl.write('%f ' % n.Point.Y)
            fl.write('\n')
            for n in z.Nodes:
                fl.write('%f ' % n.Point.Z)
            fl.write('\n')
            for f in z.Faces:
                fl.write('%f ' % f.T)
            fl.write('\n')
            for f in z.Faces:
                fl.write('%f ' % f.Hw)
            fl.write('\n')
            for f in z.Faces:
                fl.write('%f ' % f.Hi)
            fl.write('\n')
            for f in z.Faces:
                fl.write('%f ' % f.HTC)
            fl.write('\n')
            for f in z.Faces:
                fl.write('%f ' % f.Beta)
            fl.write('\n')
            for f in z.Faces:
                fl.write('%f ' % f.Tau.X)
            fl.write('\n')
            for f in z.Faces:
                fl.write('%f ' % f.Tau.Y)
            fl.write('\n')
            for f in z.Faces:
                fl.write('%f ' % f.Tau.Z)
            fl.write('\n')

            # Осталось добавить линки.
            for f in z.Faces:
                fl.write('%d %d %d\n' % (f.Nodes[0].Id + 1, f.Nodes[1].Id + 1, f.Nodes[2].Id + 1))

        fl.close()

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    pass

# --------------------------------------------------------------------------------------------------
