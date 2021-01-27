"""
GSU main functions.
"""

from functools import reduce


class Node:
    """
    Node of the grid.
    """

    def __init__(self, data, digits=10):
        """
        Constructor node.
        :param data: node data (tuple)
        :param digits: digits count for rounding coordinates
        """

        self.Data = data

        # Variable for any purpose marking.
        self.Mark = 0

        # Rounded coordinates for registration in set.
        self.RoundedCoords = round(data[0], digits), round(data[1], digits), round(data[2], digits)

    def is_near(self, n):
        """
        Check if one node is near to another.
        :param n: another node
        :return: True - if nodes are near to each other, False - otherwise
        """
        return self.RoundedCoords == n.RoundedCoords


class Face:
    """
    Face of the grid.
    """

    def __init__(self, data):
        """
        Constructor face.
        :param data: face data (tuple)
        """

        self.Data = data

        # Empty nodes links
        self.Nodes = []

    def get_nodes_marks_str(self):
        """
        Get string "am bm cm", where am, bm, cm - nodes marks.
        :return: string of nodes marks
        """

        return reduce(lambda x, y: x + ' ' + y, [str(node.Mark) for node in self.Nodes])


class Zone:
    """
    Zone of the grid.
    """

    def __init__(self, name):
        """
        Constructor.
        :param name: name of zone
        """

        self.Name = name

        # No nodes or faces in the zone yet.
        self.Nodes = []
        self.Faces = []

    def nodes_count(self):
        """
        Nodes count.
        :return: nodes count
        """

        return len(self.Nodes)

    def faces_count(self):
        """
        Faces count.
        :return: faces count
        """

        return len(self.Faces)

    def get_nodes_data_slice_str(self, i):
        """
        Get string composed from i-th elements of data of all nodes.
        :param i: index of nodes data
        :return: composed string
        """

        i_list = ['{0:.18e}'.format(node.Data[i]) for node in self.Nodes]
        i_str = ' '.join(i_list)

        return i_str

    def get_faces_data_slice_str(self, i):
        """
        Get string composed from i-th elements of data of all faces.
        :param i: index of nodes data
        :return: composed string
        """

        i_list = ['{0:.18e}'.format(face.Data[i]) for face in self.Faces]
        i_str = ' '.join(i_list)

        return i_str

    def mark_nodes(self):
        """
        Mark all nodes.
        """

        for (i, node) in enumerate(self.Nodes):
            node.Mark = i + 1


class Grid:
    """
    Grid (Surface Unstructured).
    """

    # Right variables str for load/store grids.
    VariablesStr = '"X", "Y", "Z", "T", "Hw", "Hi", "HTC", "Beta", "TauX", "TauY", "TauZ"'

    def __init__(self):
        """
        Constructor.
        """

        # Empty name.
        self.Name = ''

        # Set empty sets of nodes, faces, zones.
        self.Nodes = []
        self.Faces = []
        self.Zones = []

        # Rounded coordinates
        self.RoundedCoordsBag = set()

    def print_info(self):
        """
        Print information about grid.
        """

        print('GSU: {0}'.format(self.Name))
        print('  {0} nodes, {1} faces, {2} zones'.format(self.nodes_count(),
                                                         self.faces_count(),
                                                         self.zones_count()))

    def nodes_count(self):
        """
        Nodes count.
        :return: nodes count
        """

        return len(self.Nodes)

    def faces_count(self):
        """
        Faces count.
        :return: faces count
        """

        return len(self.Faces)

    def zones_count(self):
        """
        Zones count.
        :return: zones count
        """

        return len(self.Zones)

    def find_near_node(self, n):
        """
        Find in grid nodes collection node that is near to a given node.
        :param n: node to check
        :return: near node from grid nodes collection
        """

        # First try to find  in bag.
        if n.RoundedCoords in self.RoundedCoordsBag:
            for node in self.Nodes:
                if node.is_near(n):
                    return node
            raise Exception('We expect to find node with coordinates {0} in the grid'.format(n.RoundedCoords))

        return None

    def load(self, filename):
        """
        Load grid from file.
        :param filename: file name
        """

        # Clear all objects of the grid.
        self.Nodes.clear()
        self.Faces.clear()
        self.Zones.clear()

        # Current processed zone.
        zone = None

        # Nodes and faces to read.
        nodes_to_read = 0
        faces_to_read = 0

        # Open file and try to load it line by line.
        with open(filename, 'r') as f:
            line = f.readline()
            while line:

                if '# EXPORT MODE:' in line:

                    # Head of grid.
                    mode_line = line
                    title_line = f.readline()
                    variables_line = f.readline()

                    # Parse all and check.
                    mode = mode_line.split()[-1]
                    if mode != 'CHECK_POINT':
                        raise Exception('The loaded grid must be in CHECK_POINT mode '
                                        '({0} mode is detected).'.format(mode))
                    if 'TITLE=' not in title_line:
                        raise Exception('Wrong title line ({0}).'.format(title_line))
                    self.Name = title_line.split('=')[1][1:-2]
                    if 'VARIABLES=' not in variables_line:
                        raise Exception('Wrong variables line ({0}).'.format(variables_line))
                    variables_str = variables_line.split('=')[-1][:-1]
                    if variables_str != Grid.VariablesStr:
                        raise Exception('The loaded variables str must be {0} '
                                        '({1} variables str is detected)'.format(Grid.VariablesStr, variables_str))

                elif 'ZONE T=' in line:

                    # Create new zone.
                    zone_name = line.split('=')[-1][1:-2]
                    zone = Zone(zone_name)
                    self.Zones.append(zone)

                    # Read count of nodes and faces to read.
                    nodes_line = f.readline()
                    faces_line = f.readline()
                    packing_line = f.readline()
                    zonetype_line = f.readline()
                    varlocation_line = f.readline()
                    if 'NODES=' not in nodes_line:
                        raise Exception('Wrong nodes line ({0}).'.format(nodes_line))
                    if 'ELEMENTS=' not in faces_line:
                        raise Exception('Wrong faces line ({0}).'.format(faces_line))
                    if 'DATAPACKING=BLOCK' != packing_line[:-1]:
                        raise Exception('Wrong packing line ({0}).'.format(packing_line))
                    if 'ZONETYPE=FETRIANGLE' != zonetype_line[:-1]:
                        raise Exception('Wrong zonetype line ({0}).'.format(zonetype_line))
                    if 'VARLOCATION=([4-11]=CELLCENTERED)' != varlocation_line[:-1]:
                        raise Exception('Wrong varlocation line ({0}).'.format(varlocation_line))
                    nodes_to_read = int(nodes_line.split('=')[-1][:-1])
                    faces_to_read = int(faces_line.split('=')[-1][:-1])

                    # Read data for nodes.
                    d = []
                    for i in range(3):
                        line = f.readline()
                        d.append([float(xi) for xi in line.split()])
                    for i in range(nodes_to_read):
                        new_node = Node([d[0][i], d[1][i], d[2][i]])

                        # First we try find such node in grid nodes collection.
                        # Maybe it came from another zone already.
                        reference_node = self.find_near_node(new_node)
                        if reference_node is None:
                            # Need to isolate this.
                            self.Nodes.append(new_node)
                            self.RoundedCoordsBag.add(new_node.RoundedCoords)
                            reference_node = new_node

                        zone.Nodes.append(reference_node)

                    # Read data for faces.
                    d = []
                    for i in range(len(Grid.VariablesStr.split()) - 3):
                        line = f.readline()
                        d.append([float(xi) for xi in line.split()])
                    for i in range(faces_to_read):
                        face = Face([d[0][i], d[1][i], d[2][i], d[3][i], d[4][i], d[5][i], d[6][i], d[7][i]])
                        self.Faces.append(face)
                        zone.Faces.append(face)

                    # Read connectivity lists.
                    for i in range(faces_to_read):
                        line = f.readline()
                        zone.Faces[i].Nodes = [zone.Nodes[int(ss) - 1] for ss in line.split()]

                else:
                    raise Exception('Unexpected line : {0}.'.format(line))

                line = f.readline()
            f.close()

    def store(self, filename):
        """
        Store grid to file.
        :param filename: file name
        """

        with open(filename, 'w', newline='\n') as f:

            # Store head.
            f.write('# EXPORT MODE: CHECK_POINT\n')
            f.write('TITLE="{0}"\n'.format(self.Name))
            f.write('VARIABLES={0}\n'.format(Grid.VariablesStr))

            # Store zones.
            for zone in self.Zones:

                # Store zone head.
                f.write('ZONE T="{0}"\n'.format(zone.Name))
                f.write('NODES={0}\n'.format(zone.nodes_count()))
                f.write('ELEMENTS={0}\n'.format(zone.faces_count()))
                f.write('DATAPACKING=BLOCK\n')
                f.write('ZONETYPE=FETRIANGLE\n')
                f.write('VARLOCATION=([4-11]=CELLCENTERED)\n')

                # Write first 3 data items (X, Y, Z coordinates).
                for i in range(3):
                    f.write(zone.get_nodes_data_slice_str(i) + ' \n')

                # Write rest faces data items.
                for i in range(len(Grid.VariablesStr.split()) - 3):
                    f.write(zone.get_faces_data_slice_str(i) + ' \n')

                # Write connectivity lists.
                zone.mark_nodes()
                for face in zone.Faces:
                    f.write(face.get_nodes_marks_str() + '\n')

            f.close()
