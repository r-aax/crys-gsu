"""
GSU main functions.
"""

import random
from functools import reduce

# ======================================================================================================================


class Node:
    """
    Node of the grid.
    """

    # ------------------------------------------------------------------------------------------------------------------

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

        # Links with edges and faces.
        self.Edges = []
        self.Faces = []

    # ------------------------------------------------------------------------------------------------------------------

    def is_near(self, n):
        """
        Check if one node is near to another.
        :param n: another node
        :return: True - if nodes are near to each other, False - otherwise
        """

        return self.RoundedCoords == n.RoundedCoords

# ======================================================================================================================


class Edge:
    """
    Edge of the grid.
    """

    # ------------------------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        # Links to nodes and faces.
        self.Nodes = []
        self.Faces = []

    # ------------------------------------------------------------------------------------------------------------------

    def is_border(self):
        """
        Check if edge is border edge.
        :return: True - if edge is border edge, False - otherwise
        """

        # Border edge has only one neighbour face.
        return len(self.Faces) == 1

    # ------------------------------------------------------------------------------------------------------------------

    def is_interzones(self):
        """
        Check if edge is interzones.
        :return: True - if edge is interzones, False - otherwise
        """

        # Interzone edge has two neighbour faces from different zones.

        faces_count = len(self.Faces)

        if faces_count == 1:
            return False
        elif faces_count == 2:
            return self.Faces[0].Zone != self.Faces[1].Zone
        else:
            raise Exception('Edge cannot has {0} neighbours faces.'.format(faces_count))

    # ------------------------------------------------------------------------------------------------------------------

    def is_innerzones(self):
        """
        Check if edge is innerzones.
        :return: True - if edge is innerzones, False - otherwise
        """

        # Innerzones edge has two faces from one zone.

        faces_count = len(self.Faces)

        if faces_count == 1:
            return False
        elif faces_count == 2:
            return self.Faces[0].Zone == self.Faces[1].Zone
        else:
            raise Exception('Edge cannot has {0} neighbours faces.'.format(faces_count))


# ======================================================================================================================


class Face:
    """
    Face of the grid.
    """

    # ------------------------------------------------------------------------------------------------------------------

    def __init__(self, data):
        """
        Constructor face.
        :param data: face data (tuple)
        """

        self.Data = data

        # Links with nodes and edges.
        self.Nodes = []
        self.Edges = []

        # Link to zone (each face belongs only to one single zone).
        self.Zone = None

    # ------------------------------------------------------------------------------------------------------------------

    def get_nodes_marks_str(self):
        """
        Get string "am bm cm", where am, bm, cm - nodes marks.
        :return: string of nodes marks
        """

        return reduce(lambda x, y: x + ' ' + y, [str(node.Mark) for node in self.Nodes])

# ======================================================================================================================


class Zone:
    """
    Zone of the grid.
    """

    # ------------------------------------------------------------------------------------------------------------------

    def __init__(self, name):
        """
        Constructor.
        :param name: name of zone
        """

        self.Name = name

        # No nodes or faces in the zone yet.
        self.Nodes = []
        self.Faces = []

    # ------------------------------------------------------------------------------------------------------------------

    def nodes_count(self):
        """
        Nodes count.
        :return: nodes count
        """

        return len(self.Nodes)

    # ------------------------------------------------------------------------------------------------------------------

    def faces_count(self):
        """
        Faces count.
        :return: faces count
        """

        return len(self.Faces)

    # ------------------------------------------------------------------------------------------------------------------

    def get_nodes_data_slice_str(self, i):
        """
        Get string composed from i-th elements of data of all nodes.
        :param i: index of nodes data
        :return: composed string
        """

        i_list = ['{0:.18e}'.format(node.Data[i]) for node in self.Nodes]
        i_str = ' '.join(i_list)

        return i_str

    # ------------------------------------------------------------------------------------------------------------------

    def get_faces_data_slice_str(self, i):
        """
        Get string composed from i-th elements of data of all faces.
        :param i: index of nodes data
        :return: composed string
        """

        i_list = ['{0:.18e}'.format(face.Data[i]) for face in self.Faces]
        i_str = ' '.join(i_list)

        return i_str

    # ------------------------------------------------------------------------------------------------------------------

    def add_face(self, face):
        """
        Add face to zone (with link).
        :param face: face
        """

        self.Faces.append(face)
        face.Zone = self

    # ------------------------------------------------------------------------------------------------------------------

    def mark_nodes(self):
        """
        Mark all nodes.
        """

        for (i, node) in enumerate(self.Nodes):
            node.Mark = i + 1

# ======================================================================================================================


class Grid:
    """
    Grid (Surface Unstructured).
    """

    # ------------------------------------------------------------------------------------------------------------------

    # Right variables str for load/store grids.
    VariablesStr = '"X", "Y", "Z", "T", "Hw", "Hi", "HTC", "Beta", "TauX", "TauY", "TauZ"'

    # ------------------------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        # Empty name.
        self.Name = ''

        # Set empty sets of nodes, faces, zones.
        self.Nodes = []
        self.Edges = []
        self.Faces = []
        self.Zones = []

        # Rounded coordinates
        self.RoundedCoordsBag = set()

    # ------------------------------------------------------------------------------------------------------------------

    def print_info(self):
        """
        Print information about grid.
        """

        ec = self.edges_count()

        print('GSU: {0}'.format(self.Name))
        print('  {0} nodes, {1} edges, {2} faces, {3} zones'.format(self.nodes_count(),
                                                                    ec,
                                                                    self.faces_count(),
                                                                    self.zones_count()))

        # Edges statistics.
        border_edges_count = 0
        interzones_edges_count = 0
        innerzones_edges_count = 0
        for edge in self.Edges:
            if edge.is_border():
                border_edges_count += 1
            elif edge.is_interzones():
                interzones_edges_count += 1
            elif edge.is_innerzones():
                innerzones_edges_count += 1
            else:
                raise Exception('Unknown type of edge.')
        border_edges_p = 100.0 * border_edges_count / ec
        interzones_edges_p = 100.0 * interzones_edges_count / ec
        innerzones_edges_p = 100.0 * innerzones_edges_count / ec
        print('  edges stats: {0} border ({1:.2f}%), '
              '{2} interzones ({3:.2f}%), '
              '{4} innerzones ({5:.2f}%)'.format(border_edges_count, border_edges_p,
                                                 interzones_edges_count, interzones_edges_p,
                                                 innerzones_edges_count, innerzones_edges_p))

    # ------------------------------------------------------------------------------------------------------------------

    def nodes_count(self):
        """
        Nodes count.
        :return: nodes count
        """

        return len(self.Nodes)

    # ------------------------------------------------------------------------------------------------------------------

    def edges_count(self):
        """
        Edges count.
        :return: edges count
        """

        return len(self.Edges)

    # ------------------------------------------------------------------------------------------------------------------

    def faces_count(self):
        """
        Faces count.
        :return: faces count
        """

        return len(self.Faces)

    # ------------------------------------------------------------------------------------------------------------------

    def zones_count(self):
        """
        Zones count.
        :return: zones count
        """

        return len(self.Zones)

    # ------------------------------------------------------------------------------------------------------------------

    def mark_nodes(self):
        """
        Mark all nodes.
        """

        for (i, node) in enumerate(self.Nodes):
            node.Mark = i + 1

    # ------------------------------------------------------------------------------------------------------------------

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

    # ------------------------------------------------------------------------------------------------------------------

    def add_new_node(self, data):
        """
        Add new node with given data.
        :param data: data
        :return: node registered in self.Nodes
        """

        new_node = Node(data)
        found_node = self.find_near_node(new_node)

        if found_node is None:
            # There is no such node in the grid.
            # We have to add it.
            self.Nodes.append(new_node)
            self.RoundedCoordsBag.add(new_node.RoundedCoords)
            return new_node
        else:
            # There is already such a node in the grid.
            # Just return it.
            return found_node

    # ------------------------------------------------------------------------------------------------------------------

    def link_node_edge(node, edge):
        """
        Link node with edge.
        :param node: node
        :param edge: edge
        """

        node.Edges.append(edge)
        edge.Nodes.append(node)

    # ------------------------------------------------------------------------------------------------------------------

    def link_node_face(node, face):
        """
        Link face with node.
        :param node: node
        :param face: face
        """

        node.Faces.append(face)
        face.Nodes.append(node)

    # ------------------------------------------------------------------------------------------------------------------

    def link_edge_face(edge, face):
        """
        Link edge with face.
        :param edge: edge
        :param face: face
        """

        edge.Faces.append(face)
        face.Edges.append(edge)

    # ------------------------------------------------------------------------------------------------------------------

    def find_edge(node_a, node_b):
        """
        Find edge with given nodes.
        :param node_a: the first node
        :param node_b: the second node
        :return: edge - if it is found, None - otherwise
        """

        for edge in node_a.Edges:
            if node_b in edge.Nodes:
                return edge

        return None

    # ------------------------------------------------------------------------------------------------------------------

    def complex_link_face_node_node_edge(self, face, node_a, node_b):
        """
        Compllex link nodes with edge, and edge with face.
        :param face: face
        :param node_a: the first node
        :param node_b: th second node
        """

        # First we need to find edge.
        edge = Grid.find_edge(node_a, node_b)

        if edge is None:
            # New edge and link it.
            edge = Edge()
            self.Edges.append(edge)
            Grid.link_node_edge(node_a, edge)
            Grid.link_node_edge(node_b, edge)
            Grid.link_edge_face(edge, face)
        else:
            # Edge is already linked with nodes.
            # Link only with the face.
            Grid.link_edge_face(edge, face)

    # ------------------------------------------------------------------------------------------------------------------

    def load(self, filename):
        """
        Load grid from file.
        :param filename: file name
        """

        # Clear all objects of the grid.
        self.Nodes.clear()
        self.Edges.clear()
        self.Faces.clear()
        self.Zones.clear()

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
                        data = (d[0][i], d[1][i], d[2][i])
                        new_node = self.add_new_node(data)
                        zone.Nodes.append(new_node)

                    # Read data for faces.
                    d = []
                    for i in range(len(Grid.VariablesStr.split()) - 3):
                        line = f.readline()
                        d.append([float(xi) for xi in line.split()])
                    for i in range(faces_to_read):
                        face = Face([d[0][i], d[1][i], d[2][i], d[3][i], d[4][i], d[5][i], d[6][i], d[7][i]])
                        self.Faces.append(face)
                        zone.add_face(face)

                    # Read connectivity lists.
                    for i in range(faces_to_read):
                        line = f.readline()
                        face = zone.Faces[i]
                        nodes = [zone.Nodes[int(ss) - 1] for ss in line.split()]
                        if len(nodes) != 3:
                            raise Exception('Wrong count of face linked nodes ({0}).'.format(len(nodes)))
                        Grid.link_node_face(nodes[0], face)
                        Grid.link_node_face(nodes[1], face)
                        Grid.link_node_face(nodes[2], face)

                else:
                    raise Exception('Unexpected line : {0}.'.format(line))

                line = f.readline()
            f.close()

            # Now we need to fix rest objects links.
            for face in self.Faces:
                node_a = face.Nodes[0]
                node_b = face.Nodes[1]
                node_c = face.Nodes[2]
                self.complex_link_face_node_node_edge(face, node_a, node_b)
                self.complex_link_face_node_node_edge(face, node_a, node_c)
                self.complex_link_face_node_node_edge(face, node_b, node_c)

    # ------------------------------------------------------------------------------------------------------------------

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

    # ------------------------------------------------------------------------------------------------------------------

    def relink_nodes_to_zones(self):
        """
        Relink nodes to zones.
        """

        # Mark nodes for easy search.
        self.mark_nodes()

        for zone in self.Zones:

            # Delete old information.
            zone.Nodes.clear()

            # Bad for duplicates monitoring.
            bag = set()

            # Add nodes.
            for face in zone.Faces:
                for node in face.Nodes:
                    if node.Mark not in bag:
                        zone.Nodes.append(node)
                        bag.add(node.Mark)

    # ------------------------------------------------------------------------------------------------------------------

    def distribute_mono(self):
        """
        Create mono distribution (with 1 zone).
        """

        # Delete all zones and create one zone with name 'mono'.
        self.Zones.clear()
        zone = Zone('mono')
        self.Zones.append(zone)
        for face in self.Faces:
            zone.add_face(face)

        # Link nodes.
        self.relink_nodes_to_zones()

    # ------------------------------------------------------------------------------------------------------------------

    def distribute_random(self, count=16):
        """
        Create random distribution.
        :param count: zones count
        """

        # Delete all zones.
        # Create 'count' zones and random distribute faces between them.
        self.Zones.clear()
        for i in range(count):
            zone = Zone('random ' + str(i))
            self.Zones.append(zone)
        for face in self.Faces:
            self.Zones[random.randint(0, count - 1)].add_face(face)

        # Link nodes.
        self.relink_nodes_to_zones()

# ======================================================================================================================
