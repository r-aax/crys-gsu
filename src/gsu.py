"""
GSU main functions.
"""

import random
from functools import reduce
import utils

# ======================================================================================================================


def fun_face_cx():
    """
    Function that returns X coordinate of the face center.
    :return: function
    """

    return lambda f: f.center_point()[0]

# ----------------------------------------------------------------------------------------------------------------------


def fun_face_cy():
    """
    Function that returns Y coordinate of the face center.
    :return: function
    """

    return lambda f: f.center_point()[1]

# ----------------------------------------------------------------------------------------------------------------------


def fun_face_cz():
    """
    Function that returns Z coordinate of the face center.
    :return: function
    """

    return lambda f: f.center_point()[2]

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
        self.Mark = -1

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

        # Face mark.
        self.Mark = -1

    # ------------------------------------------------------------------------------------------------------------------

    def get_nodes_marks_str(self):
        """
        Get string "am bm cm", where am, bm, cm - nodes marks.
        :return: string of nodes marks
        """

        return reduce(lambda x, y: x + ' ' + y, [str(node.Mark + 1) for node in self.Nodes])

    # ------------------------------------------------------------------------------------------------------------------

    def center_point(self):
        """
        Get point - center of the face.
        :return: center point (tuple of coordinates)
        """

        xs = [n.Data[0] for n in self.Nodes]
        ys = [n.Data[1] for n in self.Nodes]
        zs = [n.Data[2] for n in self.Nodes]

        return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

    # ------------------------------------------------------------------------------------------------------------------

    def get_neighbour(self, edge):
        """
        Get neighbour through edge.
        :param edge: edge
        :return: neighbour
        """

        incident_faces = len(edge.Faces)

        if incident_faces == 1:
            if edge.Faces[0] != self:
                raise Exception('Error while getting face neighbour.')
            return None
        elif incident_faces == 2:
            if edge.Faces[0] == self:
                return edge.Faces[1]
            elif edge.Faces[1] == self:
                return edge.Faces[0]
            else:
                raise Exception('Error while getting face neighbour.')
        else:
            raise Exception('Wrong edge incident faces ({0}).'.format(incident_faces))

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

        # Mark.
        self.Mark = -1

        # Faces queue.
        self.FacesQueue = []

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

    def fill_queue(self, face):
        """
        Fill queue with face neighbours.
        :param face: face
        """

        for edge in face.Edges:
            nf = face.get_neighbour(edge)
            if nf is not None:
                if nf.Zone is None:
                    nf.Zone = self
                    self.FacesQueue.append(nf)

    # ------------------------------------------------------------------------------------------------------------------

    def mark_nodes(self):
        """
        Mark all nodes.
        """

        for (i, node) in enumerate(self.Nodes):
            node.Mark = i

    # ------------------------------------------------------------------------------------------------------------------

    def grow(self):
        """
        Grow zone.
        :return: 1 - if zone grows, 0 - otherwise
        """

        if not self.FacesQueue:
            # print('I can not grow : zone {0}.'.format(self.Name))
            return 0
        else:
            fc = self.FacesQueue[0]
            self.FacesQueue = self.FacesQueue[1:]
            self.add_face(fc)
            self.fill_queue(fc)
            return 1

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

    def print_info(self,
                   is_print_edges_statistics=False,
                   is_print_faces_distribution=False,
                   is_print_zones_adjacency_matrix=False):
        """
        Print information about grid.
        :param is_print_edges_statistics: flag for print edges statistics
        :param is_print_faces_distribution: flag for print faces distribution between zones
        :param is_print_zones_adjacency_matrix: flag for print zones adjacency matrix
        """

        ec = len(self.Edges)
        fc = len(self.Faces)

        print('GSU: {0}'.format(self.Name))
        print('  {0} nodes, {1} edges, '
              '{2} faces, {3} zones'.format(len(self.Nodes), ec, fc, len(self.Zones)))

        # Edges statistics.
        if is_print_edges_statistics:
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

        # Distribution faces between zones.
        if is_print_faces_distribution:
            print('  distribution faces between zones:')
            for zone in self.Zones:
                print('    {0} : {1} faces'.format(zone.Name, len(zone.Faces)))
            zones_faces_count = [len(zone.Faces) for zone in self.Zones]
            ideal_mean = len(self.Faces) / len(self.Zones)
            max_zone_faces_count = max(zones_faces_count)
            faces_distr_dev = 100.0 * (max_zone_faces_count - ideal_mean) / ideal_mean
            print('  ~ max zone faces {0}, '
                  'faces distribution deviation : {1}%'.format(max_zone_faces_count,
                                                               faces_distr_dev))

        # Distribution edges between pairs of neighbours.
        if is_print_zones_adjacency_matrix:
            m = []
            for zone in self.Zones:
                m.append([0] * len(self.Zones))
            # Calculate border lengths.
            self.mark_zones()
            for edge in self.Edges:
                if len(edge.Faces) == 2:
                    ff = edge.Faces[0]
                    sf = edge.Faces[1]
                    if ff.Zone != sf.Zone:
                        if ff.Zone is not None and sf.Zone is not None:
                            m[ff.Zone.Mark][sf.Zone.Mark] += 1
                            m[sf.Zone.Mark][ff.Zone.Mark] += 1
            # Print distribution.
            for i in range(len(self.Zones)):
                print(' '.join(['{0:5}'.format(e) for e in m[i]]))
            # Calculate stats.
            fm = utils.flatten(m)
            max_interzone_border_length = max(fm)
            print('  ~ max interzones border length : {0}'.format(max_interzone_border_length))

    # ------------------------------------------------------------------------------------------------------------------

    def mark_nodes(self):
        """
        Mark all nodes.
        """

        for (i, node) in enumerate(self.Nodes):
            node.Mark = i

    # ------------------------------------------------------------------------------------------------------------------

    def mark_faces(self):
        """
        Mark all faces.
        """

        for (i, face) in enumerate(self.Faces):
            face.Mark = i

    # ------------------------------------------------------------------------------------------------------------------

    def mark_zones(self):
        """
        Mark all zones.
        """

        for (i, zone) in enumerate(self.Zones):
            zone.Mark = i

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
                    print('LOAD: zone {0}, nodes_to_read = {1}'.format(zone_name, nodes_to_read))
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
                f.write('NODES={0}\n'.format(len(zone.Nodes)))
                f.write('ELEMENTS={0}\n'.format(len(zone.Faces)))
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

    def unlink_faces_from_zones(self):
        """
        Unlink faces from zones.
        """

        for face in self.Faces:
            face.Zone = None

    # ------------------------------------------------------------------------------------------------------------------

    def check_faces_are_linked_to_zones(self):
        """
        Check if all faces are linked to zones.
        """

        for face in self.Faces:
            if face.Zone == None:
                raise Exception('Unlinked face detected.')

    # ------------------------------------------------------------------------------------------------------------------

    def decompose_mono(self, new_name=None):
        """
        Create mono distribution (with 1 zone).
        :param new_name: grid new name
        """

        if new_name is not None:
            self.Name = new_name

        # Delete all zones and create one zone with name 'mono'.
        self.Zones.clear()
        self.unlink_faces_from_zones()
        zone = Zone('mono')
        self.Zones.append(zone)
        for face in self.Faces:
            zone.add_face(face)

        # Link nodes.
        self.relink_nodes_to_zones()
        self.check_faces_are_linked_to_zones()

    # ------------------------------------------------------------------------------------------------------------------

    def decompose_random(self, count=32, new_name=None):
        """
        Create random distribution.
        :param count: zones count
        :param new_name: grid new name
        """

        if new_name is not None:
            self.Name = new_name

        # Delete all zones.
        # Create 'count' zones and random distribute faces between them.
        self.Zones.clear()
        self.unlink_faces_from_zones()
        for i in range(count):
            zone = Zone('random ' + str(i))
            self.Zones.append(zone)
        for face in self.Faces:
            self.Zones[random.randint(0, count - 1)].add_face(face)

        # Link nodes.
        self.relink_nodes_to_zones()
        self.check_faces_are_linked_to_zones()

    # ------------------------------------------------------------------------------------------------------------------

    def decompose_linear(self, count=32, new_name=None):
        """
        Linear distribution.
        :param count: zones count
        :param new_name: grid new name
        """

        if new_name is not None:
            self.Name = new_name

        fc = len(self.Faces)
        fcpz = fc // count
        fcrm = fc % count

        # Delete all zones.
        # Create 'count' zones and random distribute faces between them.
        self.Zones.clear()
        self.unlink_faces_from_zones()
        for i in range(count):
            zone = Zone('linear ' + str(i))
            self.Zones.append(zone)

        # Distribute faces accurately.
        cur_face_i = 0
        for (i, zone) in enumerate(self.Zones):
            faces_to_add = fcpz
            if i < fcrm:
                faces_to_add += 1
            for j in range(cur_face_i, cur_face_i + faces_to_add):
                zone.add_face(self.Faces[j])
            cur_face_i += faces_to_add
        if cur_face_i != fc:
            raise Exception('Wrong linera distribution mathematics.')

        # Link nodes.
        self.relink_nodes_to_zones()
        self.check_faces_are_linked_to_zones()

    # ------------------------------------------------------------------------------------------------------------------

    def decompose_rgrow(self, count=32, new_name=None):
        """
        Distribution with random grow.
        :param count: count of zones.
        :param new_name: grid new name
        """

        if new_name is not None:
            self.Name = new_name

        # Delete all zones and links.
        self.Zones.clear()
        self.unlink_faces_from_zones()
        for i in range(count):
            zone = Zone('rgrow ' + str(i))
            self.Zones.append(zone)

        # Initial faces.
        self.mark_faces()
        self.mark_zones()
        for zone in self.Zones:
            rf = self.Faces[random.randint(0, len(self.Faces) - 1)]
            while rf.Zone is not None:
                rf = self.Faces[random.randint(0, len(self.Faces) - 1)]
            zone.add_face(rf)
            zone.FacesQueue = []
            zone.fill_queue(rf)

        # Walk and add to faces to queues.
        total_faces = count
        while total_faces < len(self.Faces):
            for i in range(3 * count):
                zi = i % count
                total_faces += self.Zones[zi].grow()

        # Links nodes.
        self.relink_nodes_to_zones()
        self.check_faces_are_linked_to_zones()

    # ------------------------------------------------------------------------------------------------------------------

    def split_zone_metric(self, zone, zl, zr, fun):
        """
        Split zone, calculate metric and roll-back.
        :param zone: zone
        :param zl: left child zone
        :param zr: right child zone
        :param fun: function for sign extraction
        :return: metric
        """

        # Now all faces are in zone.
        signs = [fun(f) for f in zone.Faces]
        signs.sort()
        blade = signs[len(signs) // 2]

        # Distribute.
        for f in zone.Faces:
            if fun(f) < blade:
                f.Zone = zl
            else:
                f.Zone = zr

        # Metric.
        r = 0
        for f in zone.Faces:
            for e in f.Edges:
                if len(e.Faces) == 2:
                    ff = e.Faces[0]
                    sf = e.Faces[1]
                    if (ff.Zone == zl and sf.Zone == zr) or (ff.Zone == zr and sf.Zone == zl):
                        r += 1

        # Roll back.
        for f in zone.Faces:
            f.Zone = zone

        return r

    # ------------------------------------------------------------------------------------------------------------------

    def split_zone(self, zone, extract_signs_funs):
        """
        Split zone.
        :param zone: zone
        :param extract_signs_funs: list of functions for extraction.
        """

        # Technical actions.
        zl = Zone(zone.Name + 'l')
        zr = Zone(zone.Name + 'r')
        self.Zones.remove(zone)
        self.Zones.append(zl)
        self.Zones.append(zr)

        # Explore metrics.
        metrics = [self.split_zone_metric(zone, zl, zr, extract_signs_funs[i])
                   for i in range(len(extract_signs_funs))]
        spl_index = metrics.index(min(metrics))

        # Split.
        signs = [extract_signs_funs[spl_index](f) for f in zone.Faces]
        signs.sort()
        blade = signs[len(signs) // 2]

        for face in zone.Faces:
            if extract_signs_funs[spl_index](face) < blade:
                zl.add_face(face)
            else:
                zr.add_face(face)

    # ------------------------------------------------------------------------------------------------------------------

    def decompose_hierarchical(self, extract_signs_funs, levels=6, new_name=None):
        """
        Hierarchical distribution with given numbers of levels.
        :param extract_signs_funs: list of functions for signs extraction
        :param levels: levels count
        :param new_name: grid new name
        """

        if new_name is not None:
            self.Name = new_name

        # Delete all zones and links.
        self.Zones.clear()
        self.unlink_faces_from_zones()

        # Check levels.
        if levels < 1:
            raise Exception('It must be at least 1 level.')

        zone = Zone('h')
        self.Zones.append(zone)
        for face in self.Faces:
            zone.add_face(face)

        for li in range(levels - 1):
            c = len(self.Zones)
            for zi in range(c):
                nm = self.Zones[0].Name
                # print('split zone {0} -> {1}, {2}.'.format(nm, nm + 'l', nm + 'r'))
                self.split_zone(self.Zones[0], extract_signs_funs)

        # Links nodes.
        self.relink_nodes_to_zones()
        self.check_faces_are_linked_to_zones()

# ======================================================================================================================
