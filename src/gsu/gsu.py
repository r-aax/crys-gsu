"""
GSU main functions.
"""

import random
import math
import utils
from geom.vect import Vect
from geom.segment import Segment
from geom.triangle import Triangle
from geom.trajectory import Trajectory
from geom.box import Box
from gsu.node import Node
from gsu.edge import Edge
from gsu.face import Face
from gsu.zone import Zone
from gsu.zones_adjacency_matrix import ZonesAdjacencyMatrix

# ==================================================================================================


def mean_nodes_point(ns):
    """
    Get mean point of nodes.

    :param ns: nodes
    :return: mean point
    """

    xs = [n.P[0] for n in ns]
    ys = [n.P[1] for n in ns]
    zs = [n.P[2] for n in ns]

    return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

# --------------------------------------------------------------------------------------------------


def fun_face_cx():
    """
    Function that returns X coordinate of the face center.

    :return: function
    """

    return lambda f: mean_nodes_point(f.Nodes)[0]

# --------------------------------------------------------------------------------------------------


def fun_face_cy():
    """
    Function that returns Y coordinate of the face center.

    :return: function
    """

    return lambda f: mean_nodes_point(f.Nodes)[1]

# --------------------------------------------------------------------------------------------------


def fun_face_cz():
    """
    Function that returns Z coordinate of the face center.

    :return: function
    """

    return lambda f: mean_nodes_point(f.Nodes)[2]

# ==================================================================================================


class Grid:
    """
    Grid (Surface Unstructured).
    """

    # ----------------------------------------------------------------------------------------------

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

        self.number_of_border_nodes = 0

    # ----------------------------------------------------------------------------------------------

    def clear(self):
        """
        Clear grid.
        """

        self.Nodes.clear()
        self.Edges.clear()
        self.Faces.clear()
        self.Zones.clear()
        self.RoundedCoordsBag.clear()

    # ----------------------------------------------------------------------------------------------

    def zones_count(self):
        """
        Get zones count.

        :return: zones count
        """

        return len(self.Zones)

    # ----------------------------------------------------------------------------------------------

    def nodes_count(self):
        """
        Get count of nodes.

        :return: nodes count
        """

        return len(self.Nodes)

    # ----------------------------------------------------------------------------------------------

    def edges_count(self):
        """
        Get count of edges.

        :return: edges count
        """

        return len(self.Edges)

    # ----------------------------------------------------------------------------------------------

    def faces_count(self):
        """
        Get count of faces.

        :return: faces count.
        """

        return len(self.Faces)

    # ----------------------------------------------------------------------------------------------

    def faces_count_in_zones(self):
        """
        Count of faces that are placed in zones.

        :return: total faces count in zones
        """

        return sum([z.faces_count() for z in self.Zones])

    # ----------------------------------------------------------------------------------------------

    def faces_count_in_fixed_zones(self):
        """
        Count of faces in fixed zones.

        :return: faces count in fixed zones
        """

        return sum([z.faces_count() for z in self.Zones if z.IsFixed])

    # ----------------------------------------------------------------------------------------------

    def get_triangles_list(self):
        """
        Get triangles list.

        :return: Triangles list.
        """

        return [f.get_triangle() for f in self.Faces]

    # ----------------------------------------------------------------------------------------------

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

        ec = self.edges_count()
        fc = self.faces_count()
        zc = len(self.Zones)

        print('GSU: {0}'.format(self.Name))
        print('  {0} nodes, {1} edges, '
              '{2} faces, {3} zones'.format(len(self.Nodes), ec, fc, zc))

        # Zones adjacency matrix.
        zam = ZonesAdjacencyMatrix(self.Edges, self.Zones)

        # Edges statistics.
        if is_print_edges_statistics:
            print('  ' + zam.edges_statistics_string())

        # Distribution faces between zones.
        if is_print_faces_distribution:
            print('  distribution faces between zones:')
            for zone in self.Zones:
                print('    {0} : {1} faces, '
                      '{2} outer edges, QF = {3:.3f}'.format(zone.Name,
                                                             zone.faces_count(),
                                                             zone.outer_edges_count(),
                                                             zone.zone_quality_factor()))
            zones_faces_count = [len(zone.Faces) for zone in self.Zones]
            zones_quality_factors = [zone.zone_quality_factor() for zone in self.Zones]
            ideal_mean = len(self.Faces) / len(self.Zones)
            max_zone_faces_count = max(zones_faces_count)
            faces_distr_dev = 100.0 * (max_zone_faces_count - ideal_mean) / ideal_mean
            print('  ~ max zone faces {0}, '
                  'faces distribution deviation : {1:.2f}%'.format(max_zone_faces_count,
                                                                   faces_distr_dev))
            print('  ~ min/max quality factor : '
                  '{0:.3f}, {1:.3f}'.format(min(zones_quality_factors),
                                            max(zones_quality_factors)))

        # Distribution edges between pairs of neighbours.
        if is_print_zones_adjacency_matrix:
            for i in range(zc + 1):
                print(' '.join(['{0:5}'.format(e) for e in zam.M[i]]))
            print('  ~ max cross-zones border length : {0}'.format(zam.max_cross_border_len()))

    # ----------------------------------------------------------------------------------------------

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
            raise Exception('We expect to find node ' \
                            'with coordinates {0} in the grid'.format(n.RoundedCoords))

        return None

    # ----------------------------------------------------------------------------------------------

    def add_node(self, n, is_merge_same_nodes):
        """
        Add node.

        :param n: node
        :param is_merge_same_nodes: flag merge same nodes
        :return: node registered in self.Nodes
        """

        found_node = self.find_near_node(n)

        if (found_node is None) or (not is_merge_same_nodes):
            # There is no such node in the grid.
            # We have to add it.
            n.GloId = len(self.Nodes)
            self.Nodes.append(n)
            self.RoundedCoordsBag.add(n.RoundedCoords)
            return n
        else:
            # There is already such a node in the grid.
            # Just return it.
            return found_node

    # ----------------------------------------------------------------------------------------------

    def add_edge(self, e, global_id=None):
        """
        Add edge to grid.

        :param e: edge
        :return: added edge
        """

        # Just add edge with global id correction.
        if not global_id is None:
            e.GloId = global_id
        else:
            e.GloId = len(self.Edges)
        self.Edges.append(e)

        return e

    # ----------------------------------------------------------------------------------------------

    def add_face(self, f, global_id=None):
        """
        Add face.

        :param f: face
        :return: added face
        """

        # Just correct global id.
        if not global_id is None:
            f.GloId = global_id
        else:
            f.GloId = len(self.Faces)
        # and add
        self.Faces.append(f)

        return f

    # ----------------------------------------------------------------------------------------------

    def set_zones_ids(self):
        """
        Set zones identifiers.
        """

        for (i, z) in enumerate(self.Zones):
            z.Id = i

    # ----------------------------------------------------------------------------------------------

    def reset_zones_ids(self):
        """
        Reset zones identifiers.
        """

        for z in self.Zones:
            z.Id = -1

    # ----------------------------------------------------------------------------------------------

    def link_node_edge(node, edge):
        """
        Link node with edge.

        :param node: node
        :param edge: edge
        """

        assert (type(node) is Node)
        assert (type(edge) is Edge)
        node.Edges.append(edge)
        edge.Nodes.append(node)

    # ----------------------------------------------------------------------------------------------

    def link_node_face(node, face):
        """
        Link face with node.

        :param node: node
        :param face: face
        """

        assert(type(node) is Node)
        assert (type(face) is Face)
        node.Faces.append(face)
        face.Nodes.append(node)

    # ----------------------------------------------------------------------------------------------

    def link_edge_face(edge, face):
        """
        Link edge with face.

        :param edge: edge
        :param face: face
        """

        assert (type(edge) is Edge)
        assert (type(face) is Face)
        # Check if it is enable to link the face with the edge.
        if len(edge.Faces) == 2:
            raise Exception('Too many faces linking with this edge ({0} - {1},'
                            'GloId = {2})'.format(edge.Nodes[0].P,
                                                  edge.Nodes[1].P,
                                                  edge.GloId))

        edge.Faces.append(face)
        face.Edges.append(edge)

    # ----------------------------------------------------------------------------------------------

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

    # ----------------------------------------------------------------------------------------------

    def complex_link_face_node_node_edge(self, face, node_a, node_b):
        """
        Complex link nodes with edge, and edge with face.

        :param face: face
        :param node_a: the first node
        :param node_b: th second node
        """

        # First we need to find edge.
        edge = Grid.find_edge(node_a, node_b)

        if edge is None:
            # New edge and link it.
            edge = Edge()
            self.add_edge(edge)
            Grid.link_node_edge(node_a, edge)
            Grid.link_node_edge(node_b, edge)
            Grid.link_edge_face(edge, face)
        else:
            # Edge is already linked with nodes.
            # Link only with the face.
            Grid.link_edge_face(edge, face)

    # ----------------------------------------------------------------------------------------------

    def bfs_path_connectivity_component(self, start, pred):
        """
        Connectivity component of BFS path from given face.

        :param start: start face
        :param pred: predicate for faces
        :return: path for one connectivity component (list of faces)
        """

        # This function uses bfs_mark field for faces.

        if start.bfs_mark or not pred(start):
            # Start face does not satisfy conditions.
            return []

        # Initiate path.
        start.bfs_mark = True
        p = [start]
        i = 0

        # Walk.
        while i < len(p):
            f = p[i]
            for e in f.Edges:
                if not e.is_border():
                    f2 = f.get_neighbour(e)
                    if not f2.bfs_mark and pred(f2):
                        f2.bfs_mark = True
                        p.append(f2)
            i = i + 1

        return p

    # ----------------------------------------------------------------------------------------------

    def get_no_bfs_mark_face(self, pred):
        """
        Get first no bfs mark face.

        :param pred: predicate for face
        :return: first face with false bfs mark
        """

        for f in self.Faces:
            if not f.bfs_mark and pred(f):
                return f

        return None

    # ----------------------------------------------------------------------------------------------

    def bfs_path(self, start, pred):
        """
        Get BFS path.

        :param start: start face
        :param pred: predicate for faces
        :return: path for whole grid (list of faces)
        """

        # This function uses bfs_mark field for faces.

        # Reset all faces bfs marks.
        for f in self.Faces:
            f.bfs_mark = False

        p = []

        # Walk while we can.
        while True:
            p = p + self.bfs_path_connectivity_component(start, pred)
            start = self.get_no_bfs_mark_face(pred)
            if start is None:
                return p

    # ----------------------------------------------------------------------------------------------

    def load(self, filename,
             is_merge_same_nodes=True):
        """
        Load grid from file.

        :param filename: file name
        :param is_merge_same_nodes: merge same nodes
        """

        variables = []
        face_variables = []
        face_variables_count = 0

        # Clear all objects of the grid.
        self.clear()

        # Open file and try to load it line by line.
        with open(filename, 'r') as f:
            line = f.readline()
            while line:

                if line[0] == '#':
                    # Comment.
                    pass
                elif 'TITLE=' in line:
                    self.Name = line.split('=')[1][1:-2]
                elif 'VARIABLES=' in line:
                    variables_str = line.split('=')[-1][:-1]
                    variables = variables_str.replace('"', '').replace(',', '').split()
                    face_variables = variables[3:]
                    face_variables_count = len(face_variables)

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
                    right_varlocation_line = 'VARLOCATION=' \
                                             '([4-{0}]=CELLCENTERED)'.format(len(variables))
                    if right_varlocation_line != varlocation_line[:-1]:
                        raise Exception('Wrong varlocation line ({0}). '
                                        'Right value is {1}'.format(varlocation_line,
                                                                    right_varlocation_line))
                    nodes_to_read = int(nodes_line.split('=')[-1][:-1])
                    # print('LOAD: zone {0}, nodes_to_read = {1}'.format(zone_name, nodes_to_read))
                    faces_to_read = int(faces_line.split('=')[-1][:-1])

                    # Read data for nodes.
                    c = []
                    for i in range(3):
                        line = f.readline()
                        c.append([float(xi) for xi in line.split()])
                    for i in range(nodes_to_read):
                        p = Vect(c[0][i], c[1][i], c[2][i])
                        node = Node(p)
                        node = self.add_node(node, is_merge_same_nodes)
                        zone.add_node(node)

                    # Read data for faces.
                    d = []
                    for i in range(face_variables_count):
                        line = f.readline()
                        d.append([float(xi) for xi in line.split()])
                    for i in range(faces_to_read):
                        face = Face(face_variables,
                                    [d[j][i] for j in range(face_variables_count)])
                        self.add_face(face)
                        zone.add_face(face)

                    # Read connectivity lists.
                    for i in range(faces_to_read):
                        line = f.readline()
                        face = zone.Faces[i]
                        nodes = [zone.Nodes[int(ss) - 1] for ss in line.split()]
                        if len(nodes) != 3:
                            raise Exception('Wrong count of ' \
                                            'face linked nodes ({0}).'.format(len(nodes)))
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

            # Relink.
            self.link_nodes_and_edges_to_zones()

    # ----------------------------------------------------------------------------------------------

    def convert_grid_stall_to_check_point(self):
        """
        Remove all Stall fields.
        """

        # Delete fields Stall, StallD, StallVX, StallVY, StallVZ.
        for f in self.Faces:
            for d in ['Stall', 'StallD', 'StallVX', 'StallVY', 'StallVZ']:
                f.Data.pop(d)

    # ----------------------------------------------------------------------------------------------

    def clean_t_hw_hi(self):
        """
        Clean T, Hw, Hi fields.
        """

        for f in self.Faces:
            f['T'] = 0.0
            f['Hw'] = 0.0
            f['Hi'] = 0.0

    # ----------------------------------------------------------------------------------------------

    def store(self, filename):
        """
        Store grid to file.

        :param filename: file name
        """

        variables = ['X', 'Y', 'Z'] + list(self.Faces[0].Data.keys())

        with open(filename, 'w', newline='\n') as f:

            # Store head.
            f.write('# crys-gsu\n')
            #if self.Name != '':
            # TITLE is not needed if empty, nut remesher craches if there is no title.
            f.write('TITLE="{0}"\n'.format(self.Name))
            f.write('VARIABLES={0}\n'.format(', '.join(['"{0}"'.format(k) for k in variables])))

            # Additional structure for calculating local identifiers
            # of the nodes for connectivity lists storing.
            loc_ids = [-1] * len(self.Nodes)

            # Store zones.
            for zone in self.Zones:

                # Store zone head.
                f.write('ZONE T="{0}"\n'.format(zone.Name))
                f.write('NODES={0}\n'.format(len(zone.Nodes)))
                f.write('ELEMENTS={0}\n'.format(len(zone.Faces)))
                f.write('DATAPACKING=BLOCK\n')
                f.write('ZONETYPE=FETRIANGLE\n')
                f.write('VARLOCATION=([4-{0}]=CELLCENTERED)\n'.format(len(variables)))

                # Write first 3 data items (X, Y, Z coordinates).
                for i in range(3):
                    f.write(zone.get_nodes_coord_slice_str(i) + ' \n')

                # Write rest faces data items.
                for e in variables[3:]:
                    f.write(zone.get_faces_data_slice_str(e) + ' \n')

                # Write connectivity lists.
                for i, node in enumerate(zone.Nodes):
                    loc_ids[node.GloId] = i
                for face in zone.Faces:
                    f.write(' '.join([str(loc_ids[n.GloId] + 1) for n in face.Nodes]) + '\n')

            f.close()

    # ----------------------------------------------------------------------------------------------

    def store_mpi(self, filename_base, ts, sf='.cry'):
        """
        Store grid for mpi program.
        As many processes count as zones count.

        :param filename_base: Base of filename.
        :param ts:            Timestamp string.
        :param sf:            Suffixes of files.
        """

        # Fixed set of variables for swim.
        variables = ['X', 'Y', 'Z', 'GloId',
                     'T', 'Hw', 'Hi', 'HTC', 'Beta', 'MImp2', 'Vd2',
                     'TauX', 'TauY', 'TauZ', 'RecoveryFactor']

        zam = ZonesAdjacencyMatrix(self.Edges, self.Zones)

        # Local identifiers.
        loc_nodes_ids = [-1] * len(self.Nodes)
        loc_edges_ids = [-1] * len(self.Edges)

        for (zi, z) in enumerate(self.Zones):
            with open('{0}_{1:05d}_{2}{3}'.format(filename_base, zi, ts, sf),
                      'w', newline='\n') as file:

                nc = len(z.Nodes)
                ec = len(z.Edges)
                fc = len(z.Faces)

                # Write head information.
                file.write('TITLE="{0}"\n'.format(z.Name))
                variables_str = ', '.join(['"{0}"'.format(x) for x in variables])
                file.write('VARIABLES={0}\n'.format(variables_str))
                file.write('MPI={0}\n'.format(zi))
                file.write('NODES={0}\n'.format(nc))
                file.write('EDGES={0}\n'.format(ec))
                file.write('CROSS-EDGES={0}\n'.format(zam.zone_cross_edges_count(zi)))
                file.write('FACES={0}\n'.format(fc))

                # Write nodes and faces data.
                file.write('NODES COORDINATES:\n')
                for i in range(3):
                    file.write(z.get_nodes_coord_slice_str(i) + ' \n')
                file.write('FACES DATA:\n')
                file.write(z.get_faces_global_ids_slice_str() + '\n')
                for e in variables[4:]:
                    file.write(z.get_faces_data_slice_str(e) + ' \n')

                # Write connectivity lists.
                for i, node in enumerate(z.Nodes):
                    loc_nodes_ids[node.GloId] = i
                for i, edge in enumerate(z.Edges):
                    loc_edges_ids[edge.GloId] = i
                file.write('FACE-NODE CONNECTIVITY LIST:\n')
                for f in z.Faces:
                    file.write(' '.join([str(loc_nodes_ids[n.GloId]) for n in f.Nodes]) + '\n')
                file.write('FACE-EDGE CONNECTIVITY LIST:\n')
                for f in z.Faces:
                    file.write(' '.join([str(loc_edges_ids[e.GloId]) for e in f.Edges]) + '\n')
                file.write('EDGE-NODE CONNECTIVITY LIST:\n')
                for e in z.Edges:
                    file.write(' '.join([str(loc_nodes_ids[n.GloId]) for n in e.Nodes]) + '\n')

                # Write MPI-borders.
                file.write('MPI-BORDERS={0}\n'.format(zam.zone_cross_borders_count(zi)))
                line = zam.zone_cross_edges_array(zi)
                for (li, le) in enumerate(line):
                    if le > 0:
                        # Border between zones with indices zi and li.
                        lz = self.Zones[li]
                        # Border between zones z and lz.
                        file.write('[{0}] '.format(li))
                        for e in z.Edges:
                            if e.is_connect_zones(z, lz):
                                file.write('{0} '.format(loc_edges_ids[e.GloId]))
                        file.write('\n')

                file.close()

    # ----------------------------------------------------------------------------------------------

    def link_nodes_and_edges_to_zones(self):
        """
        Link nodes and edges to zones.
        """

        # Clear old information.
        for z in self.Zones:
            z.Nodes.clear()
            z.Edges.clear()

        self.set_zones_ids()

        # Add nodes.
        for n in self.Nodes:
            zids = list(set([f.Zone.Id for f in n.Faces]))
            for zid in zids:
                self.Zones[zid].add_node(n)

        # Add edges.
        for e in self.Edges:
            zids = list(set([f.Zone.Id for f in e.Faces]))
            for zid in zids:
                self.Zones[zid].add_edge(e)

        self.reset_zones_ids()

    # ----------------------------------------------------------------------------------------------

    def unlink_faces_from_zones(self):
        """
        Unlink faces from zones.
        """

        for face in self.Faces:
            if not face.Zone.IsFixed:
                face.Zone = None

    # ----------------------------------------------------------------------------------------------

    def check_faces_are_linked_to_zones(self):
        """
        Check if all faces are linked to zones.
        """

        for face in self.Faces:
            if face.Zone is None:
                raise Exception('Unlinked face detected.')

    # ----------------------------------------------------------------------------------------------

    def post_decompose(self):
        """
        Post actions after decompose.
        """

        self.link_nodes_and_edges_to_zones()
        self.check_faces_are_linked_to_zones()
        self.set_zones_ids()

    # ----------------------------------------------------------------------------------------------

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

        self.post_decompose()

    # ----------------------------------------------------------------------------------------------

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

        self.post_decompose()

    # ----------------------------------------------------------------------------------------------

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

        self.post_decompose()

    # ----------------------------------------------------------------------------------------------

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
        zlc = 0
        zrc = 0
        for f in zone.Faces:
            if fun(f) < blade:
                f.Zone = zl
                zlc += 1
            else:
                f.Zone = zr
                zrc += 1

        # Prevent empty zone.
        if (zlc == 0) or (zrc == 0):
            return math.inf

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

    # ----------------------------------------------------------------------------------------------

    def split_zone(self, zone, extract_signs_funs):
        """
        Split zone.

        :param zone: zone
        :param extract_signs_funs: list of functions for extraction.
        """

        # Do not split zone if it is fixed.
        if zone.IsFixed:
            self.Zones.remove(zone)
            self.Zones.append(zone)
            return

        # Do not split zone if it is impossible.
        if len(zone.Faces) < 2:
            print('Warning : not enough faces to split zone {0} ' \
                  '({1} faces).'.format(zone.Name, len(zone.Faces)))
            # Replace to the end.
            self.Zones.remove(zone)
            self.Zones.append(zone)
            return

        # Technical actions.
        zl = Zone(zone.Name + 'l')
        zr = Zone(zone.Name + 'r')
        self.Zones.remove(zone)
        self.Zones.append(zl)
        self.Zones.append(zr)

        # Explore metrics.
        metrics = [self.split_zone_metric(zone, zl, zr, extract_signs_funs[i])
                   for i in range(len(extract_signs_funs))]
        min_metric = min(metrics)

        if min_metric is math.inf:
            print('Warning : bad metrics to split zone {0}'.format(zone.Name))
            # Replace to the end.
            self.Zones.remove(zone)
            self.Zones.append(zone)
            return

        spl_index = metrics.index(min_metric)
        extract_signs_fun = extract_signs_funs[spl_index]

        # Split.
        signs = [extract_signs_fun(f) for f in zone.Faces]
        signs.sort()
        blade = signs[len(signs) // 2]

        # First we split faces into two sets: left and right and
        # then move faces sets to zl and ar zones.
        # Otherwise we have bug inside zone.Faces iterator when we
        # remove faces from zone.Faces in the loop.
        l_faces = [f for f in zone.Faces if extract_signs_fun(f) < blade]
        r_faces = [f for f in zone.Faces if extract_signs_fun(f) >= blade]
        for f in l_faces:
            zl.add_face(f)
        for f in r_faces:
            zr.add_face(f)

    # ----------------------------------------------------------------------------------------------

    def decompose_hierarchical(self, extract_signs_funs, levels=6, new_name=None, fixed_zones=[]):
        """
        Hierarchical distribution with given numbers of levels.

        :param extract_signs_funs: list of functions for signs extraction
        :param levels: levels count
        :param new_name: grid new name
        :param fixed_zones: list of fixed zones
        """

        # Mark fixed zones.
        for z in self.Zones:
            z.IsFixed = z.Name in fixed_zones

        if new_name is not None:
            self.Name = new_name

        # Delete all zones (not fixed) and links.
        self.Zones = [z for z in self.Zones if z.IsFixed]
        self.unlink_faces_from_zones()

        # Check levels.
        if levels < 1:
            raise Exception('It must be at least 1 level.')

        zone = Zone('h')
        self.Zones.append(zone)
        for face in self.Faces:
            if face.Zone is None:
                zone.add_face(face)

        for li in range(levels - 1):
            c = len(self.Zones)
            for zi in range(c):
                nm = self.Zones[0].Name
                # print('split zone {0} -> {1}, {2}.'.format(nm, nm + 'l', nm + 'r'))
                self.split_zone(self.Zones[0], extract_signs_funs)

        self.post_decompose()


    # ----------------------------------------------------------------------------------------------

    def decompose_farhat(self, count=32, new_name=None, fz_names=[]):
        """
        Create distribution based on Farhat's algorithm.

        :param count: zones count
        :param new_name: grid new name
        :param fz_names: list of fixed zones names
        """

        # Mark fixed zones.
        for z in self.Zones:
            z.IsFixed = z.Name in fz_names

        if new_name is not None:
            self.Name = new_name

        # Keep fixed zones and create additional.
        self.Zones = [z for z in self.Zones if z.IsFixed]
        fz_count = len(self.Zones)
        self.unlink_faces_from_zones()
        for i in range(count):
            zone = Zone('farhat ' + str(i))
            self.Zones.append(zone)

        #
        # Distribute faces between zones.
        #

        fc = self.faces_count() - self.faces_count_in_fixed_zones()
        zone_sizes = [fc // count + (fc % count > i) for i in range(count)]

        # Put right number of faces into each zone.
        start = self.Faces[0]
        for i in range(count):
            z = self.Zones[i + fz_count]
            zone_size = zone_sizes[i]
            path = self.bfs_path(start, lambda f: (f.Zone is None))
            for fi in range(zone_size):
                z.add_face(path[fi])
            if len(path) > zone_size:
                start = path[zone_size]

        self.post_decompose()

    # ----------------------------------------------------------------------------------------------

    def box(self):
        """
        Get box around grid (tuple with 6 values - XYZ of the left down back point
        and XYZ of the right up front point).

        :return: tuple
        """

        xs = [n.P[0] for n in self.Nodes]
        ys = [n.P[1] for n in self.Nodes]
        zs = [n.P[2] for n in self.Nodes]

        return min(xs), min(ys), min(zs), max(xs), max(ys), max(zs)

    # ----------------------------------------------------------------------------------------------

    def move_from_mean_point(self, k=0.1):
        """
        Move zones from mean point.

        :param k: factor while zones moving
        """

        gx, gy, gz = mean_nodes_point(self.Nodes)

        for z in self.Zones:
            zx, zy, zz = mean_nodes_point(z.Nodes)
            vx, vy, vz = k * (zx - gz), k * (zy - gy), k * (zz - gz)
            for n in z.Nodes:
                p = n.P
                n.P = (p[0] + vx, p[1] + vy, p[2] + vz)

    # ----------------------------------------------------------------------------------------------

    def store_faces_calc_data(self, filename):
        """
        Store calc values to file.

        :param filename: name of file
        """

        with open(filename, 'w', newline='\n') as f:
            f.write('VARIABLES="GloId", "T", "Hw", "Hi"\n')
            for face in self.Faces:
                p = mean_nodes_point(face.Nodes)
                f.write(face.get_glo_id_t_hw_hi_str() + '\n')
            f.close()

    # ----------------------------------------------------------------------------------------------

    def load_faces_calc_data(self, filename):
        """
        Load calc data from file and write it to grid.

        :param filename: Name of file.
        """

        # Read variables from the first line.
        f = open(filename, 'r')
        line = f.readline()
        variables = line.split('=')[-1].replace('"', '').replace(',', '').split()[1:]
        line = f.readline()

        # Process all rest lines.
        while line:
            ss = line.split()
            glo_id = int(ss[0])
            values = [float(s) for s in ss[1:]]
            face = self.Faces[glo_id]
            face.Data = dict(zip(variables, values))
            line = f.readline()

        f.close()

    # ----------------------------------------------------------------------------------------------

    def mark_all_fixed_nodes(self):
        """Mark border nodes as fixed."""
        for e in self.Edges:
            if len(e.Faces) == 1:
                e.border = True
                for n in e.Nodes:
                    if not n.border:
                        n.border = True
                        self.number_of_border_nodes += 1

    # ----------------------------------------------------------------------------------------------

    def del_face(self, face):
        """

        Parameters
        ----------
        face: Face

        Returns
        -------

        """

        for node in face.Nodes:
            if face in node.Faces:
                node.Faces.remove(face)
        for edge in face.Edges:
            if face in edge.Faces:
                edge.Faces.remove(face)
        face.Zone.Faces.remove(face)
        self.Faces.remove(face)

    # ----------------------------------------------------------------------------------------------

    def del_edge(self, edge):
        """

        Parameters
        ----------
        edge

        Returns
        -------

        """

        edge.Nodes[0].Edges.remove(edge)
        edge.Nodes[1].Edges.remove(edge)
        if edge.Faces:
            for f in edge.Faces:
                f.Edges.remove(edge)
        for z in self.Zones:
            if edge in z.Edges:
                z.Edges.remove(edge)
        self.Edges.remove(edge)

    # ----------------------------------------------------------------------------------------------

    def del_node(self, node):
        """

        Parameters
        ----------
        node: Node

        Returns
        -------

        """
        if node.Edges:
            for ed in node.Edges:
                if node in ed.Nodes:
                    ed.Nodes.remove(node)
        if node.Faces:
            for face in node.Faces:
                if node in face.Nodes:
                    face.Nodes.remove(node)
        node.Faces[0].Zone.Nodes.remove(node)
        self.Nodes.remove(node)
        self.RoundedCoordsBag.remove(node.RoundedCoords)

    # ----------------------------------------------------------------------------------------------

    def divide_face(self, face, p):
        """

        Atomic grid transformation. Dividing a face by a point.

        Parameters
        ----------
        face: Face
        p: Vector

        Returns
        -------

        """

        data = face.Data.copy()
        nodes = face.Nodes
        edges = face.Edges
        zone = face.Zone
        gloid = face.GloId

        self.del_face(face)

        # get new faces
        f1 = Face(list(data.keys()), list(data.values()))
        f2 = Face(list(data.keys()), list(data.values()))
        f3 = Face(list(data.keys()), list(data.values()))

        # get nodes
        n1 = nodes[0]
        n2 = nodes[1]
        n3 = nodes[2]
        n4 = Node(p)
        node = self.add_node(n4, True)
        if id(node) == id(n4):
            zone.add_node(n4)
        else:
            n4 = node

        # get edges
        e1 = Grid.find_edge(n1, n2)
        e2 = Grid.find_edge(n2, n3)
        e3 = Grid.find_edge(n1, n3)
        e4 = Edge()
        e5 = Edge()
        e6 = Edge()

        # links
        Grid.link_node_edge(n1, e4)
        Grid.link_node_edge(n2, e5)
        Grid.link_node_edge(n3, e6)
        Grid.link_node_edge(n4, e4)
        Grid.link_node_edge(n4, e5)
        Grid.link_node_edge(n4, e6)

        Grid.link_edge_face(e1, f1)
        Grid.link_edge_face(e2, f2)
        Grid.link_edge_face(e3, f3)
        Grid.link_edge_face(e4, f1)
        Grid.link_edge_face(e4, f3)
        Grid.link_edge_face(e5, f1)
        Grid.link_edge_face(e5, f2)
        Grid.link_edge_face(e6, f2)
        Grid.link_edge_face(e6, f3)

        Grid.link_node_face(n1, f1)
        Grid.link_node_face(n1, f3)
        Grid.link_node_face(n2, f1)
        Grid.link_node_face(n2, f2)
        Grid.link_node_face(n3, f2)
        Grid.link_node_face(n3, f3)
        Grid.link_node_face(n4, f1)
        Grid.link_node_face(n4, f2)
        Grid.link_node_face(n4, f3)

        # link grid with faces
        self.add_face(f1, gloid)
        self.add_face(f2)
        self.add_face(f3)

        # link grid
        self.add_edge(e4)
        self.add_edge(e5)
        self.add_edge(e6)

        # link zone
        zone.add_face(f1)
        zone.add_face(f2)
        zone.add_face(f3)
        zone.add_edge(e4)
        zone.add_edge(e5)
        zone.add_edge(e6)

    # ----------------------------------------------------------------------------------------------

    def collapse_face(self, face):
        """

        Atomic grid transformation. Tightening the grid along the borders.

        Parameters
        ----------
        face: Face

        Returns
        -------

        """

        data = face.Data.copy()
        gloid1 = face.GloId
        zone = face.Zone

        max_edge = face.max_edge_from_face()
        if len(max_edge.Faces) == 2:
            self.del_face(face)
            self.del_edge(max_edge)
            other_face = max_edge.Faces[0]
            gloid2 = other_face.GloId
            self.del_face(other_face)

            # get new faces
            f1 = Face(list(data.keys()), list(data.values()))
            f2 = Face(list(data.keys()), list(data.values()))

            # get nodes
            n1 = max_edge.Nodes[0]
            n2 = max_edge.Nodes[1]
            n3 = [i for i in face.Nodes if id(i) != id(n1) and id(i) != id(n2)][0]
            n4 = [i for i in other_face.Nodes if id(i) != id(n1) and id(i) != id(n2)][0]

            # get edges
            face.Edges.remove(max_edge)
            e1 = [edge for edge in face.Edges for node in edge.Nodes if id(node) == id(n1)][0]
            e2 = [edge for edge in other_face.Edges for node in edge.Nodes if id(node) == id(n1)][0]
            e3 = [edge for edge in face.Edges for node in edge.Nodes if id(node) == id(n2)][0]
            e4 = [edge for edge in other_face.Edges for node in edge.Nodes if id(node) == id(n2)][0]
            e5 = Edge()
            self.add_edge(e5, max_edge.GloId)
            Grid.link_node_edge(n3, e5)
            Grid.link_node_edge(n4, e5)

            # links
            Grid.link_node_face(n1, f1)
            Grid.link_node_face(n3, f1)
            Grid.link_node_face(n4, f1)
            Grid.link_node_face(n2, f2)
            Grid.link_node_face(n3, f2)
            Grid.link_node_face(n4, f2)

            Grid.link_edge_face(e1, f1)
            Grid.link_edge_face(e2, f1)
            Grid.link_edge_face(e5, f1)
            Grid.link_edge_face(e3, f2)
            Grid.link_edge_face(e4, f2)
            Grid.link_edge_face(e5, f2)

            zone.add_face(f1)
            zone.add_face(f2)
            zone.add_edge(e5)

            self.add_face(f1, gloid1)
            self.add_face(f2, gloid2)

            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------------

    def cut_edge(self, edge, p):
        """

        Atomic grid transformation. Splitting an edge.

        Parameters
        ----------
        edge: Edge
        p: Vect

        Returns
        -------

                    n3
                   *
                  /|\
                 / | \
                /  |  \
               /   |   \
              /    |    \
             /     |     \
            / f1_1 | f1_2 \
           /       |n5     \
        n2*--------*--------* n1
          \        |        /
           \ f2_2  | f2_1  /
            \      |      /
             \     |     /
              \    |    /
               \   |   /
                \  |  /
                 \ | /
                  \|/
                   *
                    n4
        """

        if len(edge.Faces) == 2:

            # det all data and del from grid
            f1 = edge.Faces[0]
            f2 = edge.Faces[1]
            data1 = f1.Data.copy()
            gloid1 = f1.GloId
            zone = f1.Zone
            data2 = f2.Data.copy()
            gloid2 = f2.GloId

            self.del_face(f1)
            self.del_face(f2)

            self.del_edge(edge)

            # new faces
            f1_1 = Face(list(data1.keys()), list(data1.values()))
            f1_2 = Face(list(data1.keys()), list(data1.values()))
            f2_1 = Face(list(data2.keys()), list(data2.values()))
            f2_2 = Face(list(data2.keys()), list(data2.values()))

            self.add_face(f1_1, gloid1)
            self.add_face(f2_1, gloid2)
            self.add_face(f1_2)
            self.add_face(f2_2)

            # Nodes
            n1 = edge.Nodes[0]
            n2 = edge.Nodes[1]
            n3 = [i for i in f1.Nodes if id(i) != id(n1) and id(i) != id(n2)][0]
            n4 = [i for i in f2.Nodes if id(i) != id(n1) and id(i) != id(n2)][0]
            n5 = Node(p)
            node = self.add_node(n5, True)
            if id(node) == id(n5):
                zone.add_node(n5)
            else:
                n5 = node

            # Edges
            f1.Edges.remove(edge)
            f2.Edges.remove(edge)
            e13 = [ed for ed in f1.Edges for node in ed.Nodes if id(node) == id(n1)][0]
            e23 = [ed for ed in f1.Edges for node in ed.Nodes if id(node) == id(n2)][0]
            e14 = [ed for ed in f2.Edges for node in ed.Nodes if id(node) == id(n1)][0]
            e24 = [ed for ed in f2.Edges for node in ed.Nodes if id(node) == id(n2)][0]
            e15 = Edge()
            e25 = Edge()
            e35 = Edge()
            e45 = Edge()
            self.add_edge(e15, edge.GloId)
            self.add_edge(e25)
            self.add_edge(e35)
            self.add_edge(e45)

            # Links
            Grid.link_node_edge(n1, e15)
            Grid.link_node_edge(n5, e15)
            Grid.link_node_edge(n2, e25)
            Grid.link_node_edge(n5, e25)
            Grid.link_node_edge(n3, e35)
            Grid.link_node_edge(n5, e35)
            Grid.link_node_edge(n4, e45)
            Grid.link_node_edge(n5, e45)

            Grid.link_node_face(n5, f1_1)
            Grid.link_node_face(n5, f1_2)
            Grid.link_node_face(n5, f2_1)
            Grid.link_node_face(n5, f2_2)
            Grid.link_node_face(n1, f1_2)
            Grid.link_node_face(n1, f2_1)
            Grid.link_node_face(n2, f1_1)
            Grid.link_node_face(n2, f2_2)
            Grid.link_node_face(n3, f1_1)
            Grid.link_node_face(n3, f1_2)
            Grid.link_node_face(n4, f2_1)
            Grid.link_node_face(n4, f2_2)

            Grid.link_edge_face(e15, f1_2)
            Grid.link_edge_face(e15, f2_1)
            Grid.link_edge_face(e25, f1_1)
            Grid.link_edge_face(e25, f2_2)
            Grid.link_edge_face(e35, f1_2)
            Grid.link_edge_face(e35, f1_1)
            Grid.link_edge_face(e45, f2_2)
            Grid.link_edge_face(e45, f2_1)
            Grid.link_edge_face(e13, f1_2)
            Grid.link_edge_face(e23, f1_1)
            Grid.link_edge_face(e14, f2_1)
            Grid.link_edge_face(e24, f2_2)

            zone.add_face(f1_1)
            zone.add_face(f1_2)
            zone.add_face(f2_1)
            zone.add_face(f2_2)
            zone.add_edge(e15)
            zone.add_edge(e25)
            zone.add_edge(e35)
            zone.add_edge(e45)

            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------------

    def collapse_edge(self, edge):
        """

        Atomic grid transformation. Tightening along the edge.

        Parameters
        ----------
        edge: Edge

        Returns
        -------

        """

        if len(edge.Faces) == 2:

            # det all data and del from grid
            f1 = edge.Faces[0]
            f2 = edge.Faces[1]

            self.del_face(f1)
            self.del_face(f2)

            self.del_edge(edge)

            # Nodes
            n1 = edge.Nodes[0]
            n2 = edge.Nodes[1]

            vect_coord_new_p = Vect.middle_point(n1.P, n2.P)

            self.RoundedCoordsBag.remove(n1.RoundedCoords)
            n1.P = vect_coord_new_p
            n1.RoundedCoords = n1.P.rounded_coords_tuple(10)
            self.RoundedCoordsBag.add(n1.RoundedCoords)

            # edges
            f1.Edges.remove(edge)
            e1 = [ed for ed in f1.Edges for node in ed.Nodes if id(node) == id(n2)][0]
            f2.Edges.remove(edge)
            e2 = [ed for ed in f2.Edges for node in ed.Nodes if id(node) == id(n2)][0]

            if e1.Faces:
                e3 = [ed for ed in f1.Edges for node in ed.Nodes if id(node) == id(n1)][0]
                Grid.link_edge_face(e3, e1.Faces[0])
            if e2.Faces:
                e4 = [ed for ed in f2.Edges for node in ed.Nodes if id(node) == id(n1)][0]
                Grid.link_edge_face(e4, e2.Faces[0])

            self.del_edge(e1)
            self.del_edge(e2)

            n2.replacing_node(n1)
            self.del_node(n2)

            self.Faces[-1].GloId = f1.GloId
            self.Faces[-2].GloId = f2.GloId
            self.Edges[-1].GloId = edge.GloId
            self.Edges[-2].GloId = e1.GloId
            self.Edges[-3].GloId = e2.GloId
            self.Nodes[-1].GloId = n2.GloId

            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------------

    def cut_single_edge(self, edge, p):
        """

        Atomic grid transformation. Splitting a single edge.

        Parameters
        ----------
        edge: Edge
        p: Vect

        Returns
        -------

                    n3
                   *
                  /|\
                 / | \
                /  |  \
               /   |   \
              /    |    \
             /     |     \
            / f1   |   f2 \
           /       |n4     \
        n1*--------*--------* n2

        """

        if len(edge.Faces) == 1:

            # det all data and del from grid
            face = edge.Faces[0]
            data = face.Data.copy()
            gloid = face.GloId
            zone = face.Zone

            self.del_face(face)
            self.del_edge(edge)

            # new faces
            f1 = Face(list(data.keys()), list(data.values()))
            f2 = Face(list(data.keys()), list(data.values()))

            self.add_face(f1, gloid)
            self.add_face(f2)

            # Nodes
            n1 = edge.Nodes[0]
            n2 = edge.Nodes[1]
            n3 = [i for i in face.Nodes if id(i) != id(n1) and id(i) != id(n2)][0]
            n4 = Node(p)
            node = self.add_node(n4, True)
            if id(node) == id(n4):
                zone.add_node(n4)
            else:
                n4 = node

            # Edges
            face.Edges.remove(edge)
            e13 = [ed for ed in face.Edges for node in ed.Nodes if id(node) == id(n1)][0]
            e23 = [ed for ed in face.Edges for node in ed.Nodes if id(node) == id(n2)][0]
            e14 = Edge()
            e24 = Edge()
            e34 = Edge()
            self.add_edge(e14, edge.GloId)
            self.add_edge(e24)
            self.add_edge(e34)

            # Links
            Grid.link_node_edge(n4, e14)
            Grid.link_node_edge(n4, e24)
            Grid.link_node_edge(n4, e34)
            Grid.link_node_edge(n1, e14)
            Grid.link_node_edge(n2, e24)
            Grid.link_node_edge(n3, e34)

            Grid.link_node_face(n3, f1)
            Grid.link_node_face(n3, f2)
            Grid.link_node_face(n4, f1)
            Grid.link_node_face(n4, f2)
            Grid.link_node_face(n1, f1)
            Grid.link_node_face(n2, f2)

            Grid.link_edge_face(e34, f1)
            Grid.link_edge_face(e34, f2)
            Grid.link_edge_face(e13, f1)
            Grid.link_edge_face(e23, f2)
            Grid.link_edge_face(e14, f1)
            Grid.link_edge_face(e24, f2)

            zone.add_face(f1)
            zone.add_face(f2)
            zone.add_edge(e14)
            zone.add_edge(e24)
            zone.add_edge(e34)

            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------------

    def cut_edge_with_two_nodes(self, edge, p1, p2):
        """

        Atomic grid transformation. Splitting an edge with two points.

        Parameters
        ----------
        edge: Edge
        p1: Vect
        p2: Vect

        Returns
        -------
                    n3
                   *
                  /|\
                 /|| \
                / | \ \
               /  |  \ \
              /   |   \ \
             /   /    |  \
            / f1|  f2 | f3\
           /    |      \   \
        n1*-----*------*----* n2
               n4       n5
        """

        if len(edge.Faces) == 1:

            # det all data and del from grid
            face = edge.Faces[0]
            data = face.Data.copy()
            gloid = face.GloId
            zone = face.Zone

            self.del_face(face)
            self.del_edge(edge)

            # new faces
            f1 = Face(list(data.keys()), list(data.values()))
            f2 = Face(list(data.keys()), list(data.values()))
            f3 = Face(list(data.keys()), list(data.values()))

            self.add_face(f1, gloid)
            self.add_face(f2)
            self.add_face(f3)

            # Nodes
            n1 = edge.Nodes[0]
            n2 = edge.Nodes[1]
            n3 = [i for i in face.Nodes if id(i) != id(n1) and id(i) != id(n2)][0]
            n4 = Node(p1)
            node = self.add_node(n4, True)
            if id(node) == id(n4):
                zone.add_node(n4)
            else:
                n4 = node
            n5 = Node(p2)
            node = self.add_node(n5, True)
            if id(node) == id(n5):
                zone.add_node(n5)
            else:
                n5 = node

            # Edges
            face.Edges.remove(edge)
            e13 = [ed for ed in face.Edges for node in ed.Nodes if id(node) == id(n1)][0]
            e23 = [ed for ed in face.Edges for node in ed.Nodes if id(node) == id(n2)][0]
            e14 = Edge()
            e25 = Edge()
            e45 = Edge()
            e34 = Edge()
            e35 = Edge()
            self.add_edge(e14, edge.GloId)
            self.add_edge(e25)
            self.add_edge(e45)
            self.add_edge(e34)
            self.add_edge(e35)

            # Links
            Grid.link_node_edge(n1, e14)
            Grid.link_node_edge(n2, e25)
            Grid.link_node_edge(n3, e34)
            Grid.link_node_edge(n3, e35)
            Grid.link_node_edge(n4, e14)
            Grid.link_node_edge(n4, e34)
            Grid.link_node_edge(n4, e45)
            Grid.link_node_edge(n5, e25)
            Grid.link_node_edge(n5, e35)
            Grid.link_node_edge(n5, e45)

            Grid.link_node_face(n1, f1)
            Grid.link_node_face(n2, f3)
            Grid.link_node_face(n3, f1)
            Grid.link_node_face(n3, f2)
            Grid.link_node_face(n3, f3)
            Grid.link_node_face(n4, f1)
            Grid.link_node_face(n4, f2)
            Grid.link_node_face(n5, f2)
            Grid.link_node_face(n5, f3)

            Grid.link_edge_face(e14, f1)
            Grid.link_edge_face(e13, f1)
            Grid.link_edge_face(e34, f1)
            Grid.link_edge_face(e34, f2)
            Grid.link_edge_face(e45, f2)
            Grid.link_edge_face(e35, f2)
            Grid.link_edge_face(e25, f3)
            Grid.link_edge_face(e35, f3)
            Grid.link_edge_face(e23, f3)

            zone.add_face(f1)
            zone.add_face(f2)
            zone.add_face(f3)
            zone.add_edge(e14)
            zone.add_edge(e25)
            zone.add_edge(e45)
            zone.add_edge(e34)
            zone.add_edge(e35)

            return True

        elif len(edge.Faces) == 2:

            # det all data and del from grid
            face1 = edge.Faces[0]
            face2 = edge.Faces[1]
            data1 = face1.Data.copy()
            gloid1 = face1.GloId
            zone = face1.Zone
            data2 = face2.Data.copy()
            gloid2 = face2.GloId

            self.del_face(face1)
            self.del_face(face2)
            self.del_edge(edge)

            # new faces
            f1 = Face(list(data1.keys()), list(data1.values()))
            f2 = Face(list(data1.keys()), list(data1.values()))
            f3 = Face(list(data1.keys()), list(data1.values()))
            f4 = Face(list(data2.keys()), list(data2.values()))
            f5 = Face(list(data2.keys()), list(data2.values()))
            f6 = Face(list(data2.keys()), list(data2.values()))

            self.add_face(f1, gloid1)
            self.add_face(f2, gloid2)
            self.add_face(f3)
            self.add_face(f4)
            self.add_face(f5)
            self.add_face(f6)

            # Nodes
            n1 = edge.Nodes[0]
            n2 = edge.Nodes[1]
            n3 = [i for i in face1.Nodes if id(i) != id(n1) and id(i) != id(n2)][0]
            n6 = [i for i in face2.Nodes if id(i) != id(n1) and id(i) != id(n2)][0]
            n4 = Node(p1)
            node = self.add_node(n4, True)
            if id(node) == id(n4):
                zone.add_node(n4)
            else:
                n4 = node
            n5 = Node(p2)
            node = self.add_node(n5, True)
            if id(node) == id(n5):
                zone.add_node(n5)
            else:
                n5 = node

            # Edges
            face1.Edges.remove(edge)
            face2.Edges.remove(edge)
            e13 = [ed for ed in face1.Edges for node in ed.Nodes if id(node) == id(n1)][0]
            e23 = [ed for ed in face1.Edges for node in ed.Nodes if id(node) == id(n2)][0]
            e16 = [ed for ed in face2.Edges for node in ed.Nodes if id(node) == id(n1)][0]
            e26 = [ed for ed in face2.Edges for node in ed.Nodes if id(node) == id(n2)][0]
            e14 = Edge()
            e25 = Edge()
            e45 = Edge()
            e34 = Edge()
            e35 = Edge()
            e46 = Edge()
            e56 = Edge()
            self.add_edge(e14, edge.GloId)
            self.add_edge(e25)
            self.add_edge(e45)
            self.add_edge(e34)
            self.add_edge(e35)
            self.add_edge(e46)
            self.add_edge(e56)

            # Links
            Grid.link_node_edge(n1, e14)
            Grid.link_node_edge(n2, e25)
            Grid.link_node_edge(n3, e34)
            Grid.link_node_edge(n3, e35)
            Grid.link_node_edge(n4, e14)
            Grid.link_node_edge(n4, e34)
            Grid.link_node_edge(n4, e45)
            Grid.link_node_edge(n4, e46)
            Grid.link_node_edge(n5, e25)
            Grid.link_node_edge(n5, e35)
            Grid.link_node_edge(n5, e45)
            Grid.link_node_edge(n5, e56)
            Grid.link_node_edge(n6, e46)
            Grid.link_node_edge(n6, e56)

            Grid.link_node_face(n1, f1)
            Grid.link_node_face(n1, f4)
            Grid.link_node_face(n2, f3)
            Grid.link_node_face(n2, f6)
            Grid.link_node_face(n3, f1)
            Grid.link_node_face(n3, f2)
            Grid.link_node_face(n3, f3)
            Grid.link_node_face(n4, f1)
            Grid.link_node_face(n4, f2)
            Grid.link_node_face(n4, f4)
            Grid.link_node_face(n4, f5)
            Grid.link_node_face(n5, f2)
            Grid.link_node_face(n5, f3)
            Grid.link_node_face(n5, f5)
            Grid.link_node_face(n5, f6)
            Grid.link_node_face(n6, f4)
            Grid.link_node_face(n6, f5)
            Grid.link_node_face(n6, f6)

            Grid.link_edge_face(e14, f1)
            Grid.link_edge_face(e13, f1)
            Grid.link_edge_face(e34, f1)
            Grid.link_edge_face(e34, f2)
            Grid.link_edge_face(e45, f2)
            Grid.link_edge_face(e35, f2)
            Grid.link_edge_face(e25, f3)
            Grid.link_edge_face(e35, f3)
            Grid.link_edge_face(e23, f3)
            Grid.link_edge_face(e16, f4)
            Grid.link_edge_face(e14, f4)
            Grid.link_edge_face(e46, f4)
            Grid.link_edge_face(e45, f5)
            Grid.link_edge_face(e46, f5)
            Grid.link_edge_face(e56, f5)
            Grid.link_edge_face(e26, f6)
            Grid.link_edge_face(e25, f6)
            Grid.link_edge_face(e56, f6)

            zone.add_face(f1)
            zone.add_face(f2)
            zone.add_face(f3)
            zone.add_face(f4)
            zone.add_face(f5)
            zone.add_face(f6)
            zone.add_edge(e14)
            zone.add_edge(e25)
            zone.add_edge(e45)
            zone.add_edge(e34)
            zone.add_edge(e35)
            zone.add_edge(e46)
            zone.add_edge(e56)

            return True

        else:
            return False

    # ----------------------------------------------------------------------------------------------

    def cut_two_edges_with_two_nodes(self, edge1, edge2, p1, p2):
        """

        Parameters
        ----------
        edge1: Edge
        edge2: Edge
        p1: Vect
        p2: Vect

        Returns
        -------

        """

        pass

# ==================================================================================================
