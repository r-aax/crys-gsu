"""
GSU main functions.
"""

import random
import math
import utils

# ==================================================================================================


def rounded_point(p, digits=10):
    """
    Get rounded point for constructing bags of points.
    :param p: point
    :param digits: accuracy
    :return: rounded point
    """

    return round(p[0], digits), round(p[1], digits), round(p[2], digits)

# --------------------------------------------------------------------------------------------------


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


class Node:
    """
    Node of the grid.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, p):
        """
        Constructor node.
        :param p: node point (tuple of coordinates)
        """

        # Global identifier (in grid numeration).
        self.GloId = -1

        self.P = p

        # Rounded coordinates for registration in set.
        self.RoundedCoords = rounded_point(p)

        # Links with edges and faces.
        self.Edges = []
        self.Faces = []

    # ----------------------------------------------------------------------------------------------

    def is_near(self, n):
        """
        Check if one node is near to another.
        :param n: another node
        :return: True - if nodes are near to each other, False - otherwise
        """

        return self.RoundedCoords == n.RoundedCoords

# ==================================================================================================


class Edge:
    """
    Edge of the grid.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        # Global identifier (in grid numeration).
        self.GloId = -1

        # Links to nodes and faces.
        self.Nodes = []
        self.Faces = []

    # ----------------------------------------------------------------------------------------------

    def is_border(self):
        """
        Check if edge is border edge.
        :return: True - if edge is border edge, False - otherwise
        """

        # Border edge has only one neighbour face.
        return len(self.Faces) == 1

    # ----------------------------------------------------------------------------------------------

    def is_cross(self):
        """
        Check if edge is cross-zones.
        :return: True - if edge is cross-zones, False - otherwise
        """

        # Cross-zone edge has two neighbour faces from different zones.

        faces_count = len(self.Faces)

        if faces_count == 1:
            return False
        elif faces_count == 2:
            return self.Faces[0].Zone != self.Faces[1].Zone
        else:
            raise Exception('Edge cannot has {0} neighbours faces.'.format(faces_count))

    # ----------------------------------------------------------------------------------------------

    def is_inner(self):
        """
        Check if edge is inner.
        :return: True - if edge is inner, False - otherwise
        """

        # Inner edge has two faces from one zone.

        faces_count = len(self.Faces)

        if faces_count == 1:
            return False
        elif faces_count == 2:
            return self.Faces[0].Zone == self.Faces[1].Zone
        else:
            raise Exception('Edge cannot has {0} neighbours faces.'.format(faces_count))

    # ----------------------------------------------------------------------------------------------

    def is_outer(self):
        """
        Check if edge is outer for its zone.
        :return: True - if edge is outer, False - otherwise
        """

        return not self.is_inner()

    # ----------------------------------------------------------------------------------------------

    def is_connect_zones(self, z0, z1):
        """
        Check if edge connect two given zones.
        :param z0: first zone
        :param z1: second zone
        :return: True - if edge connects two given zones, False - otherwise
        """

        if len(self.Faces) != 2:
            return False

        fz0, fz1 = self.Faces[0].Zone, self.Faces[1].Zone

        return ((z0 == fz0) and (z1 == fz1)) or ((z0 == fz1) and (z1 == fz0))

    # ----------------------------------------------------------------------------------------------

    def is_adjacent_with(self, e):
        """
        Check if edge adjacent with another edge.
        :param e: another edge
        :return: True - if edges are adjacent, False - otherwise
        """

        a0, a1 = self.Nodes[0], self.Nodes[1]
        b0, b1 = e.Nodes[0], e.Nodes[1]

        return (a0 == b0) or (a0 == b1) or (a1 == b0) or (a1 == b1)


    # ----------------------------------------------------------------------------------------------

    def get_ids(self):
        """
        Get 4 ids:
            first node id,
            second node id,
            min zone face id,
            max zone face id.
        :return:
        """

        first_node_id = self.Nodes[0].GloId
        second_node_id = self.Nodes[1].GloId
        ff = self.Faces[0]
        sf = self.Faces[1]
        fzi = ff.Zone.Id
        szi = sf.Zone.Id
        min_zone_id = min(fzi, szi)
        max_zone_id = max(fzi, szi)
        if fzi < szi:
            min_zone_face_id = ff.GloId
            max_zone_face_id = sf.GloId
        else:
            min_zone_face_id = sf.GloId
            max_zone_face_id = ff.GloId

        return (first_node_id, second_node_id, min_zone_id, max_zone_id, min_zone_face_id, max_zone_face_id)

# ==================================================================================================


class Face:
    """
    Face of the grid.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, data):
        """
        Constructor face.
        :param data: face data (tuple)
        """

        # Global identifier (in grid numeration).
        self.GloId = -1

        self.Data = data

        # Links with nodes and edges.
        self.Nodes = []
        self.Edges = []

        # Link to zone (each face belongs only to one single zone).
        self.Zone = None

    # ----------------------------------------------------------------------------------------------

    def get_t(self):
        """
        Get temperature value.
        :return: temperature value
        """

        return self.Data[0]

    # ----------------------------------------------------------------------------------------------

    def set_t(self, t):
        """
        Set temperature value to face data.
        :param t: temperature
        """

        self.Data[0] = t

    # ----------------------------------------------------------------------------------------------

    def get_hw(self):
        """
        Get water height value.
        :return: water height
        """

        return self.Data[1]

    # ----------------------------------------------------------------------------------------------

    def set_hw(self, hw):
        """
        Set water height value to face data.
        :param hw: water height
        """

        self.Data[1] = hw

    # ----------------------------------------------------------------------------------------------

    def get_hi(self):
        """
        Get ice height value.
        :return: ice height
        """

        return self.Data[2]

    # ----------------------------------------------------------------------------------------------

    def set_hi(self, hi):
        """
        Set ice height value to face data.
        :param hi: ice height
        """

        self.Data[2] = hi

    # ----------------------------------------------------------------------------------------------

    def get_beta(self):
        """
        Get Beta value.
        :return: Beta value
        """

        return self.Data[4]

    # ----------------------------------------------------------------------------------------------

    def set_beta(self, beta):
        """
        Set Beta value.
        :param beta: Beta value
        """

        self.Data[4] = beta

    # ----------------------------------------------------------------------------------------------

    def get_mimp2(self):
        """
        Get MImp2 value.
        :return: MImp2 value
        """

        return self.Data[5]

    # ----------------------------------------------------------------------------------------------

    def set_mimp2(self, mimp2):
        """
        Set MImp2 value.
        :param mimp2: MImp2 value
        """

        self.Data[5] = mimp2

    # ----------------------------------------------------------------------------------------------

    def get_vd2(self):
        """
        Get Vd2 value.
        :return: Vd2 value
        """

        return self.Data[6]

    # ----------------------------------------------------------------------------------------------

    def set_vd2(self, vd2):
        """
        Set Vd2 value.
        :param vd2: Vd2 value
        """

        self.Data[6] = vd2

    # ----------------------------------------------------------------------------------------------

    def get_glo_id_t_hw_hi_str(self):
        """
        Get string with global identifier, temperature, and water and ice heights.
        :return: string
        """

        # Random data for t, hw, hi.
        a = [self.get_t() + random.random(),
             self.get_hw() + random.random(),
             self.get_hi() + random.random()]
        i_list = [str(self.GloId)] + ['{0:.18e}'.format(ai) for ai in a]
        i_str = ' '.join(i_list)

        return i_str

    # ----------------------------------------------------------------------------------------------

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

    # ----------------------------------------------------------------------------------------------

    def unlink_from_zone(self):
        """
        Unlink face from zone.
        """

        if self.Zone is None:
            return

        # Face is linked.
        # Unlink it.
        self.Zone.Faces.remove(self)
        self.Zone = None

    # ----------------------------------------------------------------------------------------------

    def get_center(self):
        """
        Get center point.
        :return: center point
        """

        a, b, c, = self.Nodes[0].P, self.Nodes[1].P, self.Nodes[2].P

        return ((a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0)

# ==================================================================================================


class Zone:
    """
    Zone of the grid.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, name):
        """
        Constructor.
        :param name: name of zone
        """

        self.Id = -1

        self.Name = name

        # No nodes or faces in the zone yet.
        self.Nodes = []
        self.Edges = []
        self.Faces = []

        # Fixed zone flag.
        self.IsFixed = False

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

    def outer_edges_count(self):
        """
        Get count of outer edges.
        :return: Outer edges count.
        """

        return len([e for e in self.Edges if e.is_outer()])

    # ----------------------------------------------------------------------------------------------

    def faces_count(self):
        """
        Get count of faces.
        :return: faces count.
        """

        return len(self.Faces)

    # ----------------------------------------------------------------------------------------------

    def zone_quality_factor(self):
        """
        Zone quality factor.
        :return: Zone quality factor.
        """

        return math.sqrt(self.faces_count()) / self.outer_edges_count()

    # ----------------------------------------------------------------------------------------------

    def get_nodes_coord_slice_str(self, i):
        """
        Get string composed from i-th coord of all nodes.
        :param i: index of nodes coord
        :return: composed string
        """

        i_list = ['{0:.18e}'.format(node.P[i]) for node in self.Nodes]
        i_str = ' '.join(i_list)

        return i_str

    # ----------------------------------------------------------------------------------------------

    def get_faces_data_slice_str(self, i):
        """
        Get string composed from i-th elements of data of all faces.
        :param i: index of nodes data
        :return: composed string
        """

        i_list = ['{0:.18e}'.format(face.Data[i]) for face in self.Faces]
        i_str = ' '.join(i_list)

        return i_str

    # ----------------------------------------------------------------------------------------------

    def get_faces_global_ids_slice_str(self):
        """
        Get string composed from global identifiers of all faces.
        :return: composed string
        """

        i_list = [str(face.GloId) for face in self.Faces]
        i_str = ' '.join(i_list)

        return i_str

    # ----------------------------------------------------------------------------------------------

    def add_node(self, n):
        """
        Add node to zone.
        :param n: node
        :return: added node
        """

        # Just add node.
        self.Nodes.append(n)

        return n

    # ----------------------------------------------------------------------------------------------

    def add_edge(self, e):
        """
        Add edge to zone.
        :param e: edge
        :return: added edge
        """

        # Just add egde.
        self.Edges.append(e)

        return e

    # ----------------------------------------------------------------------------------------------

    def add_face(self, f):
        """
        Add face to zone (with link).
        :param f: face
        :return: added face
        """

        # If face is already link to some zome,
        # we have to unlink it first.
        if f.Zone is not None:
            f.unlink_from_zone()

        # Just add and set link to the zone.
        f.Zone = self
        self.Faces.append(f)

        return f

    # ----------------------------------------------------------------------------------------------

    def capture_nearest_face(self):
        """
        Capture face nearest to zone (breadth-first search).
        """

        for f in self.Faces:
            for e in f.Edges:
                if not e.is_border():
                    f2 = f.get_neighbour(e)
                    if f2.Zone is None:
                        self.add_face(f2)
                        return

# ==================================================================================================


class ZonesAdjacencyMatrix:
    """
    Matrix of zones adjacency.
    For example if there is 3 zones matrix should be the following:
      | i00 c01 c02 br0 |
      | c01 i11 c12 br1 |
      | c02 c12 i22 br2 |
      | br0 br1 br2   0 |
    where cxy - count of cross edges between x-th and y-th zones,
          ixx - count of inner edges for x-th zone,
          brx - count of border edges for x-th zone.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, es, zs):
        """
        Constructor.
        :param es: edges list
        :param zs: zones list
        """

        # Init size and zero matrix.
        zc = len(zs)
        self.ZonesCount = zc
        self.M = []
        # Do not copy arrays because we have to create arrays, not references.
        for i in range(zc + 1):
            self.M.append([0] * (zc + 1))

        # Calculate for each edge.
        for e in es:
            fc = len(e.Faces)
            if fc == 1:
                f0 = e.Faces[0]
                z0 = zs.index(f0.Zone)
                self.inc_border(z0)
            elif fc == 2:
                f0, f1 = e.Faces[0], e.Faces[1]
                z0, z1 = zs.index(f0.Zone), zs.index(f1.Zone)
                self.inc(z0, z1)
            else:
                raise Exception('Wrong edge faces count ({0}).'.format(fc))

    # ----------------------------------------------------------------------------------------------

    def inc(self, i, j):
        """
        Increment matrix element value.
        :param i: first zone index
        :param j: second zone index
        """

        self.M[i][j] += 1

        if i != j:
            self.M[j][i] += 1

    # ----------------------------------------------------------------------------------------------

    def inc_border(self, i):
        """
        Increment value of border edges count.
        :param i: zone number
        """

        self.inc(i, self.ZonesCount)

    # ----------------------------------------------------------------------------------------------

    def edges_statistics(self):
        """
        Get edges statistics.
        Statistics is a tuple with following elements:
          ec - full edges count
          bec - border edges count
          iec - inner edges count
          cec - cross edges  count
          becp - border edges count percent
          iecp - inner edges count percent
          cecp - cross edges count percent
        :return: tuple
        """

        ec = 0
        bec = 0
        iec = 0
        cec = 0

        # Count all lines without the last one.
        for i in range(self.ZonesCount):
            line = self.M[i]
            # Border edges count for this zone is in the last row.
            # Inner edges count for this zone is on the main diagonal of the matrix.
            # All values between these two cells are cross edges.
            bec += line[self.ZonesCount]
            iec += line[i]
            cec += sum(line[(i + 1):self.ZonesCount])

        # Total count and percents.
        ec = bec + iec + cec
        becp, iecp, cecp = 100.0 * bec / ec, 100.0 * iec / ec, 100.0 * cec / ec

        return ec, bec, iec, cec, becp, iecp, cecp

    # ----------------------------------------------------------------------------------------------

    def edges_statistics_string(self):
        """
        String of edges statistics:
        :return: string
        """

        ec, bec, iec, cec, becp, iecp, cecp = self.edges_statistics()

        return 'edges stats: ' \
               '{0} border ({1:.2f}%), ' \
               '{2} inner ({3:.2f}%), ' \
               '{4} cross ({5:.2f}%)'.format(bec, becp, iec, iecp, cec, cecp)

    # ----------------------------------------------------------------------------------------------

    def zone_cross_edges_array(self, zi):
        """
        Array with count of cross edges.
        :param zi: zone index
        :return: array with cross edges count.
        """

        line = self.M[zi]
        line2 = line[:]
        line2[zi] = 0
        line2[self.ZonesCount] = 0

        return line2

    # ----------------------------------------------------------------------------------------------

    def zone_max_cross_border_len(self, zi):
        """
        Maximum border length for given zone.
        :param zi: zone index
        :return: max zone border length
        """

        return max(self.zone_cross_edges_array(zi))

    # ----------------------------------------------------------------------------------------------

    def max_cross_border_len(self):
        """
        Max value of cross zones border lengths.
        :return: max border length
        """

        return max([self.zone_max_cross_border_len(zi) for zi in range(self.ZonesCount)])

    # ----------------------------------------------------------------------------------------------

    def zone_cross_edges_count(self, zi):
        """
        Get zone cross-edges count.
        :param zi: zone index
        :return: count of cross-edges for this zone
        """

        return sum(self.zone_cross_edges_array(zi))

    # ----------------------------------------------------------------------------------------------

    def zone_cross_borders_count(self, zi):
        """
        Get borders count for given zone.
        :param zi: zone index
        :return: borders count
        """

        return len([x for x in self.zone_cross_edges_array(zi) if x > 0])

# ==================================================================================================


class EdgesChain:
    """
    Chain of edges.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, g, zi, zj):
        """
        Constructor.
        :param g: grid
        :param zi: first zone index
        :param zj: second zone index
        """

        self.Edges = []
        self.Subchains = []

        for e in g.Edges:
            if e.is_connect_zones(g.Zones[zi], g.Zones[zj]):
                self.add_edge(e)

        self.split_edges_into_subchains()

    # ----------------------------------------------------------------------------------------------

    def add_edge(self, e):
        """
        Add edge into chain.
        :param e: edge
        """

        self.Edges.append(e)

    # ----------------------------------------------------------------------------------------------

    def edges_count(self):
        """
        Get edges count.
        :return: edges count
        """

        return len(self.Edges)

    # ----------------------------------------------------------------------------------------------

    def is_empty(self):
        """
        Check if edges chain is empty.
        :return: True - if is empty, False - otherwise.
        """

        return self.edges_count() == 0

    # ----------------------------------------------------------------------------------------------

    def split_edges_into_subchains(self):
        """
        Sort edges into right chain.
        """

        def try_to_add_subchain_to_subchains(subchains, n):
            for subchain in subchains:
                if n[-1].is_adjacent_with(subchain[0]):
                    r = list(range(len(n)))
                    r.reverse()
                    for i in r:
                        subchain.insert(0, n[i])
                    return True
                if n[0].is_adjacent_with(subchain[-1]):
                    r = range(len(n))
                    for i in r:
                        subchain.append(n[i])
                    return True
                return False

        # Empty subchains and
        # set each edge into itw own subchain.
        self.Subchains = [[e] for e in self.Edges]

        # Try to connect subchains while it is possible.
        while True:
            count_before = len(self.Subchains)
            rearranged_subchains = []
            for s in self.Subchains:
                if not try_to_add_subchain_to_subchains(rearranged_subchains, s):
                    rearranged_subchains.append(s)
            self.Subchains = rearranged_subchains
            count_after = len(self.Subchains)
            if count_after == count_before:
                break

        # Replace edges by ids.
        self.Subchains = [[e.get_ids() for e in subchain] for subchain in self.Subchains]

    # ----------------------------------------------------------------------------------------------

    def get_edges_ids(self):
        """
        Get edges identifiers.
        :return: list of identifiers
        """

        return [e.GloId for e in self.Edges]

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

        # Mode.
        self.Mode = ''

        # Variables.
        self.VariablesStr = ''
        self.Variables = []
        self.FaceVariablesCount = 0

        # Set empty sets of nodes, faces, zones.
        self.Nodes = []
        self.Edges = []
        self.Faces = []
        self.Zones = []

        # Rounded coordinates
        self.RoundedCoordsBag = set()

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

    def random_face(self):
        """
        Get random face.
        :return: random face
        """

        fc = self.faces_count()
        ind = random.randint(0, fc - 1)

        return self.Faces[ind]

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

    def add_edge(self, e):
        """
        Add edge to grid.
        :param e: edge
        :return: added edge
        """

        # Just add edge with global id correction.
        e.GloId = len(self.Edges)
        self.Edges.append(e)

        return e

    # ----------------------------------------------------------------------------------------------

    def add_face(self, f):
        """
        Add face.
        :param f: face
        :return: added face
        """

        # Just correct global id and add.
        f.GloId = len(self.Faces)
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

        node.Edges.append(edge)
        edge.Nodes.append(node)

    # ----------------------------------------------------------------------------------------------

    def link_node_face(node, face):
        """
        Link face with node.
        :param node: node
        :param face: face
        """

        node.Faces.append(face)
        face.Nodes.append(node)

    # ----------------------------------------------------------------------------------------------

    def link_edge_face(edge, face):
        """
        Link edge with face.
        :param edge: edge
        :param face: face
        """

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

        # Clear all objects of the grid.
        self.clear()

        # Open file and try to load it line by line.
        with open(filename, 'r') as f:
            line = f.readline()
            while line:

                if '# EXPORT_MODE=' in line:

                    # Head of grid.
                    mode_line = line
                    title_line = f.readline()
                    variables_line = f.readline()

                    # Parse all and check.
                    self.Mode = mode_line.split('=')[-1][:-1]
                    if 'TITLE=' not in title_line:
                        raise Exception('Wrong title line ({0}).'.format(title_line))
                    self.Name = title_line.split('=')[1][1:-2]
                    if 'VARIABLES=' not in variables_line:
                        raise Exception('Wrong variables line ({0}).'.format(variables_line))
                    self.VariablesStr = variables_line.split('=')[-1][:-1]
                    self.Variables = self.VariablesStr.replace('"', '').replace(',', '').split()
                    self.FaceVariablesCount = len(self.Variables) - 3

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
                                             '([4-{0}]=CELLCENTERED)'.format(len(self.Variables))
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
                        p = [c[0][i], c[1][i], c[2][i]]
                        node = Node(p)
                        node = self.add_node(node, is_merge_same_nodes)
                        zone.add_node(node)

                    # Read data for faces.
                    d = []
                    for i in range(self.FaceVariablesCount):
                        line = f.readline()
                        d.append([float(xi) for xi in line.split()])
                    for i in range(faces_to_read):
                        face = Face([d[j][i] for j in range(self.FaceVariablesCount)])
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

        # Correct grid if there is no MImp2, Vd2 fields present.
        # crys#141
        if not '"MImp2", "Vd2"' in self.VariablesStr:
            self.VariablesStr = self.VariablesStr.replace('"Beta", ', '"Beta", "MImp2", "Vd2", ')
            gi = self.Variables.index('Beta') + 1
            li = gi - 3
            # Сначала вставляем 'Vd2', а затем в ту же позицию 'MImp2'.
            self.Variables.insert(gi, 'Vd2')
            self.Variables.insert(gi, 'MImp2')
            self.FaceVariablesCount = self.FaceVariablesCount + 2
            for f in self.Faces:
                # Ставим в нужную позицию два нуля.
                f.Data.insert(li, 0.0)
                f.Data.insert(li, 0.0)

    # ----------------------------------------------------------------------------------------------

    def get_variable_index(self, variable_name):
        """
        Get index of variable by its name.
        :param variable_name: variable name
        :return: variable index
        """

        return self.Variables.index(variable_name)

    # ----------------------------------------------------------------------------------------------

    def convert_grid_stall_to_check_point(self):
        """
        Convert grid from STALL mode to CHCEK_POINT mode.
        """

        if self.Mode != 'STALL':
            raise Exception('convertion works only with STALL mode')

        self.Mode = 'CHECK_POINT'
        self.Variables = self.Variables[:-5]
        self.VariablesStr = '"X", "Y", "Z", ' \
                            '"T", "Hw", "Hi", "HTC", ' \
                            '"Beta", "MImp2", "Vd2", ' \
                            '"TauX", "TauY", "TauZ"'
        self.FaceVariablesCount -= 5

        for f in self.Faces:
            f.Data = f.Data[:-5]

    # ----------------------------------------------------------------------------------------------

    def clean_t_hw_hi(self):
        """
        Clean T, Hw, Hi fields.
        """

        for f in self.Faces:
            f.set_t(0.0)
            f.set_hw(0.0)
            f.set_hi(0.0)

    # ----------------------------------------------------------------------------------------------

    def store(self, filename):
        """
        Store grid to file.
        :param filename: file name
        """

        with open(filename, 'w', newline='\n') as f:

            # Store head.
            f.write('# EXPORT_MODE={0}\n'.format(self.Mode))
            f.write('TITLE="{0}"\n'.format(self.Name))
            f.write('VARIABLES={0}\n'.format(self.VariablesStr))

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
                f.write('VARLOCATION=([4-{0}]=CELLCENTERED)\n'.format(len(self.Variables)))

                # Write first 3 data items (X, Y, Z coordinates).
                for i in range(3):
                    f.write(zone.get_nodes_coord_slice_str(i) + ' \n')

                # Write rest faces data items.
                for i in range(self.FaceVariablesCount):
                    f.write(zone.get_faces_data_slice_str(i) + ' \n')

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
        :param filename_base: base of filename
        "param ts: timestamp string
        :param sf: suffixes of files
        """

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
                variables = self.Variables[:3] + ['GloId'] + self.Variables[3:]
                variables = ['"{0}"'.format(x) for x in variables]
                variables_str = ', '.join(variables)
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
                for i in range(self.FaceVariablesCount):
                    file.write(z.get_faces_data_slice_str(i) + ' \n')

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
                print('split zone {0} -> {1}, {2}.'.format(nm, nm + 'l', nm + 'r'))
                self.split_zone(self.Zones[0], extract_signs_funs)

        self.post_decompose()

    # ----------------------------------------------------------------------------------------------

    def each_zone_capture_nearest_face(self):
        """
        Each zone capture nearest face.
        """

        for z in self.Zones:
            z.capture_nearest_face()

    # ----------------------------------------------------------------------------------------------

    def align_cross_borders(self):
        """
        Print cross borders.
        """

        def subchain_hi_faces_dbl(s):
            dbl = 0
            for i in range(1, len(s)):
                if s[i][4] == s[i - 1][4]:
                    dbl = dbl + 1
            return dbl

        def subchain_lo_faces_dbl(s):
            dbl = 0
            for i in range(1, len(s)):
                if s[i][5] == s[i - 1][5]:
                    dbl = dbl + 1
            return dbl

        def subchain_hi_lo_faces_dbl(s):
            return (subchain_hi_faces_dbl(s), subchain_lo_faces_dbl(s))

        def need_align(dbl):
            return (dbl[0] > 0) and (dbl[1] > 0)

        def align_hi(s, g):
            for i in range(1, len(s)):
                if s[i][4] == s[i - 1][4]:
                    # Index i - 1. Delete next and pass face to lo zone.
                    s.pop(i)
                    ids = s[i - 1]
                    s[i - 1] = (0, 0, 0, 0, -i, ids[4])
                    face = g.Faces[ids[4]]
                    # old_zone = g.Zones[ids[2]]
                    new_zone = g.Zones[ids[3]]
                    # old_zone.Faces.remove(face)
                    new_zone.add_face(face)
                    return

        def align_lo(s, g):
            for i in range(1, len(s)):
                if s[i][5] == s[i - 1][5]:
                    # Index i - 1. Delete next and pass face to hi zone.
                    s.pop(i)
                    ids = s[i - 1]
                    s[i - 1] = (0, 0, 0, 0, -i, ids[5])
                    face = g.Faces[ids[5]]
                    # old_zone = g.Zones[ids[3]]
                    new_zone = g.Zones[ids[2]]
                    # old_zone.Faces.remove(face)
                    new_zone.add_face(face)
                    return

        zc = self.zones_count()

        for zi in range(zc):
            for zj in range(zi + 1, zc):
                ch = EdgesChain(self, zi, zj)

                if ch.is_empty():
                    continue

                # print('Edges between zones {0} and {1}'.format(zi, zj))
                # print('  {0}'.format(ch.get_edges_ids()))
                for sch in ch.Subchains:
                    # print('    sch:')
                    # for ids in sch:
                    #     print('      {0}'.format(ids))

                    # Correct subchain.
                    dbl = subchain_hi_lo_faces_dbl(sch)
                    # print('DBL before {0}'.format(dbl))
                    while need_align(dbl):
                        align_hi(sch, self)
                        if subchain_lo_faces_dbl(sch) > 0:
                            align_lo(sch, self)
                        dbl = subchain_hi_lo_faces_dbl(sch)
                    # print('DBL after {0}'.format(dbl))

        self.post_decompose()

    # ----------------------------------------------------------------------------------------------

    def decompose_pressure(self, count=32, new_name=None,
                           fz_names=[], is_align=False):
        """
        Create distribution based on pressure algorithm.
        :param count: zones count
        :param new_name: grid new name
        :param fz_names: list of fixed zones names
        :param is_align: need for cross borders align
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
            zone = Zone('pressure ' + str(i))
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

        if is_align:
            self.align_cross_borders()

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

        Warning!
        This function rewrite VARIABLES list of grid.

        :param filename: name of file
        """

        with open(filename, 'r') as f:
            line = f.readline()

            while line:

                if 'VARIABLES=' in line:

                    # Set new variables.

                    variables_xyz = '"X", "Y", "Z"'
                    variables_oth = line.split('=')[-1][9:-1]

                    # We know only basic modes.
                    if variables_oth == '"T", "Hw", "Hi"':
                        self.Mode = 'BASIC'
                    elif variables_oth == '"T", "Hw", "Hi", "HTC", ' \
                                          '"Beta", "MImp2", "Vd2", ' \
                                          '"TauX", "TauY", "TauZ"':
                        self.Mode = 'CHECK_POINT'
                    elif variables_oth == '"T", "Hw", "Hi", "HTC", ' \
                                          '"Beta", "MImp2", "Vd2", ' \
                                          '"TauX", "TauY", "TauZ", ' \
                                          '"Stall", "StallD", "StallVX", "StallVY", "StallVZ"':
                        self.Mode = 'STALL'
                    elif variables_oth == '"MassImpinged", "WaterFilmHeight", ' \
                                          '"CurrentIceGrowth", ' \
                                          '"TotalIceGrowth", "CurrentMassEvaporation", ' \
                                          '"SurfaceTemperature", ' \
                                          '"FilmVx", "FilmVy", "FilmVz", "FilmV", ' \
                                          '"IcesolConvectiveFlux", ' \
                                          '"EvapHeatFlux", "IceThickness", "HTC"':
                        self.Mode = 'CIAM'
                    else:
                        raise Exception('unknown export mode')

                    self.VariablesStr = variables_xyz + ', ' + variables_oth
                    self.Variables = self.VariablesStr.replace('"', '').replace(',', '').split()
                    self.FaceVariablesCount = len(self.Variables) - 3

                else:
                    ss = line.split()
                    glo_id = int(ss[0])
                    ss = ss[1:]
                    face = self.Faces[glo_id]
                    face.Data = [float(x) for x in ss]

                line = f.readline()

            f.close()

# ==================================================================================================
