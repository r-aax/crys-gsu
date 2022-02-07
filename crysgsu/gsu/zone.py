"""
Zone realization.
"""

# ==================================================================================================


class Zone:
    """
    Zone of the grid.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, name):
        """
        Constructor.

        :param name: Name of zone.
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

        :return: Nodes count.
        """

        return len(self.Nodes)

    # ----------------------------------------------------------------------------------------------

    def edges_count(self):
        """
        Get count of edges.

        :return: Edges count.
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

        :return: Faces count.
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

        :param i: Index of nodes coord.
        :return:  Composed string.
        """

        i_list = ['{0:.18e}'.format(node.P[i]) for node in self.Nodes]
        i_str = ' '.join(i_list)

        return i_str

    # ----------------------------------------------------------------------------------------------

    def get_faces_data_slice_str(self, e):
        """
        Get string composed from i-th elements of data of all faces.

        :param e: Data element.
        :return:  Composed string.
        """

        i_list = ['{0:.18e}'.format(face[e]) for face in self.Faces]
        i_str = ' '.join(i_list)

        return i_str

    # ----------------------------------------------------------------------------------------------

    def get_faces_global_ids_slice_str(self):
        """
        Get string composed from global identifiers of all faces.

        :return: Composed string.
        """

        i_list = [str(face.GloId) for face in self.Faces]
        i_str = ' '.join(i_list)

        return i_str

    # ----------------------------------------------------------------------------------------------

    def add_node(self, n):
        """
        Add node to zone.

        :param n: Node.
        :return:  Added node.
        """

        # Just add node.
        self.Nodes.append(n)

        return n

    # ----------------------------------------------------------------------------------------------

    def add_edge(self, e):
        """
        Add edge to zone.

        :param e: Edge.
        :return:  Added edge.
        """

        # Just add egde.
        self.Edges.append(e)

        return e

    # ----------------------------------------------------------------------------------------------

    def add_face(self, f):
        """
        Add face to zone (with link).

        :param f: Face.
        :return:  Added face.
        """

        # If face is already link to some zone,
        # we have to unlink it first.
        if f.Zone is not None:
            f.unlink_from_zone()

        # Just add and set link to the zone.
        f.Zone = self
        self.Faces.append(f)

        return f

# ==================================================================================================
