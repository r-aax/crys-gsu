"""
Zone realization.
"""
import numpy as np

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

        return np.sqrt(self.faces_count()) / self.outer_edges_count()

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

    def coords_np_array(self):
        """Return coordiantes of nodes as [Nx3] numpy array.
        
        Return
        ------
        np array
        """
        x, y, z = list(), list(), list()
        for n in self.Nodes:
            x.append(n.P.X)
            y.append(n.P.Y)
            z.append(n.P.Z)
        return np.array([x, y, z], dtype=np.double).T

    def variables_np_array(self):
        """Return variables of nodes as [Fx12] numpy array.
        
        Return
        ------
        np array
        """
        gloids, ts, hws = list(), list(), list()
        his, htcs, betas = list(), list(), list()
        mimp2s = np.zeros((len(self.Faces)), dtype=np.float)
        vd2s =  np.zeros((len(self.Faces)), dtype=np.float)
        tauxs, tauys, tauzs = list(), list(), list()
        recoveryfactors =  np.zeros((len(self.Faces)), dtype=np.float)

        for f in self.Faces:
            gloids.append(f.GloId)
            ts.append(f['T'])
            hws.append(f['Hw'])
            his.append(f['Hi'])
            htcs.append(f['HTC'])
            betas.append(f['Beta'])
            tauxs.append(f['TauX'])
            tauys.append(f['TauY'])
            tauzs.append(f['TauZ'])
        return np.array([gloids, ts, hws, his, htcs, betas, mimp2s, vd2s, tauxs, tauys, tauzs, recoveryfactors], dtype=np.double).T

    def connectivity_list_faces_nodes(self):
        """Return connectivity list of faces and nodes.

        Return
        ------
        np array
        """
        connectivity_list = list()
        for f in self.Faces:
            n1, n2, n3 = f.Nodes[0].GloId, f.Nodes[1].GloId, f.Nodes[2].GloId
            assert not n1 == n2 == n3
            connectivity_list.append([n1, n2, n3])
        return np.array(connectivity_list, dtype=np.float)

    def get_variable(self, name):
        """Return numpy array of variable values.
        
        Parameters
        ----------
          name : str
            name of the variable

        Return
        ------
          numpy array
        """
        return np.array([f[name] for f in self.Faces])
    
    def get_real_face(self, e):
        """Get real face of en edge.
        
        Parameters
        ----------
          e : int
            local id of the edge
        
        Returns
        -------
          f : Face
            incident face that belongs to the zone
        """
        assert len(self.Edges[e].Faces) == 2

        f1 = self.Edges[e].Faces[0]
        f2 = self.Edges[e].Faces[1]

        if f1.Zone == self.Id:
            return f1
        elif f2.Zone == self.Id:
            return f2
        else:
            raise RuntimeError()

    def get_ghost_face(self, e):
        """Get real face of en edge.
        
        Parameters
        ----------
          e : int
            local id of the edge
        
        Returns
        -------
          f : Face
            incident face that belongs to the zone
        """
        assert len(self.Edges[e].Faces) == 2

        f1 = self.Edges[e].Faces[0]
        f2 = self.Edges[e].Faces[1]

        if f1.Zone == self.Id:
            return f2
        elif f2.Zone == self.Id:
            return f1
        else:
            raise RuntimeError()

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
