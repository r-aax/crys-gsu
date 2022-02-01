"""
Node realization.
"""

from geom.vect import Vect

# ==================================================================================================


class Node:
    """
    Node of the grid.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, p: Vect):
        """
        Constructor node.

        :param p: Node point (Vector).
        """

        # Global identifier (in grid numeration).
        self.GloId = -1

        self.P = p

        # Rounded coordinates for registration in set.
        self.RoundedCoords = self.P.rounded_coords_tuple(10)

        # Links with edges and faces.
        self.Edges = []
        self.Faces = []

        self.border = False

    # ----------------------------------------------------------------------------------------------

    def is_near(self, n):
        """
        Check if one node is near to another.

        :param n: Another node
        :return:  True - if nodes are near to each other,
                  False - otherwise.
        """

        return self.RoundedCoords == n.RoundedCoords

    # ----------------------------------------------------------------------------------------------

    def __eq__(self, other):
        """

        Parameters
        ----------
        other: Another node

        Returns: True - if nodes are near to each other,
                 False - otherwise.
        -------

        """

        return self.is_near(other)

    # ----------------------------------------------------------------------------------------------

    def replacing_node(self, n):
        """

        Parameters
        ----------
        n: other Node

        Returns: the starting point will be replaced and consolidation by the point n
        -------

        """

        n.Edges += self.Edges
        n.Faces += self.Faces
        for ed in self.Edges:
            ed.Nodes.remove(self)
            ed.Nodes += [n]
        for face in self.Faces:
            face.Nodes.remove(self)
            face.Nodes += [n]

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
