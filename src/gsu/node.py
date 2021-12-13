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

    # ----------------------------------------------------------------------------------------------

    def is_near(self, n):
        """
        Check if one node is near to another.

        :param n: Another node
        :return:  True - if nodes are near to each other,
                  False - otherwise.
        """

        return self.RoundedCoords == n.RoundedCoords

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
