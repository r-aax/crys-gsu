"""
Edge realization.
"""

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
        :return: True - if edge is border edge,
                 False - otherwise.
        """

        # Border edge has only one neighbour face.
        return len(self.Faces) == 1

    # ----------------------------------------------------------------------------------------------

    def is_cross(self):
        """
        Check if edge is cross-zones.
        :return: True - if edge is cross-zones,
                 False - otherwise.
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
        :return: True - if edge is inner,
                 False - otherwise.
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
        :return: True - if edge is outer,
                 False - otherwise.
        """

        return not self.is_inner()

    # ----------------------------------------------------------------------------------------------

    def is_connect_zones(self, z0, z1):
        """
        Check if edge connect two given zones.
        :param z0: First zone.
        :param z1: Second zone.
        :return:   True - if edge connects two given zones,
                   False - otherwise.
        """

        if len(self.Faces) != 2:
            return False

        fz0, fz1 = self.Faces[0].Zone, self.Faces[1].Zone

        return ((z0 == fz0) and (z1 == fz1)) or ((z0 == fz1) and (z1 == fz0))

    # ----------------------------------------------------------------------------------------------

    def is_adjacent_with(self, e):
        """
        Check if edge adjacent with another edge.
        :param e: Another edge.
        :return:  True - if edges are adjacent,
                  False - otherwise.
        """

        a0, a1 = self.Nodes[0], self.Nodes[1]
        b0, b1 = e.Nodes[0], e.Nodes[1]

        return (a0 == b0) or (a0 == b1) or (a1 == b0) or (a1 == b1)

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
