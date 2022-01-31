"""
Face realization.
"""

import random
from geom.vect import Vect
from geom.triangle import Triangle

# ==================================================================================================


class Face:
    """
    Face of the grid.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, variables, values):
        """
        Constructor face.

        :param variables: Variables list.
        :param values:    Values list.
        """

        # Global identifier (in grid numeration).
        self.GloId = -1

        # Create face data as a dictionary.
        self.Data = dict(zip(variables, values))

        # Links with nodes and edges.
        self.Nodes = []
        self.Edges = []

        # Link to zone (each face belongs only to one single zone).
        self.Zone = None

    # ----------------------------------------------------------------------------------------------

    def __getitem__(self, item):
        """
        Get data field.

        :param item: Variable name.
        :return:     Value.
        """

        return self.Data.get(item, 0.0)

    # ----------------------------------------------------------------------------------------------

    def __setitem__(self, key, value):
        """
        Set data field.

        :param key:   Key value.
        :param value: Value.
        """

        self.Data[key] = value

    # ----------------------------------------------------------------------------------------------

    def get_glo_id_t_hw_hi_str(self):
        """
        Get string with global identifier, temperature, and water and ice heights.

        :return: String.
        """

        # Random data for t, hw, hi.
        a = [self['T'] + random.random(),
             self['Hw'] + random.random(),
             self['Hi'] + random.random()]
        i_list = [str(self.GloId)] + ['{0:.18e}'.format(ai) for ai in a]
        i_str = ' '.join(i_list)

        return i_str

    # ----------------------------------------------------------------------------------------------

    def get_neighbour(self, edge):
        """
        Get neighbour through edge.

        :param edge: Edge.
        :return:     Neighbour.
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

        :return: Center point.
        """

        a, b, c, = self.Nodes[0].P, self.Nodes[1].P, self.Nodes[2].P

        return ((a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0)

    # ----------------------------------------------------------------------------------------------

    def get_triangle(self):
        """
        Get triangle.

        :return: Triangle.
        """

        t = Triangle(Vect.from_iterable(self.Nodes[0].P),
                     Vect.from_iterable(self.Nodes[1].P),
                     Vect.from_iterable(self.Nodes[2].P))

        # Back reference for fixing parent object for triangle.
        t.BackRef = self

        return t

    # ----------------------------------------------------------------------------------------------

    def max_edge_from_face(self):
        """

        Returns: edge with max len
        -------

        """

        max_edge = self.Edges[0]
        max_len_edge = max_edge.len_edge()
        for edge in self.Edges[1:]:
            len_edge = edge.len_edge()
            if len_edge > max_len_edge:
                max_len_edge = len_edge
                max_edge = edge

        return max_edge

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
