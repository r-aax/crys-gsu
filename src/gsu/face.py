"""
Face realization.
"""

import random

# ==================================================================================================


class Face:
    """
    Face of the grid.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, data):
        """
        Constructor face.
        :param data: Face data.
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
        :return: Temperature value.
        """

        return self.Data[0]

    # ----------------------------------------------------------------------------------------------

    def set_t(self, t):
        """
        Set temperature value to face data.
        :param t: Temperature.
        """

        self.Data[0] = t

    # ----------------------------------------------------------------------------------------------

    def get_hw(self):
        """
        Get water height value.
        :return: Water height.
        """

        return self.Data[1]

    # ----------------------------------------------------------------------------------------------

    def set_hw(self, hw):
        """
        Set water height value to face data.
        :param hw: Water height.
        """

        self.Data[1] = hw

    # ----------------------------------------------------------------------------------------------

    def get_hi(self):
        """
        Get ice height value.
        :return: Ice height.
        """

        return self.Data[2]

    # ----------------------------------------------------------------------------------------------

    def set_hi(self, hi):
        """
        Set ice height value to face data.
        :param hi: Ice height.
        """

        self.Data[2] = hi

    # ----------------------------------------------------------------------------------------------

    def get_beta(self):
        """
        Get Beta value.
        :return: Beta value.
        """

        return self.Data[4]

    # ----------------------------------------------------------------------------------------------

    def set_beta(self, beta):
        """
        Set Beta value.
        :param beta: Beta value.
        """

        self.Data[4] = beta

    # ----------------------------------------------------------------------------------------------

    def get_mimp2(self):
        """
        Get MImp2 value.
        :return: MImp2 value.
        """

        return self.Data[5]

    # ----------------------------------------------------------------------------------------------

    def set_mimp2(self, mimp2):
        """
        Set MImp2 value.
        :param mimp2: MImp2 value.
        """

        self.Data[5] = mimp2

    # ----------------------------------------------------------------------------------------------

    def get_vd2(self):
        """
        Get Vd2 value.
        :return: Vd2 value.
        """

        return self.Data[6]

    # ----------------------------------------------------------------------------------------------

    def set_vd2(self, vd2):
        """
        Set Vd2 value.
        :param vd2: Vd2 value.
        """

        self.Data[6] = vd2

    # ----------------------------------------------------------------------------------------------

    def get_glo_id_t_hw_hi_str(self):
        """
        Get string with global identifier, temperature, and water and ice heights.
        :return: String.
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


# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
