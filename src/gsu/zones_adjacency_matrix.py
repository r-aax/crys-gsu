"""
Zones adjacency matrix realization.
"""

# ==================================================================================================

class ZonesAdjacencyMatrix:
    """
    Matrix of zones adjacency.
    For example if there is 3 zones matrix should be the following::

      | i00 c01 c02 br0 |
      | c01 i11 c12 br1 |
      | c02 c12 i22 br2 |
      | br0 br1 br2   0 |

    where

    * cxy - count of cross edges between x-th and y-th zones,
    * ixx - count of inner edges for x-th zone,
    * brx - count of border edges for x-th zone.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, es, zs):
        """
        Constructor.

        :param es: Edges list.
        :param zs: Zones list.
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

        :param i: First zone index.
        :param j: Second zone index.
        """

        self.M[i][j] += 1

        if i != j:
            self.M[j][i] += 1

    # ----------------------------------------------------------------------------------------------

    def inc_border(self, i):
        """
        Increment value of border edges count.

        :param i: Zone number.
        """

        self.inc(i, self.ZonesCount)

    # ----------------------------------------------------------------------------------------------

    def edges_statistics(self):
        """
        Get edges statistics.
        Statistics is a tuple with following elements:

        * ec - full edges count
        * bec - border edges count
        * iec - inner edges count
        * cec - cross edges  count
        * becp - border edges count percent
        * iecp - inner edges count percent
        * cecp - cross edges count percent

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
        :return: String.
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

        :param zi: Zone index.
        :return:   Array with cross edges count.
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

        :param zi: Zone index.
        :return:   Max zone border length.
        """

        return max(self.zone_cross_edges_array(zi))

    # ----------------------------------------------------------------------------------------------

    def max_cross_border_len(self):
        """
        Max value of cross zones border lengths.

        :return: Max border length.
        """

        return max([self.zone_max_cross_border_len(zi) for zi in range(self.ZonesCount)])

    # ----------------------------------------------------------------------------------------------

    def zone_cross_edges_count(self, zi):
        """
        Get zone cross-edges count.

        :param zi: Zone index.
        :return:   Count of cross-edges for this zone.
        """

        return sum(self.zone_cross_edges_array(zi))

    # ----------------------------------------------------------------------------------------------

    def zone_cross_borders_count(self, zi):
        """
        Get borders count for given zone.

        :param zi: Zone index.
        :return:   Borders count.
        """

        return len([x for x in self.zone_cross_edges_array(zi) if x > 0])

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
