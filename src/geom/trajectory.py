"""
Trajectory realization.
Trajectory is a sequence of points.
"""

import math

# ==================================================================================================


class Trajectory:
    """
    Trajectory realization.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self, start):
        """
        Constructor.
        Init empty trajectory (with empty list of points).
        :param start: Start point.
        """

        self.Points = [start]

    # ----------------------------------------------------------------------------------------------

    def __str__(self):
        """
        String representation.
        :return: String.
        """

        return 'Trajectory [{0} - {1}, {2} points]'.format(self.Points[0],
                                                           self.Points[-1],
                                                           len(self.Points))

    # ----------------------------------------------------------------------------------------------

    def add_point(self, p):
        """
        Add point.
        :param p: Point.
        """

        self.Points.append(p)

    # ----------------------------------------------------------------------------------------------

    def points_count(self):
        """
        Get points count.
        :return: Points count.
        """

        return len(self.Points)

    # ----------------------------------------------------------------------------------------------

    def segments_count(self):
        """
        Get segments count.
        :return: Segments count.
        """

        return self.points_count() - 1

    # ----------------------------------------------------------------------------------------------

    def dump(self, f, zone_name):
        """
        Dump trajectory as separate zone to file.
        :param f: File.
        :param zone_name: Name of zone.
        """

        pc = self.points_count()

        f.write('ZONE T="{0}"\n'.format(zone_name))
        f.write('NODES={0}\n'.format(pc))
        f.write('ELEMENTS={0}\n'.format(pc - 1))
        f.write('DATAPACKING=BLOCK\n')
        # Paraview doesn't work with FELINESEG zone type, so imitate it as FETRIANGLE.
        f.write('ZONETYPE=FETRIANGLE\n')
        # Paraview demands VARLOCATION line.
        f.write('VARLOCATION=([4-3]=CELLCENTERED)\n')
        f.write(' '.join([str(p[0]) for p in self.Points]) + '\n')
        f.write(' '.join([str(p[1]) for p in self.Points]) + '\n')
        f.write(' '.join([str(p[2]) for p in self.Points]) + '\n')
        for ii in range(pc - 1):
            f.write('{0} {1} {1}\n'.format(ii + 1, ii + 2))

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
