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

# ==================================================================================================


if __name__ == '__main__':
    pass

# ==================================================================================================
