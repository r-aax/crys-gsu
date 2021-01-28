"""
Main module.
"""

import gsu

if __name__ == '__main__':

    print('Testing crys-gsu main : begin.')

    test = 'air_inlet_2'
    g = gsu.Grid()
    g.load('grids/' + test + '.dat')
    g.print_info()

    """
    print('Create initial version.')
    g.store('grids/' + test + '_origin.dat')
    """

    """
    print('Create mono version.')
    g.distribute_mono()
    g.print_info()
    g.store('grids/' + test + '_mono.dat')
    """

    print('Create random version.')
    g.distribute_random()
    g.print_info()
    g.store('grids/' + test + '_random.dat')

    print('Testing crys-gsu main : ok.')
