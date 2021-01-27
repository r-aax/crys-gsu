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
    g.store('grids/' + test + '_new.dat')
    print('Testing crys-gsu main : ok.')
