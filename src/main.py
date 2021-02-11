"""
Main module.
"""

import gsu

if __name__ == '__main__':

    print('Testing crys-gsu main : begin.')

    test = 'bunny'
    g = gsu.Grid()
    g.load('grids/' + test + '.dat')
    g.print_info()

    # Distributions.
    g.distribute_mono(new_name=test + ' mono')
    g.print_info()
    g.store('grids/' + test + '_mono.dat')
    #
    g.distribute_random(new_name=test + ' random')
    g.print_info()
    g.store('grids/' + test + '_random.dat')
    #
    g.distribute_linear(new_name=test + ' linear')
    g.print_info()
    g.store('grids/' + test + '_linear.dat')
    #
    g.distribute_rgrow(new_name=test + ' rgrow')
    g.print_info()
    g.store('grids/' + test + '_rgrow.dat')
    #
    g.distribute_hierarchical([lambda f: f.center_point()[0],
                               lambda f: f.center_point()[1],
                               lambda f: f.center_point()[2]],
                              new_name=test + ' hierarchical')
    g.print_info()
    g.store('grids/' + test + '_hierarchical.dat')

    print('Testing crys-gsu main : ok.')
