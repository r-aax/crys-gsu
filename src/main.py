"""
Main module.
"""

import gsu

if __name__ == '__main__':

    print('Testing crys-gsu main : begin.')

    test = 'wing_50'
    g = gsu.Grid()
    g.load('grids/' + test + '.dat')
    #g.print_info()

    '''
    print('Create initial version.')
    g.store('grids/' + test + '_origin.dat')
    '''

    '''
    print('Create mono version.')
    g.distribute_mono()
    g.print_info()
    g.store('grids/' + test + '_mono.dat')
    '''

    print('Create random version.')
    g.distribute_random()
    g.print_info()
    g.store('grids/' + test + '_random.dat')

    print('Create linear version.')
    g.distribute_linear()
    g.print_info()
    g.store('grids/' + test + '_linear.dat')

    '''
    print('Create uniform version.')
    g.distribute_uniform(lambda f: f.center_point()[0])
    g.print_info()
    g.store('grids/' + test + '_uniform.dat')
    '''

    print('Create rgrow version.')
    g.distribute_rgrow()
    g.print_info()
    g.store('grids/' + test + '_rgrow.dat')

    print('Create hierarchical version.')
    g.distribute_hierarchical([lambda f: f.center_point()[0],
                               lambda f: f.center_point()[1],
                               lambda f: f.center_point()[2]])
    g.print_info()
    g.store('grids/' + test + '_hierarchical.dat')

    print('Testing crys-gsu main : ok.')
