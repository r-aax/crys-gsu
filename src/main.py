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
    '''
    g.decompose_mono(new_name=test + ' mono')
    g.print_info()
    g.store('grids/' + test + '_mono.dat')
    '''
    #
    '''
    g.decompose_random(new_name=test + ' random')
    g.print_info()
    g.store('grids/' + test + '_random.dat')
    '''
    #
    '''
    g.decompose_linear(new_name=test + ' linear')
    g.print_info()
    g.store('grids/' + test + '_linear.dat')
    '''
    #
    '''
    g.decompose_rgrow(new_name=test + ' rgrow')
    g.print_info()
    g.store('grids/' + test + '_rgrow.dat')
    '''
    #
    g.decompose_hierarchical([gsu.fun_face_cx(), gsu.fun_face_cy(), gsu.fun_face_cz()],
                             levels=6,
                             new_name=test + ' hierarchical')
    g.print_info(is_print_edges_statistics=True,
                 is_print_faces_distribution=True,
                 is_print_zones_adjacency_matrix=True)
    g.store('grids/' + test + '_hierarchical.dat')
    g.store_mpi('grids/' + test + '_mpi')

    print('Testing crys-gsu main : ok.')
