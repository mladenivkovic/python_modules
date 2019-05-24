#!/usr/bin/env python3


#===============================================================
# Compute the volume of a particle in various ways and pray that
# you get the same results bae
#===============================================================


import numpy as np
import matplotlib.pyplot as plt


import meshless as ms



#---------------------------
# initialize variables
#---------------------------


# temp during rewriting
srcfile = './snapshot_for_hopkins.hdf5'    # swift output file
ptype = 'PartType0'                        # for which particle type to look for
pcoord = [0.5, 0.5]                        # coordinates of particle to work for





#========================
def main():
#========================
    

    x, y, h, rho, m, ids, npart = ms.read_file(srcfile, ptype)
    pind = ms.find_index(x, y, h, pcoord)
    nbors = ms.find_neighbours(pind, x, y, h)

    npart = x.shape[0]


    kernels = ['cubic_spline']
    
    for kernel in kernels:

        H = ms.get_H(h, kernel)

        # compute all psi_k(x_l) for all l, k
        psi_k_at_l = np.zeros((npart, npart), dtype=np.float128)

        for k in range(npart):
            for l in range(npart):
                psi_k_at_l[k,l] = ms.psi(x[l], y[l], x[k], y[k], H[l], kernel)

        neighbours = [[] for i in x]
        omega = np.zeros(npart, dtype=np.float128)

        for l in range(npart):

            # find and store all neighbours;
            neighbours[l] = ms.find_neighbours(l, x, y, H)

            # compute normalisation omega for all particles
            # needs psi_k_at_l to be computed already
            omega[l] =  np.sum(psi_k_at_l[:, l])
            # omega_k = sum_l W(x_k - x_l) = sum_l psi_l(x_k) as it is currently stored in memory



        V_i = ms.V(pind, m, rho)
        V_tot_array = np.sum(m/rho)

        comp_V = 1/omega[pind]
        V_tot_comp = np.sum(1.0/omega)
        print('{0:20} | {1:20} | {2:20} | {3:20}'.format(' ', '1/omega', 'm/rho', 'relative difference'))
        for i in range(20+3*20 + 3*3):
            print('-', end='')
        print()
        print('{0:20} | {1:20.10f} | {2:20.10f} | {3:20.3E}'.format("Single particle:", comp_V, V_i, (comp_V - V_i)/V_i))
        print('{0:20} | {1:20.10f} | {2:20.10f} | {3:20.3E}'.format("Total Volume:", V_tot_comp, V_tot_array, (V_tot_comp-V_tot_array)/V_tot_array))





if __name__ == '__main__':
    main()

