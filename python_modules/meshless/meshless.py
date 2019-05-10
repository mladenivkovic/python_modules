#!/usr/bin/env python3

#===========================================
#
# A module containing common routines for 
# the meshless effective area visualisiation
# with 2d datasets
#
#===========================================


from meshlessio import read_file
from kernels import W, kernels
from particles import find_index, find_neighbours, V


import numpy as np


#=====================================
def Aij_Hopkins(pind, x, y, h, m, rho):
#=====================================
    """
    Compute A_ij as defined by Hopkins 2015
    i, j: particle indices for x, y, h, m arrays
    x, y, h, m: full data arrays as read in from hdf5 file
    """

    nbors = find_neighbours(pind, x, y, h)

    xj = x[nbors]
    yj = y[nbors]
    hj = h[nbors]

    #-------------------------------------------------------
    # Part 1: For particle at x_i (Our chosen particle)
    #-------------------------------------------------------

    # compute psi_j(x_i)
    psi_j = compute_psi(x[pind], y[pind], xj, yj, h[pind])

    # normalize psi_j
    omega_xi =  (np.sum(psi_j) + psi(x[pind], y[pind], x[pind], y[pind], h[pind]))
    psi_j /= omega_xi
    psi_j = np.float64(psi_j)

    # compute B_i
    B_i = get_matrix(x[pind], y[pind], xj, yj, psi_j)

    # compute psi_tilde_j(x_i)
    psi_tilde_j = np.empty((len(nbors), 2), dtype=np.float)
    for i, n in enumerate(nbors):
        dx = np.array([xj[i]-x[pind], yj[i]-y[pind]])
        psi_tilde_j[i] = np.dot(B_i, dx) * psi_j[i]



    #---------------------------------------------------------------------------
    # Part 2: values of psi/psi_tilde of particle i at neighbour positions x_j
    #---------------------------------------------------------------------------

    psi_i = np.zeros(len(nbors), dtype=np.float128)         # psi_i(xj)
    psi_tilde_i = np.empty((len(nbors), 2), dtype=np.float)  # psi_tilde_i(x_j)

    for i,n in enumerate(nbors):
        # first compute all psi(xj) from neighbour's neighbours to get weight omega
        nneigh = find_neighbours(n, x, y, h)
        xk = x[nneigh]
        yk = y[nneigh]
        for j, nn in enumerate(nneigh):
            psi_k = compute_psi(x[n], y[n], xk, yk, h[n])
            if nn == pind: # store psi_i, which is the psi for the particle whe chose at position xj; psi_i(xj)
                psi_i[i] = psi_k[j]
    
        omega_xj = (np.sum(psi_k) + psi(x[n], y[n], x[n], y[n], h[n]))

        psi_i[i]/= omega_xj


        # now compute B_j^{\alpha \beta}
        B_j = get_matrix(x[n], y[n], xk, yk, h[n])

        # get gradient
        dx = np.array([x[pind]-x[n], y[pind]-y[n]])
        psi_tilde_i[i] = np.dot(B_j, dx) * np.float64(psi_i[i])



    #-------------------------------
    # Part 3: Compute A_ij, x_ij
    #-------------------------------

    A_ij = np.empty((len(nbors),2), dtype = np.float)
    x_ij = np.empty((len(nbors),2), dtype = np.float)

    for i,n in enumerate(nbors):
        A_ij[i] = V(pind, m, rho)*psi_tilde_j[i] - V(n, m, rho)*psi_tilde_i[i]

        hfact = h[pind]/(h[pind]+h[n])
        x_ij[i] = np.array([x[pind]-hfact*(x[pind]-x[n]), y[pind]-hfact*(y[pind]-y[n])])



    return A_ij, x_ij





#=============================================================
def compute_psi(xi, yi, xj, yj, h, kernel='cubic_spline'):
#=============================================================
    """
    Compute all psi_j(x_i)
    xi, yi: floats
    xj, yj: arrays of neighbour's positions
    h: float
    """

    # psi_j(x_i)
    psi_j = np.zeros(xj.shape[0], dtype=np.float128)

    for i in range(xj.shape[0]):
        psi_j[i] = psi(xi, yi, xj[i], yj[i], h[i], kernel)

    return psi_j








#====================================================
def psi(x, y, xi, yi, hi, kernel='cubic_spline'):
#====================================================
    """
    UNNORMALIZED Volume fraction at position x of some particle
    with coordinates xi, yi
    ind: neighbour index in x/y/h array
    """
    q = np.float128(np.sqrt((x - xi)**2 + (y - yi)**2)/h)

    return W(q, hi, kernel)








#=============================================
def get_matrix(xi, yi, xj, yj, psi_j):
#=============================================
    """
    Get B_i ^{alpha beta}
    xi, yi: floats; Evaluate B at this position
    xj, yj: arrays; Neighbouring points
    psi_j:  array;  volume fraction of neighbours at position x_i; psi_j(x_i)
    """

    E00 = np.sum((xj-xi)**2 * psi_j)
    E01 = np.sum((xj-xi)*(yj-yi) * psi_j)
    E11 = np.sum((yj-yi)**2 * psi_j)
          
    E = np.matrix([[E00, E01], [E01, E11]])

    B = E.getI()
    return B


