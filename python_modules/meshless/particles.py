#!/usr/bin/env python3

#===============================
# Particle related methods
#===============================




#===============================
def find_index(x, y, h, pcoord):
#===============================
    """
    Find the index in the read-in arrays where
    the particle with coordinates of your choice is
    """

    tolerance = 1e-6

    for i in range(x.shape[0]):
        if abs(x[i]-pcoord[0]) < tolerance and abs(y[i] - pcoord[1]) < tolerance:
            pind = i
            break

    print("got pind", pind, 'with x=', x[pind], 'y=', y[pind], 'h=', h[pind])
 

    return pind





#============================================
def find_neighbours(ind, x, y, h, fact=2):
#============================================
    """
    Find indices of all neighbours within fact*h (where kernel != 0)
    """


    x0 = x[ind]
    y0 = y[ind]
    fhsq = h[ind]*h[ind]*fact*fact
    neigh = []

    for i in range(x.shape[0]):
        if i==ind:
            continue

        dist = (x[i]-x0)**2 + (y[i]-y0)**2
        if dist < fhsq:
            neigh.append(i)

    return neigh



#===================
def V(ind, m, rho):
#===================
    """
    Volume estimate for particle with index ind
    """

    return m[ind]/rho[ind]



