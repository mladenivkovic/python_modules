#!/usr/bin/env python3

#===============================
# Particle related methods
#===============================




#===============================
def find_index():
#===============================
    """
    Find the index in the read-in arrays where
    the particle with coordinates of your choice is
    """

    import settings as s
    

    for i in range(s.npart):
        if abs(s.x[i]-s.pcoord[0]) < s.tolerance and abs(s.y[i] - s.pcoord[1]) < s.tolerance:
            s.pind = i
            break

    if s.verbose:
        print("got pind", s.pind, 'with x=', s.x[s.pind], 'y=', s.y[s.pind], 'h=', s.h[s.pind])
 
    #TODO: remove return
    return s.pind



