#!/usr/bin/env python3

#======================
# Contains IO routines
#======================



#==========================================
def read_file():
#==========================================
    """
    Read swift output hdf5 file.
    """

    import h5py
    import settings as s

    f = h5py.File(s.srcfile)

    s.x = f[s.ptype]['Coordinates'][:,0]
    s.y = f[s.ptype]['Coordinates'][:,1]
    s.h = f[s.ptype]['SmoothingLength'][:]
    s.rho = f[s.ptype]['Density'][:]
    s.m = f[s.ptype]['Masses'][:]

    s.npart = s.x.shape[0]

    f.close()

    # TODO: remove returns
    return s.x, s.y, s.h, s.rho, s.m, s.npart



