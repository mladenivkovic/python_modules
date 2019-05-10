#!/usr/bin/env python3

#======================
# Contains IO routines
#======================



#==========================================
def read_file(srcfile, ptype):
#==========================================
    """
    Read swift output hdf5 file.
    """

    import h5py

    f = h5py.File(srcfile)

    x = f[ptype]['Coordinates'][:,0]
    y = f[ptype]['Coordinates'][:,1]
    h = f[ptype]['SmoothingLength'][:]
    rho = f[ptype]['Density'][:]
    m = f[ptype]['Masses'][:]

    npart = x.shape[0]

    f.close()

    return x, y, h, rho, m, npart



