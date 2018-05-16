import os

import numpy as np

from hwaves import hwf_density

tests_path = os.path.join(os.getcwd(),'tests')

nx = 30
ny = 30
nz = 30
dx = 0.2
dy = 0.2
dz = 0.2

def test_cartesian_density():
    ijk_xyz_PV = hwf_density.cartesian_density(nx,ny,nz,dx,dy,dz)
    fpath = os.path.join(tests_path,'cartesian_density_test.dat')
    hwf_density.write_cartesian(ijk_xyz_PV,fpath)
    assert(os.path.exists(fpath))

