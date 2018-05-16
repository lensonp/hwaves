from hwaves import hwf_density

nx = 10
ny = 10
nz = 10
dx = 0.02
dy = 0.02
dz = 0.02


ijk_xyz_PV = hwf_density.cartesian_density(nx,ny,nz,dx,dy,dz)

