import os

import numpy as np
from scipy import special
from scipy.special import sph_harm 
from scipy.special import genlaguerre 
from scipy.misc import factorial as fact

from hwf import psi_xyz, bohr_rad_A, elec_mass_amu

def cartesian_density(n,l,m,nx,ny,nz,dx,dy,dz,Z=1,n=1,l=0):
    x = range(nx)*dx - (nx-1)*dx/2
    y = range(ny)*dy - (ny-1)*dy/2
    z = range(nz)*dz - (nz-1)*dz/2
    vol = dx * dy * dz
    xabs = np.abs(x)
    yabs = np.abs(y)
    zabs = np.abs(z)
    Zx = Z*x
    Zy = Z*y
    Zz = Z*z
    # generalized laguerre poly for n,l:
    lagpoly = genlaguerre(n-l-1,2*l+1)
    # wavefunction value array
    psi_cart = np.array([psi_xyz(xi,yabs,zabs,Z,lagpoly) for xi in xabs])
    # evaluate at voxel center
    P_cart = np.array(np.abs(psi_cart*np.conj(psi_cart)),dtype=float)
    # multiply by voxel volume
    PV_cart = np.array(P_cart*vol,dtype=float)
    # TODO: a more sophisticated voxel integration? 
    ijk_xyz_PV = []
    for ix,Zxi in enumerate(Zx_A):
        for iy,Zyi in enumerate(Zy_A):
            for iz,Zzi in enumerate(Zz_A):
                ijk_xyz_PV.append([x_idx[ix],y_idx[iy],z_idx[iz],Zxi,Zyi,Zzi,PV_cart[ix,iy,iz]])
    return ijk_xyz_PV

def spherical_density(n,l,m,r_max_A=3,nr=100,ntheta=72,nphi=36,Z=1):

    dtheta = 360./ntheta
    dphi = 180./nphi
    r_step_A = float(r_max_A)/nr 

    # these arrays define the voxel centers
    th_deg = np.arange(float(dtheta)/2,360,dtheta)
    print(len(th_deg))
    ph_deg = np.arange(float(dphi)/2,180,dphi)
    print(len(ph_deg))
    r_A = np.arange(float(r_step_A)/2,r_max_A,r_step_A)
    print(len(r_A))

    # various conversions etc
    th_rad = th_deg*np.pi/180 
    ph_rad = ph_deg*np.pi/180
    th_grid, ph_grid = np.meshgrid(th_rad,ph_rad)
    r_a0 = r_A/bohr_rad_A   
    Zr_A = Z*r_A
    Zr_a0 = Z*r_a0          
    Zr_step_A = Z*r_step_A
    
    # generalized laguerre poly for n,l:
    lagpoly = genlaguerre(n-l-1,2*l+1)
    # evaluate this at r
    lag_of_r = lagpoly(2*Zr_a0/n)

    # bohr_rad_A has units of Angstrom --> Rnl has units of Angstrom^(-3/2)
    Rnl = np.array(
    np.sqrt( (2*Z/float(n*bohr_rad_A))**3 * fact(n-l-1) / (2*n*fact(n+l)) )
    * np.exp(-1*Zr_a0/n)
    * (2*Zr_a0/n)**l 
    * lag_of_r
    )

    # spherical harmonic for m,l:
    Ylm = sph_harm(m,l,th_grid,ph_grid)
   
    # total wavefunction: psi(r,theta,psi) = Rnl(r)*Ylm(theta,psi)
    psi = np.zeros(shape=(len(r_A),len(th_deg),len(ph_deg)),dtype=complex) 
    V_r_th_ph = np.zeros(shape=(len(r_A),len(th_deg),len(ph_deg)),dtype=float) 
    print 'computing in spherical coords...'
    for ir in range(len(r_A)):
        psi[ir,:,:] = Ylm.T * Rnl[ir]
        Zri_center = Zr_A[ir]
        V_r_th_ph[ir,:,:] = float(4)/3*np.pi*((Zri_center+Zr_step_A/2)**3 - (Zri_center-Zr_step_A/2)**3) \
                            * np.ones(shape=(len(th_deg),len(ph_deg)),dtype=float)

    P_spherical = np.array(np.abs(psi*np.conj(psi)),dtype=float)
    PV_spherical = np.array(P_spherical*V_r_th_ph,dtype=float)

    r_th_ph_PV = []
    for ir,Zri in zip(range(len(Zr_A)),Zr_A):
        for ith,thi in zip(range(len(th_deg)),th_deg):
            for iph,phi in zip(range(len(ph_deg)),ph_deg):
                r_th_ph_PV.append([Zri,thi,phi,PV_spherical[ir,ith,iph]]) 
    return r_th_ph_PV


