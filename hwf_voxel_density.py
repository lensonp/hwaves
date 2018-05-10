import os

import numpy as np
from scipy import special
from scipy.special import sph_harm 
from scipy.special import genlaguerre 
from scipy.misc import factorial as fact
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
from masstable import Table
import periodictable as pt
from mayavi import mlab

bohr_rad_A = 0.529177       #Angstrom
elec_mass_amu = 5.485799090 #amu

def psi_xyz(x,y,z,Z,lag=None,n=1,l=0):
    if not lag:
        lag = genlaguerre(n-l-1,2*l+1)
    r = np.sqrt(x**2+y**2+z**2)
    Zr_a0 = Z*r/bohr_rad_A
    lag_of_r = lag(2*Zr_a0/n)
    Rnl = np.sqrt( (2*Z/float(n*bohr_rad_A))**3 * fact(n-l-1) / (2*n*fact(n+l)) ) * np.exp(-1*Zr_a0/n) * (2*Zr_a0/n)**l * lag_of_r
    th = np.arctan(y/x)
    ph = np.arctan(np.sqrt(x**2+y**2)/z)
    th_grid, ph_grid = np.meshgrid(th,ph)
    Ylm = sph_harm(m,l,th_grid,ph_grid) 
    return Rnl * Ylm

def plot_3d(n,l,m,Z=1,savedir='.'):
    th_step_deg = 5
    th_max_deg = 360
    ph_step_deg = 5
    ph_max_deg = 180

    #r_max_A = 2 
    #r_step_A = 0.4
    #r_step_A = 0.2
    #r_step_A = 0.1
    #r_step_A = 0.08
    #r_step_A = 0.06
    #r_step_A = 0.04

    r_max_A = 4 
    #r_step_A = 0.8
    #r_step_A = 0.4
    #r_step_A = 0.2
    #r_step_A = 0.16
    #r_step_A = 0.12
    r_step_A = 0.08

    # these arrays define the voxel centers
    th_deg = np.arange(float(th_step_deg)/2,th_max_deg,th_step_deg)
    ph_deg = np.arange(float(ph_step_deg)/2,ph_max_deg,ph_step_deg)
    r_A = np.arange(float(r_step_A)/2,r_max_A,r_step_A)

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
        V_r_th_ph[ir,:,:] = float(4)/3*np.pi*((Zri_center+Zr_step_A/2)**3 - (Zri_center-Zr_step_A/2)**3) * np.ones(shape=(len(th_deg),len(ph_deg)),dtype=float)

    P_spherical = np.array(np.abs(psi*np.conj(psi)),dtype=float)
    PV_spherical = np.array(P_spherical*V_r_th_ph,dtype=float)

    r_th_ph_PV = []
    for ir,Zri in zip(range(len(Zr_A)),Zr_A):
        for ith,thi in zip(range(len(th_deg)),th_deg):
            for iph,phi in zip(range(len(ph_deg)),ph_deg):
                r_th_ph_PV.append([Zri,thi,phi,PV_spherical[ir,ith,iph]]) 
    r_th_ph_PV = np.array(r_th_ph_PV,dtype=float)

    #x_A = np.arange(-1*r_max_A+float(r_step_A)/2, r_max_A, r_step_A)
    #y_A = np.arange(-1*r_max_A+float(r_step_A)/2, r_max_A, r_step_A)
    #z_A = np.arange(-1*r_max_A+float(r_step_A)/2, r_max_A, r_step_A)
    x_A = np.hstack([-1*r_A[::-1],r_A])
    y_A = np.hstack([-1*r_A[::-1],r_A])
    z_A = np.hstack([-1*r_A[::-1],r_A])
    r_idx = np.array(range(len(r_A)),dtype=int)
    x_idx = np.hstack([-1*r_idx[::-1]-1,r_idx+1])
    y_idx = np.hstack([-1*r_idx[::-1]-1,r_idx+1])
    z_idx = np.hstack([-1*r_idx[::-1]-1,r_idx+1])
    #import pdb; pdb.set_trace()
    xabs_A = np.abs(x_A)
    yabs_A = np.abs(y_A)
    zabs_A = np.abs(z_A)
    Zx_A = Z*x_A
    Zy_A = Z*y_A
    Zz_A = Z*z_A

    print 'computing in cartesian coords...'
    psi_cart = np.array([psi_xyz(xi,yabs_A,zabs_A,Z,lagpoly) for xi in xabs_A])
    P_cart = np.array(np.abs(psi_cart*np.conj(psi_cart)),dtype=float)
    PV_cart = np.array(P_cart*r_step_A**3,dtype=float)
    x_y_z_PV = []
    ijk_xyz_PV = []
    for ix,Zxi in zip(range(len(x_idx)),Zx_A):
        for iy,Zyi in zip(range(len(y_idx)),Zy_A):
            for iz,Zzi in zip(range(len(z_idx)),Zz_A):
    #            x_y_z_PV.append([Zxi,Zyi,Zzi,PV_cart[ix,iy,iz]])
                i_ix = x_idx[ix]
                j_iy = y_idx[iy]
                k_iz = z_idx[iz]
                ijk_xyz_PV.append([i_ix,j_iy,k_iz,Zxi,Zyi,Zzi,PV_cart[ix,iy,iz]])
    #x_y_z_PV = np.array(x_y_z_PV,dtype=float) 
    ijk_xyz_PV = np.array(ijk_xyz_PV,dtype=float) 

    print 'writing files...'    
    fspher = file(savedir+'/P_r_theta_phi.dat','a')
    np.savetxt(fspher,r_th_ph_PV,header='Z*r (Angstrom), theta (deg), phi (deg), voxel density')
    #fcart = file(savedir+'/P_x_y_z.dat','a')
    #np.savetxt(fcart,x_y_z_PV,header='Z*x (Angstrom), Z*y (Angstrom), Z*z (Angstrom), voxel density')
    fijk = file(savedir+'/ijk_xyz_P_boxsize{}A_{}points.dat'.format(2*r_max_A,np.shape(ijk_xyz_PV)[0]),'w')
    np.savetxt(fijk,ijk_xyz_PV,header='i_x, j_y, k_z, Z*x (Angstrom), Z*y (Angstrom), Z*z (Angstrom), voxel density',
        fmt=['%i','%i','%i','%.4f','%.4f','%.4f','%.8e'])

if __name__ == '__main__':
    Z = 1
    N = 0
    outdir = 'plots_Z{}'.format(Z)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    #nlm_run = [(1,0,0)]
    nlm_run = [(2,0,0)]
    #    (2,0,0),(2,1,-1),(2,1,0),(2,1,1),
    #    (3,0,0),(3,1,-1),(3,1,0),(3,1,1),
    #    (3,2,-2),(3,2,-1),(3,2,0),(3,2,1),(3,2,2)]

    for n,l,m in nlm_run:
        print '---------------------'
        print 'Z = {}, n = {}, l = {}, m = {}'.format(Z,n,l,m)
        nlmdir = outdir+'/n{}l{}m{}'.format(n,l,m)
        if not os.path.exists(nlmdir):
            os.mkdir(nlmdir)
        print 'generating voxel probability map for n={}, l={}, m={}...'.format(n,l,m)
        plot_3d(n,l,m,Z,nlmdir)
        

