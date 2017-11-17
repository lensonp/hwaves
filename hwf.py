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

def spherical_harmonic(theta,phi,l,m):
    """Compute values of a spherical harmonic function.

    Parameters
    ----------
    theta : array
        Theta values (spherical coordinates)
        at which the spherical harmonic will be computed
    phi : array
        Phi values (spherical coordinates)
        at which the spherical harmonic will be computed
    l : int
        Angular momentum quantum number
    m : int
        Magnetic quantum number

    Returns
    -------
    Ylm : array
        Array of complex-valued spherical harmonics 
        computed at numpy.meshgrid(`theta`,`phi`) 
    """
    th, ph = np.meshgrid(theta,phi)
    Ylm = sph_harm(m,l,th,ph)
    return Ylm

def plot_spherical_harmonic(theta,phi,l,m):
    Ylm = spherical_harmonic(theta,phi,l,m)
    th, ph = np.meshgrid(theta,phi)
    sph_amp = np.abs(Ylm * np.conj(Ylm))
    x_amp = sph_amp * np.sin(ph) * np.cos(th)
    y_amp = sph_amp * np.sin(ph) * np.sin(th)
    z_amp = sph_amp * np.cos(ph)

    fig = plt.figure()
    #ax = Axes3D(fig)
    ax = fig.add_subplot(111,projection='3d')
    ax.plot_surface(x_amp, y_amp, z_amp,
    rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
    linewidth=0, antialiased=False, alpha=0.5)
    plt.show()
    #fig.savefig(savedir+'/Y{}{}.png'.format(l,m))

def radial_wf(r_A,n,l,Z_prot,N_neut):
    """Get wavefunction values wrt radial distance from the nucleus.

    Parameters
    ----------
    r_A : array of float
        array of radial points (in Angstroms)
        at which the wavefunction will be computed
    n : int
        principal quantum number
    l : int
        angular momentum quantum number
    Z_prot : int
        number of protons in the nucleus
    N_neut : int
        number of neutrons in the nucleus

    Returns
    ------- 
    Rnl : array
        array of complex wavefunction values
        at all of the input points `r_A`
    """
    
    # TODO: 
    # use nucmt to determine nuclear mass from Z+N
    # nucmt gives mass excess in energy units
    # M_nuc = Z + N + (mass_excess)/c^2
    #atmass = Get atomic mass of either most common or weighted-average isotope
    #nucmt = Table('AME2003')
    #mu = (reduced mass of nucleus-electron system, approx equal to m_e) 
    #a_mu = elec_mass_amu / mu * bohr_rad_A  
    a_mu = bohr_rad_A       # a_mu in units of Angstrom

    r_a_mu = r_A/a_mu       # radius in units of a_mu
    Zr_a_mu = Z*r_a_mu      # Z*r in units of a_mu

    # get generalized laguerre for n,l
    lagpoly = genlaguerre(n-l-1,2*l+1)
    lag_of_r = lagpoly(2*Zr_a_mu/n)

    # a_mu has units of Angstrom --> Rnl has units of Angstrom^(-3/2)
    Rnl = np.array(
    np.sqrt( (2*Z/float(n*a_mu))**3 * fact(n-l-1) / (2*n*fact(n+l)) )
    * np.exp(-1*Zr_a_mu/n)
    * (2*Zr_a_mu/n)**l 
    * lag_of_r)

    # Rnlsqr has units Angstrom^(-3)
    Rnlsqr = Rnl * Rnl

    # Pnl has units Angstrom^(-1)
    Pnl = Rnlsqr * 4 * np.pi * r_A**2  
    #Pnlsum = np.sum(Pnl) * max(r_A) / len(r_A)
    #print 'integral of Pnl: {}'.format(Pnlsum)

    return Rnl

def plot_radial_wf(r_A,n,l,Z_prot,N_neut):

    Rnl = radial_wf(r_A,n,l,Z_prot,N_neut)
    Rnlsqr = Rnl * Rnl
    Pnl = Rnlsqr * 4 * np.pi * r_A**2 

    fig = plt.figure()
    ax = fig.gca()
    ax.plot(Z*r_A,Rnl)
    ax.set_xlabel('Z*r (Angstrom)')
    ax.set_ylabel('Rnl')
    #fig.savefig(savedir+'/R{}{}.png'.format(n,l))
    #np.savetxt(savedir+'/R{}{}.csv'.format(n,l),np.array([Z*r_A,Rnl]),header='Z*r (Angstrom), Rnl') 

    fig = plt.figure()
    ax = fig.gca()
    ax.plot(Z*r_A,Rnlsqr)
    ax.set_xlabel('Z*r (Angstrom)')
    ax.set_ylabel('Rnl**2')
    #fig.savefig(savedir+'/R{}{}squared.png'.format(n,l))
    #np.savetxt(savedir+'/R{}{}squared.csv'.format(n,l),np.array([Z*r_A,Rnlsqr]),header='Z*r (Angstrom), Rnl*Rnl') 

    fig = plt.figure()
    ax = fig.gca()
    ax.plot(Z*r_A,Pnl)
    ax.set_xlabel('Z*r (Angstrom)')
    ax.set_ylabel('Pnl')
    #fig.savefig(savedir+'/P{}{}.png'.format(n,l))
    #np.savetxt(savedir+'/P{}{}.csv'.format(n,l),np.array([Z*r_A,Pnl]),header='Z*r (Angstrom), Rnl*Rnl*4*pi*r*r') 

    plt.show()

