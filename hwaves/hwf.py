import os

import numpy as np
from scipy.special import sph_harm 
from scipy.special import genlaguerre 
from scipy.misc import factorial as fact
#from masstable import Table
#import periodictable as pt

bohr_rad_A = 0.529177       #Angstrom
elec_mass_amu = 5.485799090 #amu

def spherical_harmonic(theta,phi,l=0,m=0):
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

def radial_wf(r_A,Z,N_neut,n,l):
    """Get wavefunction values wrt radial distance from the nucleus.

    Parameters
    ----------
    r_A : array of float
        array of radial points (in Angstroms)
        at which the wavefunction will be computed
    Z : int
        number of protons in the nucleus
    N_neut : int
        number of neutrons in the nucleus
    n : int
        principal quantum number
    l : int
        angular momentum quantum number

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

    # Rnlsqr has units Angstrom^(-3): density per volume 
    Rnlsqr = Rnl * Rnl

    # Pnl has units Angstrom^(-1): spherical-shell integrated density per radius
    Pnl = Rnlsqr * 4 * np.pi * r_A**2  
    #Pnlsum = np.sum(Pnl) * max(r_A) / len(r_A)
    #print 'integral of Pnl: {}'.format(Pnlsum)

    return Rnl

def psi_xyz(x,y,z,Z,lag=None,n=1,l=0,m=0):
    # TODO: is this vectorized? 
    if not lag:
        lag = genlaguerre(n-l-1,2*l+1)
    r = np.sqrt(x**2+y**2+z**2)
    Zr_a0 = Z*r/bohr_rad_A
    lag_of_r = lag(2*Zr_a0/n)
    Rnl = np.sqrt( (2*Z/float(n*bohr_rad_A))**3 * fact(n-l-1) / (2*n*fact(n+l)) ) \
            * np.exp(-1*Zr_a0/n) * (2*Zr_a0/n)**l * lag_of_r
    th = np.arctan(y/x)
    ph = np.arctan(np.sqrt(x**2+y**2)/z)
    th_grid, ph_grid = np.meshgrid(th,ph)
    Ylm = sph_harm(m,l,th_grid,ph_grid) 
    return Rnl * Ylm

