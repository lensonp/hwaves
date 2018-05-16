import os

import numpy as np
import matplotlib.pyplot as plt

from hwf import spherical_harmonic, radial_wf

bohr_rad_A = 0.529177       #Angstrom
elec_mass_amu = 5.485799090 #amu

def plot_spherical_harmonic(theta,phi,l,m):
    Ylm = spherical_harmonic(theta,phi,l,m)
    th, ph = np.meshgrid(theta,phi)
    sph_amp = np.abs(Ylm * np.conj(Ylm))
    x_amp = sph_amp * np.sin(ph) * np.cos(th)
    y_amp = sph_amp * np.sin(ph) * np.sin(th)
    z_amp = sph_amp * np.cos(ph)

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot_surface(x_amp, y_amp, z_amp,
    rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
    linewidth=0, antialiased=False, alpha=0.5)
    plt.show()

def plot_radial_wf(r_A,n,l,Z_prot,N_neut):

    Rnl = radial_wf(r_A,Z_prot,N_neut,n,l)
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

