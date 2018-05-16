import os

import numpy as np
import matplotlib.pyplot as plt

from .hwf import spherical_harmonic, radial_wf

bohr_rad_A = 0.529177       #Angstrom
elec_mass_amu = 5.485799090 #amu

def plot_spherical_harmonic(theta,phi,l=0,m=0):
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
    return fig


def plot_radial_wf(r_A,Z=1,N_neut=0,n=1,l=0,showplot=False):

    Rnl = radial_wf(r_A,Z,N_neut,n,l)
    Rnlsqr = Rnl * Rnl
    Pnl = Rnlsqr * 4 * np.pi * r_A**2 

    fig = plt.figure()
    ax = fig.gca()
    ax.plot(Z*r_A,Rnl)
    ax.set_xlabel('Z*r (Angstrom)')
    ax.set_ylabel('Rnl')

    fig = plt.figure()
    ax = fig.gca()
    ax.plot(Z*r_A,Rnlsqr)
    ax.set_xlabel('Z*r (Angstrom)')
    ax.set_ylabel('Rnl**2')

    fig = plt.figure()
    ax = fig.gca()
    ax.plot(Z*r_A,Pnl)
    ax.set_xlabel('Z*r (Angstrom)')
    ax.set_ylabel('Pnl')

    if showplot:
        plt.show()
    return fig

def plot_isosurface():
    # TODO: take volumetric data, plot a surface
    pass

def plot_scatter():
    # TODO: take volumetric data, plot a 3d scatter
    pass




