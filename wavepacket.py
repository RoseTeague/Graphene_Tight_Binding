"""
Module to create a 2-Dimensional Gaussian Wavepacket at t=0
"""
import numpy as np
import matplotlib.pyplot as plt


# Define initial parameter space
# 2D grid centered at (0,0) with points equally spaced by delta
xmin = -5*10**(-9)
xmax = 5*10**(-9)
ymin = -5*10**(-9)
ymax = 5*10**(-9)
delta = 10**(-11)
Nx = int((xmax-xmin)/delta)
Ny = int((ymax-ymin)/delta)
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)

# Define scale of the Wavepacket
kx = 5*10**(9)
ky = 5*10**(9)
s = 10**(-9) #Width of wavepacket
sigma = 0 #Phase

# Define Gaussian wavepacket
def Psi():
    """
    ===========================================================================
    Creation of a 2D Gaussian wavepacket centered at the origin
    ===========================================================================

    Parameters
    -----------

    Returns
    -----------
    Psi : ndarray, complex
        A Nx x Ny matrix, with the value of the wavefunction defined at each point
    """
    Psi = np.zeros((Nx,Ny),dtype = complex)
    for l in range(0,Nx):
        for m in range(0,Ny):
            Psi[l][m] = np.exp(-0.5*((x[l]/s)**2+(y[m]/s)**2))*np.exp((kx*x[l]+ky*y[m])*1j)
    return Psi
