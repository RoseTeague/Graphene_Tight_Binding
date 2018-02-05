"""
Module to create a 2-Dimensional Gaussian Wavepacket at t=0
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

#Parameters
a = 1.42*10**(-10)

def Crystals(a, Nx, Ny):
    """
    ===============================================================================
    Function to define the crystal sample for graphene
    ===============================================================================

    Parameters
    -----------
    a : float
        Lattice size parameter
    Nx : integer
        Number of atoms along the x-direction
    Ny : integer
        Numver of atoms along the y-direction

    Returns
    ----------
    coord : array
        list of atomic locations (x,y)
    """

    # Define the graphene Hexagonal structure
    deltax = 3*a
    deltay = np.sqrt(3)/2*a

    # Define initial parameter space
    # 2D grid centered at (0,0) with points equally spaced by delta
    xmin = -Nx/4*deltax
    xmax = Nx/4*deltax
    ymin = -Ny/2*deltay
    ymax = Ny/2*deltay

    # Compose the x-locations of atoms as the sum of 2 regularly spaced rows
    x1 = np.arange(xmin,xmax,deltax)
    x2 = x1-a

    # Final range of x- and y- atomic positions
    x = np.sort(x1.tolist()+x2.tolist())
    y = np.arange(ymin,ymax,deltay)

    # Rewrite atomic positions as a list of coordinates, with alternating rows
    # shifted by +/- 3a/2
    coord=[]
    sign = 1
    for j in y:
        for i in x:
            coord.append([i,j])
        sign = sign*(-1)
        x = x+sign*3*a/2

    # Convert coordinate list to array and plot atomic positions
    coord = np.array(coord)
    xpoint,ypoint = coord.T
    #plt.plot(xpoint,ypoint,color = 'gray',marker = '.',linestyle = 'None', markersize = '0.2')
    #plt.show()
    return x, y, coord

def Psi2(s, sigma, kx, ky, coord, m, n, plot=False):
    """
    ===========================================================================
    Creation of a 2D Gaussian wavepacket centered at the origin
    ===========================================================================

    Parameters
    -----------
    s : float
        Width of gaussian Wavepacket
    sigma : float
        Phase of guassian Wavepacket
    kx, ky : float
        wavenumbers along x and y directions
    coord : array
        List of atomic locations (x,y)
    m  : Number of atoms along x
    n  : Number of atoms along y


    Returns
    -----------
    Psi : ndarray, complex
    A m x n matrix, with the value of the wavefunction defined at each point
    """

    # Calculate value of Gaussian at each atomic location
    Psi = []
    for i in coord[2]:
        Psi.append(np.exp(-0.5*((i[0]/s)**2+(i[1]/s)**2))*np.exp((kx*i[0]+ky*i[1])*1j))

    # Normalise the wavefunction
    # A = scipy.integrate.simps(scipy.integrate.simps(Psi))
    # Psi = Psi/sqrt(A)
    Psi = Psi / np.sqrt(np.pi/s)
    Psi = np.reshape(Psi, (n*m,1))

    if plot:
        #Plot the absolute square of the wavefunction
        psi2 = np.array(np.abs(Psi)**2)
        psi2 = np.reshape(psi2,(n,m))
        plt.imshow(psi2)
        plt.show()
        #plt.savefig('Gaussian on Graphene')

    return np.array(Psi)

if __name__ == "__main__":
    a = 1.42*10**(-10)
    points = Crystal(a, 100, 100)
    psi = Psi(10*a, 0, 10/a, 10/a, points, 100, 100, plot = True)
