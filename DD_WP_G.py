"""
Module to create a 2-Dimensional Gaussian Wavepacket at t=0
"""
import numpy as np

def Crystal(m=20, n=10):
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

    a = 1.42#working in angstroms

    #might be able to ignore these ... as in we could write them as one ... 
    d_x = a*0.5
    d_y = 0.5*a*np.sqrt(3)

    #this shall be a list of all the possible x and y values ...
    X = np.zeros((m,n))
    Y = np.zeros((n,m))

    y_1 = np.arange(n)
    Y[:,] = - d_y*y_1.T
    Y = Y.reshape((n*m,1))

    x_1 = np.arange(m)
    dx = np.zeros((1,m))
    dx[0,0:m:2] = d_x*0.5
    dx[0,1:m:2] = -d_x*0.5
    X[:,] = (d_x + a)*x_1
    X[0:m:2,:] = X[0:m:2,:] + dx
    X[1:m:2,:] = X[1:m:2,:] - dx
    X = X.T.reshape((n*m,1))

    X_0 = 0.5*m*(d_x + a)
    Y_0 = -0.5*n*d_y

    return X, Y, X_0, Y_0

def Psi(s, kx, ky, m, n):
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

    X, Y, X_0, Y_0 = Crystal(m,n)

    # Calculate value of Gaussian at each atomic location
    Psi = np.zeros((n*m,1),dtype=complex)

    Psi = (np.exp(-0.5*(((X - X_0)/s)**2 + ((Y - Y_0)/s)**2))*np.exp((kx*X+ky*Y)*1j))/np.sqrt(np.pi*s)

    return Psi

if __name__ == "__main__":
    Crystal(m=10, n=10)
    #Psi(10, 1/10, 1/10, 100, 100)
