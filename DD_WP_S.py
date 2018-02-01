"""
Module to create a 2-Dimensional Gaussian Wavepacket at t=0
"""
import numpy as np

def Crystal(m, n):
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

    #Define the general square/rectangular lattice

    #Parameters
    a = 1.42
    b = 1.42

    X = np.zeros((m,n))
    Y = np.zeros((n,m))

    y_1 = np.arange(n)
    y_1 = y_1.reshape((n,1))
    Y[:,0:m] = - a*y_1
    Y = Y.reshape((n*m,1))
    Y_mn = Y.reshape((n,m))
    Y_T = Y_mn.T
    Y = Y_T.reshape((n*m,1))

    #print(Y)

    x_1 = np.arange(m)
    x_1 = x_1.reshape((m,1))
    X[:,] = - b*x_1
    X = X.reshape((n*m,1))

    #print(X)

    #well no wonder it is the same ... I have not done anything differently ...
    #It does not actually matter though ...
    #It will when they are not symmetric though!

    X_0 = 0.5*m*b
    Y_0 = -0.5*n*a

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
    Psi = np.zeros((n*m,1), dtype=complex)

    Psi = (np.exp(-0.5*(((X - X_0)/s)**2 + ((Y - Y_0)/s)**2))*np.exp((kx*X+ky*Y)*1j))/np.sqrt(np.pi*s)

    return np.array(Psi)

if __name__ == "__main__":
    points = Crystal(5, 10)
    #Psi(10*a, 0, 10/a, 10/a, points, 100, 100, plot = True)
