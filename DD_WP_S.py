"""
Module to create a 2-Dimensional Gaussian Wavepacket at t=0
"""
import numpy as np
import matplotlib.pyplot as plt

def Crystal(m, n):
    """
    ============================================================================
    Function to create a square/rectangular crystal
    ============================================================================

    Inputs
    -----------
    m : integer
        Number of atoms along the x-direction

    n : integer
        Numver of atoms along the y-direction

    Parameters
    -----------
    a : float
        Lattice size parameter of graphene

    Returns
    ----------
    X,
    Y,
    X_0,
    Y_0

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"

    #Lattice parameters in the y and x direction, respectively.
    a = 1.42
    b = 1.42

    #Creating initial arrays for populating the x and y values of each carbon atom.
    X = np.zeros((m,n))
    Y = np.zeros((n,m))

    #Creating a list from 0 to n for each atom in a column. This is reshaped
    #such that further operations can be performed.
    y_1 = np.arange(n)
    y_1 = y_1.reshape((n,1))

    #Each column has regularly spaced atoms.
    Y[:,0:m] = - a*y_1

    #Reshaping for further calculations
    Y_T = Y.T
    Y = Y_T.reshape((n*m,1))

    #Creating a list from 0 to n for each atom in a column. This is reshaped
    #such that further operations can be performed.
    x_1 = np.arange(m)
    x_1 = x_1.reshape((m,1))
    X[:,] = - b*x_1

    #Reshaping for further calculations
    X = X.reshape((n*m,1))

    #Defining initial position of wave packet.
    X_0 = X[int(n*m/2 - m/2)-1]
    Y_0 = Y[int(n*m/2 - m/2)-1]

    return X, Y, X_0, Y_0


def Psi(s, kx, ky, m, n,pos):
    """
    ===========================================================================
    Creation of a 2D Gaussian wavepacket
    ===========================================================================

    Inputs
    -----------
    s : float
        Width of gaussian Wavepacket
    sigma : float
        Phase of guassian Wavepacket
    kx, ky : float
        wavenumbers along x and y directions

    m  : Number of atoms along x
    n  : Number of atoms along y


    Returns
    -----------
    Psi : ndarray (mxn,1), complex
        a matrix with the value of the wavefunction defined at each atom

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    #also need to assert that the other numbers are real ...

    #calling Crystal function to import all of the positions of atoms and
    #the initial wave packet position
    X, Y, X_0, Y_0 = pos

    #Defining and empty array of complex type for populating
    Psi = np.zeros((n*m,1),dtype=complex)

    #Calculating normalised wave packet
    Psi = (np.exp(-0.5*(((X - X_0)/s)**2 + ((Y - Y_0)/s)**2))*np.exp((kx*X + ky*Y)*1j))/np.sqrt(4*np.pi*s)

    return Psi

if __name__ == "__main__":
    pos = Crystal(100,100)
    Psi = Psi(14, 0.1, 0.1, 100, 100, pos)
    X,Y,X0,Y0 = pos
    pd = np.abs(Psi)**2
    plt.contourf(X.reshape((100,100)),Y.reshape((100,100)),pd.reshape((100,100)),100, cmap = 'gnuplot')
    plt.plot(X,Y,'bo',markersize = 0.2)
    plt.plot(X0,Y0,'ro',markersize = 0.2)
    plt.show()
