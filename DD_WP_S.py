"""
Module to create a 2-Dimensional Gaussian Wavepacket
"""
import numpy as np
import matplotlib.pyplot as plt


def Crystal(m=10, n=10):
    """
    ============================================================================
                    Function to create a square/rectangular crystal
    ============================================================================
    The square/rectangular crystal crystal that is created here has armchairs
    in the x direction and zig zags in the y direction. Each carbon atom has a 
    unique label {n,m}, where n refers to the row number and m to the column number.

    The coordinates of each carbon atom are created in sets of columns, since
    there is an intrinsic regularity that can be exploited in vector form.

    In later calculations, we require a vector form of all the coordinates. This
    vector is arrange as a column vector composed of each column of carbon atoms.

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
    X - arrary (n*m,1),
        vector with all of the x positions
        
    Y - arrary (n*m,1),
        vector with all of the y positions
    X_0  - float,
        initial x position for gaussian wave packet
        
    Y_0 - float,
        initial y position for gaussian wave packet

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
    Initial Gaussian wave packet distributed on the sites of each carbon atom.
    
    Takes the positions from Crystal and caculations the wavefunction on each
    atom.

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

    #Calling Crystal function to import all of the positions of atoms and
    #the initial wave packet position
    X, Y, X_0, Y_0 = pos

    #Defining and empty array of complex type for populating
    Psi = np.zeros((n*m,1),dtype=complex)

    #Calculating normalised wave packet
    Psi = (np.exp(-0.5*(((X - X_0)/s)**2 + ((Y - Y_0)/s)**2))*np.exp((kx*X + ky*Y)*1j))/np.sqrt(4*np.pi*s)

    return Psi

if __name__ == "__main__":
    n = 100
    m = 100
    kx = 1
    ky = 1
    s = 2

    
    pos = Crystal(n,m)
    Psi = Psi(s, kx, ky, m, n, pos)
    X,Y,X0,Y0 = pos
    pd = np.abs(Psi)**2
    
    plt.contourf(X.reshape((100,100)),Y.reshape((100,100)),pd.reshape((100,100)),100, cmap = 'gnuplot')
    plt.plot(X,Y,'bo',markersize = 0.2)
    plt.plot(X0,Y0,'ro',markersize = 0.2)
    plt.show()
