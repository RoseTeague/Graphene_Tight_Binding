"""
Module to create a 2-Dimensional Gaussian Wavepacket at t=0
"""
import numpy as np

def Crystal(m=20, n=10):
    """
    ============================================================================
    Function to create a graphene crystal
    ============================================================================

    The crystal of graphene that is created here is the same as in Refs. 1 and 2.
    There are armchairs in the x direction and zig zags in the y direction. Each
    carbon atom has a unique label {n,m}, where n refers to the row number and
    m to the column number.

    The coordinates of each carbon atom are created in sets of columns, since
    there is an intrinsic regularity that can be exploited in vector form.

    In later calculations, we require a vector form of all the coordinates. This
    vector is arrange as a column vector composed of each column of carbon atoms.


    Inputs
    -----------
    m : integer
        Number of atoms along the x-direction

    n : integer
        Number of atoms along the y-direction

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

    References
    ----------
    1.

    2.
    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"

    #Bond length of C-C. Working in angstroms. This is consitently done throughout.
    a = 1.42

    #Defining changes in x and y between carbon atoms. The d_x here is the distance
    #between carbon atoms in the same column.
    d_x = a*0.5
    d_y = 0.5*a*np.sqrt(3)

    #Creating initial arrays for populating the x and y values of each carbon atom.
    X = np.zeros((m,n))
    Y = np.zeros((n,m))

    #Creating a list from 0 to n for each carbon atom in a column. This is reshaped
    #such that further operations can be performed.
    y_1 = np.arange(n)
    y_1 = y_1.reshape((n,1))

    #Each column of carbon atoms have the same y values and they are all equally
    #spaced in the schosen supercell. Populating each row with the carbon atom
    #number multiplied by the spacing between them.
    Y[:,0:m] = - d_y*y_1
    Y_T = Y.T
    Y = Y_T.reshape((n*m,1))

    #Creating a list from 0 to m for each carbon atom in a row. This is reshaped
    #such that further operations can be performed.
    x_1 = np.arange(m)
    x_1 = x_1.reshape((m,1))

    #The carbons in each column are not all the same, but in a zig zag.
    #To start, we define the middle down a column and then displace carbon atoms
    #either side of the middle.
    dx = np.zeros((1,n))
    dx[0,0:n:2] = d_x*0.5
    dx[0,1:n:2] = -d_x*0.5
    X[:,] = (d_x + a)*x_1
    X[0:m:2,:] = X[0:m:2,:] + dx
    X[1:m:2,:] = X[1:m:2,:] - dx

    #Then need to reshape to the required output form.
    X = X.reshape((n*m,1))

    #Defining initial position of wave packet.
    X_0 = 0.5*m*(d_x + a)
    Y_0 = -0.5*n*d_y

    #need to write some unit tests for this ...

    return X, Y, X_0, Y_0

def Psi(s, kx, ky, m, n):
    """
    ===========================================================================
    Creation of a 2D Gaussian wavepacket
    ===========================================================================

    Initial Gaussian wave packet distributed on the sites of each carbon atom.

    Inputs
    -----------
    s : float
        Width of Gaussian Wavepacket ... is it though?

    kx, ky : float
        wavenumbers along x and y directions

    m  : integer
        Number of atoms along x

    n  : integer
        Number of atoms along y

    Returns
    -----------
    Psi : ndarray (mxn,1), complex
        a matrix with the value of the wavefunction defined at each carbon atom

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    #also need to assert that the other numbers are real ...

    #calling Crystal function to import all of the positions of carbon atoms and
    #the initial wave packet position
    X, Y, X_0, Y_0 = Crystal(m,n)

    #Defining and empty array of complex type for populating
    Psi = np.zeros((n*m,1),dtype=complex)

    #Calculating normalised wave packet
    Psi = (np.exp(-0.5*(((X - X_0)/s)**2 + ((Y - Y_0)/s)**2))*np.exp((kx*X + ky*Y)*1j))/np.sqrt(4*np.pi*s)

    return Psi

if __name__ == "__main__":
    Crystal(10,10)
