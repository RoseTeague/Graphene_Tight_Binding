"""Module for 1 D problem
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import constants
from scipy import sparse
import scipy.linalg
from scipy.sparse.linalg import expm_multiply
from scipy.sparse import dia_matrix


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
    X[:,] = b*x_1

    #Reshaping for further calculations
    X = X.reshape((n*m,1))

    #Defining initial position of wave packet.
    X_0 = 0.5*m*b
    Y_0 = -0.5*n*a

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

    #calling Crystal function to import all of the positions of carbon atoms and
    #the initial wave packet position
    X, Y, X_0, Y_0 = Crystal(m, n)

    #Defining and empty array of complex type for populating
    Psi = np.zeros((n*m,1),dtype=complex)

    #Calculating normalised wave packet
    Psi = (np.exp(-0.5*((((X - X_0)/s)**2) + ((Y - Y_0)/s)**2))*np.exp((kx*X + ky*Y)*1j))/np.sqrt(4*np.pi*s)

    return Psi


def TBH(n=10,m=10,dt=0.1e-15):
    """
    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"

    N = n*m
    #Hopping integral in eV puts matrix elements in with sensible numbers
    HI = -2.7

    #Constants that populate the off diagonal elements of the hamiltonian
    H_1 = (HI*1j*constants.e*dt)/constants.hbar

    #Constructing the set of tridiagonal matrices for H_m
    #H = np.zeros((3,n),dtype=complex)

    H_d = np.full(N,0)
    H_cd_u = np.full(N-1,-H_1)#.reshape((n*m-1,1))
    H_cd_l = np.full(N-1,-H_1)

    #need to check these ... might not actually need two of these ...
    H_cd_u[n-1:N:n] = 0
    H_cd_l[n-1:N:n] = 0
    #print(H_cd_u)
    #need to set every other one of these to zero ...

    H_rd = np.full(n*m-n,-H_1)

    H = np.diag(H_d, 0) + np.diag(H_rd, -n) + np.diag(H_rd, n) + np.diag(H_cd_l, -1) + np.diag(H_cd_u, 1)
    H = sparse.csr_matrix(H)

    return H

def TB_solver_2D(n, m, dt, DT):
    """Split operator technique for propogation of wave packets with square
        lattice tight-binding hamiltonian.

    """

    Ns = round(DT/dt)+1


    wvf2D = Psi(5,1, 1, m,n)

    H  = TBH(n,m,dt=0.1e-15)

    points = Crystal(m, n)

    for i in range(Ns):

        wvf2D = expm_multiply(H, wvf2D)

    wvf_c = np.conjugate(wvf2D)

    pd = np.multiply(wvf_c, wvf2D)
    pd = np.reshape(pd,(n,m))

    plt.contourf(points[0].reshape((n,m)),points[1].reshape((n,m)),pd)#,cmap='RdGy'
    plt.show()
    
    return wvf2D

if __name__ == "__main__":
    #Crystal(5,5)
    #TBH(5,5,dt=0.1e-15)
    n = 100
    m = 100
    dt=0.1e-15
    DT = 1e-15
    wvf2=TB_solver_2D(n,m,dt,DT)
    print(wvf2D)
