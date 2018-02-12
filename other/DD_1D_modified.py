"""Module for 1 D problem
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import constants
from scipy import sparse
import scipy.linalg
from scipy.sparse.linalg import expm_multiply


def Crystal(n=10, m=0):
    """
    ============================================================================
    Function to create a graphene crystal
    ============================================================================

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"

    #Bond length of C-C. Working in angstroms. This is consitently done throughout.
    a = 1.42

    #Creating initial arrays for populating the x and y values of each carbon atom.
    X = np.zeros((n,1))

    #Creating a list from 0 to n for each carbon atom in a column. This is reshaped
    #such that further operations can be performed.
    x_1 = np.arange(n)
    x_1 = x_1.reshape((n,1))

    #Each column of carbon atoms have the same y values and they are all equally
    #spaced in the schosen supercell. Populating each row with the carbon atom
    #number multiplied by the spacing between them.
    X = a*x_1

    X_0 = a*n*0.5

    #need to write some unit tests for this ...

    return X, X_0

def Psi(s, kx, ky, n, m, pos):
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

    #calling Crystal function to import all of the positions of carbon atoms and
    #the initial wave packet position
    X, X_0 = pos

    #Defining and empty array of complex type for populating
    Psi = np.zeros((n,1),dtype=complex)

    #Calculating normalised wave packet
    Psi = np.sqrt(0.5/(np.pi*s**2))*(np.exp(-0.5*(((X - X_0)/s)**2))*np.exp(kx*X*1j))

    return Psi


def TBH(DP,n=10,m=0,dt=0.1e-15,V=False):
    """
    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"

    #Hopping integral in eV puts matrix elements in with sensible numbers
    HI = -2.7

    #Constants that populate the off diagonal elements of the hamiltonian
    H_1 = (HI*1j*constants.e*dt)/constants.hbar

    #Constructing the set of tridiagonal matrices for H

    #The matrix for n eigenstates
    N=np.arange(n,dtype=complex)

    #Parameters
    d=2*np.pi
    delta=1
    F=0.05


    M = np.full(n-1,-H_1)

    #The final Hamiltonian
    H = d*F*np.diag(N, 0) + (delta/4)*np.diag(M, -1) + np.diag(M, 1)


    return H

def TB_solver_S(n,m,pos,wfc,H,T,dt,video):
    """Split operator technique for propogation of wave packets with square
        lattice tight-binding hamiltonian.

    """

    Ns = int(T/dt)
    #n = 100
    #m=0

    #wvf = Psi(5, 1, n)

    #H  = TBH(n,dt=0.1e-15)

    for i in range(Ns):

        wfc = expm_multiply(H, wfc)

        if video:
            wfc_c = np.conjugate(wfc)

            pd = np.multiply(wfc_c, wfc)
            plt.close()
            plt.plot(np.arange(n), pd, 'b-')
            plt.savefig('Images/'+str(i))

    wfc_c = np.conjugate(wfc)

    pd = np.multiply(wfc_c, wfc)

    plt.close()
    plt.plot(np.arange(n), pd, 'b')
    plt.show()


