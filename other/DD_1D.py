#!/usr/bin/env python3
"""Module for 1D problem
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

def Crystal(n=10):
    """
    ============================================================================
                            Function to a 1D crystal
    ============================================================================

    Function that creates a 1D lattice.

    Inputs
    ------

    n - integer,
        Number of atoms along the y-direction

    Returns
    --------
    X - arrary (n,1),
        vector with all of the x positions

    X_0  - float,
        initial x position for gaussian wave packet

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"

    #Bond length of 1D lattice. Working in angstroms.
    #This is consitently done throughout.
    a = 1.42

    #Creating initial arrays for populating the x values of each atom.
    X = np.zeros((n,1))

    #Creating a list from 0 to n for each atom in a column. This is reshaped
    #such that further operations can be performed.
    x_1 = np.arange(n)
    x_1 = x_1.reshape((n,1))

    #Populating 1D positions array
    X = a*x_1

    #Finding initial position
    X_0 = a*n*0.5

    return X, X_0

def Psi(s, kx, n, pos):
    """
    ===========================================================================
                        Creation of a 2D Gaussian wavepacket
    ===========================================================================

    Initial Gaussian wave packet distributed on the sites of a 1D lattice.

    Inputs
    -----------
    s - float,
        Width of gaussian Wavepacket

    kx - float,
        wavenumbers along x and y directions

    n - integer,
        Number of atoms along y

    pos - array,
        vector of all atomic positions in square lattice

    Returns
    -----------
    Psi - ndarray (mxn,1), complex
        a matrix with the value of the wavefunction defined at each carbon atom

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(s) is int or float, "The width of the wave packet must be numeric"
    assert type(kx) is int or float, "The initial wave vector must be numeric"
    assert type(pos) is tuple, "The vector of lattice positions and initial position of wave packet comes in a tuple"

    #Assigning crystal points and initial position
    X, X_0 = pos

    #Defining and empty array of complex type for populating
    Psi = np.zeros((n,1),dtype=complex)

    #Calculating normalised wave packet
    Psi = np.sqrt(0.5/(np.pi*s**2))*(np.exp(-0.5*(((X - X_0)/s)**2))*np.exp(kx*X*1j))

    return Psi

def TBH(n=10,dt=0.1e-15):
    """
    ============================================================================
                    Full Tight Binding Hamiltonian of 1D Lattice
    ============================================================================

    Function that creates full hamiltonian of a 1D lattice to apply in the time
    propogation operator.

    It is advised that the following paper is consulted for guidance: H. J.
    Korsch and K. Rapedius, Computations in quantum mechanics made easy,
    Eur. Phys. J., 2016, 37, 055410.

    Inputs
    ------
    n - int,
        number of rows of carbon atoms

    dt - float or int,
        time step in seconds

    Parameters
    ----------
    HI - float,
        hopping integral of carbon in graphene

    Returns
    -------
    H - sparse matrix,
        Hamiltonian of graphene saved in sparse form

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(dt) is float or int, "The time step must be numeric"

    #Hopping integral in eV puts matrix elements in sensible numbers
    HI = -2.7

    #Constants that populate the off diagonal elements of the hamiltonian
    H_1 = (HI*1j*constants.e*dt)/constants.hbar

    #Constructing the set of tridiagonal matrices for H

    M = np.full(n-1,-H_1)

    #The final Hamiltonian
    H = sparse.csr_matrix(np.diag(M, -1) + np.diag(M, 1))

    return H

def TB_solver_S(n,wfc,H,T,dt,video):
    """
    ============================================================================
                Time Propogation Operator Acting on Wave Packet
    ============================================================================

    Applying time popogation operator to 1D wave packet. The function takes the
    full hamiltonian of the system and the initial wave packet. The solver takes
    the dot product of the exponential of the hamiltonian and the wave packet
    to find the wave packet at the next time step.

    It is advised that the following paper is consulted for guidance: H. J.
    Korsch and K. Rapedius, Computations in quantum mechanics made easy,
    Eur. Phys. J., 2016, 37, 055410.

    Inputs
    ------
    n - int,
        number of rows of lattice sites

    wfc - arrary (n,1),
        wavefunction at each atom

    H - sparse array,
        hamiltonian of 1D lattice

    T - float or int,
        duration of calculation in seconds

    dt - float or int,
        time step in seconds

    video - boolean,
        determines if a video is produced

    Returns
    -------
    pd - array,
        vector of probability density on each atomic site

    video - mp4,
        movie of wave packet propogation

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(dt) is float or int, "The time step must be numeric"
    assert type(T) is float or int, "The duration of the calculation must be numeric"
    assert type(video) is bool, "video must be a boolean"

    #determining the number of time steps
    Ns = int(T/dt)

    #looping over all time steps
    for i in range(Ns):

        #calculating wave packet at next time step
        wfc = expm_multiply(H, wfc)

        #producing video
        if video:

            wfc_c = np.conjugate(wfc)

            pd = np.multiply(wfc_c, wfc)
            plt.close()
            plt.plot(np.arange(n), pd, 'b-')
            plt.savefig('Images/'+str(i))

    #calculating conjugate of wave packet
    wfc_c = np.conjugate(wfc)

    #calculating probability density of wave packet
    pd = np.multiply(wfc_c, wfc)

    #plotting
    plt.close()
    plt.plot(np.arange(n), pd, 'b')
    plt.show()

    return pd
