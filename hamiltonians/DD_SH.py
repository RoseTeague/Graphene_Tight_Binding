#!/usr/bin/env python3
"""Module for Tight-Binding Hamiltonian for square lattice
"""

import numpy as np
from scipy import constants

def TBH(DP,n=10,m=10,dt=0.1e-15,V=False):
    """
    ============================================================================
          Tight Binding Hamiltonian of Square Lattice in Tridiagonal Form
    ============================================================================

    Function that creates Hamiltonian of square lattice for the split operator
    technique. These matrices are prepared in tridiagonal form for efficient
    calculations in the split operator technique.

    In a square lattice the hopping up and down a column is identical to hopping
    left and right along a row. Hence, the tridiagonal matrices have the same
    form for the off diagonal elements, but can differ in the diagonal elements
    in the presence of an external potential.

    The diagonal elements contain the information of the atoms in the square
    lattice. The reference energy of the orbtial is taken to be zero, for
    simplicity. If there is an external potential on an atom in the lattice,
    it will appear on the diagonal. If there is no external potential, then
    thers is only a 1 along the diagonal in the split operator technique.

    The off diagonal elements describe the hopping between atoms in the lattice.
    Hopping along columns and rows is the same in a square lattice, since each
    atom has hopping contributions up, down, left and right. Hence, there are
    blocks of hopping integrals followed by a zero to describe the communication
    between nearest neighbours.

    Inputs
    ------
    n - int,
        number of rows of carbon atoms

    m - int,
        number of columns of carbon atoms

    dt - float or int,
        time step in seconds

    V - string,
        determines wha type of external potential to include

    Parameters
    ----------
    HI - float,
        hopping integral of square lattice

    Returns
    -------
    TH1P - 3x(n*m) matrix,
        tridiagonal matrix for H_m of split operator technique. Positive version
        in the linear equation.

    TH1N - 3x(n*m) matrix,
        tridiagonal matrix for H_m of split operator technique. Negative version
        in the linear equation.

    TH2P - 3x(n*m) matrix,
        tridiagonal matrix for H_n of split operator technique. Positive version
        in the linear equation.

    TH2N - 3x(n*m) matrix,
        tridiagonal matrix for H_n of split operator technique. Negative version
        in the linear equation.

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert type(V) is str, "The specification for the potential must be a string"
    assert n % 2 == 0, "The Hamiltonian can only be constructed for an even number of rows"
    assert m % 2 == 0, "The Hamiltonian can only be constructed for an even number of columns"
    assert type(dt) is float or int, "The time step must be numeric"

    #Total number of carbon atoms.
    N = n*m

    #Hopping integral in eV puts matrix elements in with sensible numbers.
    HI = -2.7

    #Constants that populate the off diagonal elements of the hamiltonian.
    H_1 = (HI*1j*constants.e*dt*0.25)/constants.hbar
    H_2 = (HI*1j*constants.e*dt*0.5)/constants.hbar

    #Constant that multiplies the disorder potential on the diagonal element.
    H_V = 1j*constants.e*dt*0.5/constants.hbar

    #Constructing the set of tridiagonal matrices for H_m.
    TH1P = np.zeros((3,N),dtype=complex)
    TH1N = np.zeros((3,N),dtype=complex)

    #Setting diagonal elements of TH1P and TH1N. If there is a disorder potenital
    #then it is included here.
    if V == 'one dimensional' or V == 'two dimensional':

        TH1P[1,:] = 1 + 0.25*H_V*DP.reshape((1,N))
        TH1N[1,:] = 1 - 0.25*H_V*DP.reshape((1,N))

    else:

        TH1P[1] = 1
        TH1N[1] = 1

    #Setting off-diagonal elements. Most of these are equal to the hopping integral
    #term, so we set all of them to that initially.
    TH1P[0] = H_1
    TH1N[0] = -H_1
    TH1P[2] = H_1
    TH1N[2] = -H_1

    #Then reset the ones that should be zero because different columns do not
    #talk to eachother in this matrix.
    TH1P[0,0:N:n] = 0
    TH1N[0,0:N:n] = 0
    TH1P[2,n-1:N:n] = 0
    TH1N[2,n-1:N:n] = 0

    #Construct the set of tridiagonal matrices for H_n.
    TH2P = np.zeros((3,N),dtype=complex)
    TH2N = np.zeros((3,N),dtype=complex)

    #Setting diagonal elements of TH2P and TH2N. If there is a disorder potenital
    #then it is included here.
    if V == 'one dimensional' or V == 'two dimensional':

        DP = DP.reshape((n,m)).T.reshape((N,1))
        TH2P[1,:] = 1 + H_V*DP.reshape((1,N))
        TH2N[1,:] = 1 - H_V*DP.reshape((1,N))

    else:

        TH2P[1] = 1
        TH2N[1] = 1

    #Setting off-diagonal elements of m matrix in the same way as n.
    TH2P[0] = H_2
    TH2N[0] = -H_2
    TH2P[2] = H_2
    TH2N[2] = -H_2

    TH2P[0,0:N:m] = 0
    TH2N[0,0:N:m] = 0
    TH2P[2,m-1:N:m] = 0
    TH2N[2,m-1:N:m] = 0

    #return matrix forms of hamiltonians for split operator technique.
    return TH1P, TH1N, TH2P, TH2N