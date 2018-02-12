#!/usr/bin/env python3
"""Module for 1 D problem
"""

import numpy as np
from scipy import constants
from scipy import sparse

def FTBH(DP, n=10,m=10,dt=0.1e-15, V='no'):
    """
    ============================================================================
                    Full Tight Binding Hamiltonian of Graphene
    ============================================================================

    Function that creates full hamiltonian of graphene to apply in the time
    propogation operator.

    It is advised that the following paper is consulted for guidance: H. J.
    Korsch and K. Rapedius, Computations in quantum mechanics made easy,
    Eur. Phys. J., 2016, 37, 055410.

    The main diagonal of the hamiltonian contains the on-site energy, which has
    been set to zero for simplicity, and the disorder potential, if present. The
    off-diagonal elements contain hopping integarls. There are a total of 4
    diagonals off of the main diagonal. The two imediately adjacent to the
    main diagonal describe the hopping up and down a column of carbon atoms in
    graphene. While the diagonals displaced by n from the main are the ones
    associated with hopping left and right along a row of carbon atoms.

    Inputs
    ------
    DP - array,
        disorder potential imported from DD_DP_S

    n - int,
        number of rows of carbon atoms

    m - int,
        number of columns of carbon atoms

    dt - float or int,
        time step in seconds

    V - string,
        determines what type of external potential should be included

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
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert type(V) is str, "The specification for the potential must be a string"
    assert n % 2 == 0, "The Hamiltonian can only be constructed for an even number of rows"
    assert m % 2 == 0, "The Hamiltonian can only be constructed for an even number of columns"
    assert type(dt) is float or int, "The time step must be numeric"

    #Total number of carbon atoms.
    N = n*m

    #Hopping integral in eV puts matrix elements in with sensible numbers.
    HI = -2.7

    #Constants that populate the off diagonal elements of the hamiltonian
    H_1 = (HI*1j*constants.e*dt)/constants.hbar
    H_V = 1j*constants.e*dt/constants.hbar

    #Constructing the main diagonal.
    H_d = np.full(N,0,dtype=complex)

    if V == 'one dimensional' or V == 'two dimensional':

        H_d = - H_V*DP

    #Off diagonal elements adjacent to main diagonal
    H_cd = np.full(N-1,-H_1,dtype=complex)
    H_cd[n-1:N:n] = 0

    #Off diagonal elements n awayy from main diagonal
    H_rd = np.full(n*m-n,-H_1)
    H_rd[1:N:2] = 0

    #Constructing full hamiltonian and saving in sparse form
    H = sparse.csr_matrix(np.diag(H_d, 0) + np.diag(H_rd, -n) + np.diag(H_rd, n) + np.diag(H_cd, -1) + np.diag(H_cd, 1))

    return H
