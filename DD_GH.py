#!/usr/bin/env python3
"""Module for Tight-Binding Hamiltonian of Graphene
"""

import numpy as np
from scipy import constants

def TBH(DP,n=10,m=10,dt=0.1e-15,V='no'):
    """
    ============================================================================
            Tight Binding Hamiltonian of Graphene in Tridiagonal Form
    ============================================================================

    Function that creates Hamiltonian of graphene. These matrices are prepared
    in tridiagonal form for efficient calculations in the split operator technique.

    It is advised that the following paper is consulted for guidance: A. Chaves,
    L. Covaci, Kh. Yu. Rakhimov, G. A. Farias and F. M. Peeters, Wave packet
    dynamics and valley filter in strained graphene, Phys. Rev. B, 82, 2010,
    205430. DOI:https://doi.org/10.1103/PhysRevB.82.205430

    The full graphene Hamiltonian is split into two contriubtions: one from
    hopping up and down columns of carbon atoms (TH1N and TH1P), and another
    from hopping left and right along rows of carbon atoms (TH2N and TH2P).
    Both of these can be written in a tridiagonal form if the wavefunction
    is reshuffeled between calculations.

    The main diagonal of these tridagonal matrices describe the on-site energy
    and any external potentials that are present. The reference energy of the
    p-orbital of carbon in graphene is set to zero here for simplicity. In the
    absence of a potential, then, there is only a 1 along the diagonal.

    The off diagonal elements contain the hopping integral contributions. For
    columns it is relatively simple: each carbon atom, apart from the first and
    last in the column, has a carbon atom above and below it. Hence, there
    are blocks of hopping integral terms, of length n-1 on the off diagonal.
    Between these blocks there is a zero since the last element of one column
    does not talk to the first element of the next column. For the rows the
    structure of the off diagonal elements appear in these blocks again, but
    every other element of the matrix is zero, owing to the hexagonal structure
    of graphene.

    Inputs
    ------
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
    H_V = (1j*constants.e*dt*0.5)/constants.hbar

    #Constructing the set of tridiagonal matrices for H_m
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

    #Setting diagonal elements.
    if V == 'one dimensional' or V == 'two dimensional':

        DP = DP.reshape((n,m)).T.reshape((N,1))

        TH2P[1,:] = 1 + H_V*DP.reshape((1,N))
        TH2N[1,:] = 1 - H_V*DP.reshape((1,N))

    else:

        TH2P[1] = 1
        TH2N[1] = 1

    #Setting off-diagonal elements of m matrix

    #Counters for initial and final position for off-diagonal elements
    ip = 1
    fp = m

    #looping over each row.
    for i in range(n):

        if i % 2 == 0:

            #Populate off-diagonal elements for odd rows.
            TH2P[0,ip:fp:2] = H_2
            TH2N[0,ip:fp:2] = -H_2
            TH2P[2,ip-1:fp-1:2] = H_2
            TH2N[2,ip-1:fp-1:2] = -H_2

        else:

            #Populate off-diagonal elements for even rows.
            TH2P[0,ip+1:fp-1:2] = H_2
            TH2N[0,ip+1:fp-1:2] = -H_2
            TH2P[2,ip:fp-1:2] = H_2
            TH2N[2,ip:fp-1:2] = -H_2

        #Update initial and final points.
        ip += m
        fp += m

    #return matrix forms of hamiltonians for split operator technique.
    return TH1P, TH1N, TH2P, TH2N

if __name__ == "__main__":
    from DD_WP_G import Crystal
    from DD_DP_G import oneDdisorderpotential

    n = 4
    m = 4
    lc = 5

    pos = Crystal(m,n)
    DP = oneDdisorderpotential(m,n,lc,pos)

    TBH(DP, 5, 5, dt=0.1e-15, V='one dimensional')
