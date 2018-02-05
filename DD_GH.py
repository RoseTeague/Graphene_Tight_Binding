"""Module for Tight-Binding Hamiltonian of Graphene
"""

import numpy as np
from scipy import constants
from scipy import sparse
import scipy.linalg
from DD_DP_G import oneDdisorderpotential

def TBH(n=10,m=10,dt=0.1e-15,lc=10,V=False):
    """Function that creates Hamiltonians of graphene for the split operator
    technique. The number of rows and columns is all that is required to
    construct the matrices. These matrices are prepared in tridiagonal form
    for efficient calculations in the split operator technique.

    ... need to actually explain what is going on ...

    Inputs
    ------
    n - int, number of rows of carbon atoms

    m - int, number of columns of carbon atoms

    dt - float, time step in seconds

    V - boolean, determines if there is an external potential

    Parameters
    ----------
    HI - float, hopping integral of carbon

    Returns
    -------
    TH1P - 3x(n*m) matrix, tridiagonal matrix for H_m of split operator technique.
            Positive version in the linear equation.

    TH1N - 3x(n*m) matrix, tridiagonal matrix for H_m of split operator technique
            Negative version in the linear equation.

    TH2P - 3x(n*m) matrix, tridiagonal matrix for H_n of split operator technique
            Positive version in the linear equation.

    TH2N - 3x(n*m) matrix, tridiagonal matrix for H_n of split operator technique
            Negative version in the linear equation.

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert n % 2 == 0, "The Hamiltonian can only be constructed for an even number of rows"
    assert m % 2 == 0, "The Hamiltonian can only be constructed for an even number of columns" 

    #Total number of carbon atoms
    N = n*m

    #Hopping integral in eV puts matrix elements in with sensible numbers
    HI = -2.7

    #Constants that populate the off diagonal elements of the hamiltonian
    H_1 = (HI*1j*constants.e*dt*0.25)/constants.hbar
    H_2 = (HI*1j*constants.e*dt*0.5)/constants.hbar

    #Constant that multiplies the disorder potential on the diagonal element
    H_V = (1j*constants.e*dt*0.5)/constants.hbar

    #Constructing the set of tridiagonal matrices for H_m
    TH1P = np.zeros((3,N),dtype=complex)
    TH1N = np.zeros((3,N),dtype=complex)

    #Setting diagonal elements of TH1P and TH1N
    if V:

        #call the potential function ...
        DP = oneDdisorderpotential(m,n)

        #Initial counters are required to populate Hamiltonian
        ip = 0
        fp = n

        #Need to loop over each column of carbon atoms to set the disorder potential
        for j in range(m):

            #The even and odd column numbers are populated from the disorder potential
            #differently, because the atoms are in a zig zag.
            if j % 2 == 0:

                TH1P[1,ip:fp:2] = 1 + 0.25*H_V*DP[2*j+1,0]
                TH1N[1,ip:fp:2] = 1 - 0.25*H_V*DP[2*j+1,0]

                TH1P[1,ip+1:fp+1:2] = 1 + 0.25*H_V*DP[2*j,0]
                TH1N[1,ip+1:fp+1:2] = 1 - 0.25*H_V*DP[2*j,0]

            else:

                TH1P[1,ip:fp:2] = 1 + 0.25*H_V*DP[2*j,0]
                TH1N[1,ip:fp:2] = 1 - 0.25*H_V*DP[2*j,0]

                TH1P[1,ip+1:fp+1:2] = 1 + 0.25*H_V*DP[2*j+1,0]
                TH1N[1,ip+1:fp+1:2] = 1 - 0.25*H_V*DP[2*j+1,0]

            #conting to the next set of atoms
            ip += n
            fp += n

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

    #Construct the set of tridiagonal matrices for H_n
    TH2P = np.zeros((3,N),dtype=complex)
    TH2N = np.zeros((3,N),dtype=complex)

    #Setting diagonal elements
    if V:

        ip = 0
        fp = m

        for j in range(n):

            if j % 2 == 0:

                TH2P[1,ip:fp] = 1 + H_V*DP[0:2*m:2,0]
                TH2N[1,ip:fp] = 1 - H_V*DP[0:2*m:2,0]

            else:

                TH2P[1,ip:fp] = 1 + H_V*DP[1:2*m:2,0]
                TH2N[1,ip:fp] = 1 - H_V*DP[1:2*m:2,0]

            ip += m
            fp += m

    else:

        TH2P[1] = 1
        TH2N[1] = 1

    #Setting off-diagonal elements of m matrix

    #Counters for initial and final position for off-diagonal elements
    ip = 1
    fp = m

    #looping over each row
    #also need to add a satement in for if there are an odd or even number of rows ... ?
    #because the number of non-zero elements changes ...

    for i in range(n):

        if i % 2 == 0:

            #Populate off-diagonal elements for odd rows
            TH2P[0,ip:fp:2] = H_2
            TH2N[0,ip:fp:2] = -H_2
            TH2P[2,ip-1:fp-1:2] = H_2
            TH2N[2,ip-1:fp-1:2] = -H_2

        else:

            #Populate off-diagonal elements for even rows
            TH2P[0,ip+1:fp-1:2] = H_2
            TH2N[0,ip+1:fp-1:2] = -H_2
            TH2P[2,ip:fp-1:2] = H_2
            TH2N[2,ip:fp-1:2] = -H_2

        #Update initial and final points
        ip += m
        fp += m

    #also need to work out if we can do periodic boundary conditions ...

    #return matrix forms of hamiltonians for split operator technique
    return TH1P, TH1N, TH2P, TH2N

if __name__ == "__main__":
    TBH(5,5,dt=0.1e-15,lc=10,V=True)
