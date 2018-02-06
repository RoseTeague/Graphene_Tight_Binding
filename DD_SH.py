"""Module for Tight-Binding Hamiltonian for square lattice
"""

import numpy as np
from scipy import constants
from DD_DP_S import oneDdisorderpotential

def TBH(DP,n=10,m=10,dt=0.1e-15,V=False):
    """Function that creates Hamiltonians of square lattice for the split operator
    technique. The number of rows and columns is all that is required to
    construct the matrices. These matrices are prepared in tridiagonal form
    for efficient calculations in the split operator technique.

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


    For further details of matrices of split operator technique see ...
    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"

    #Total number of carbon atoms
    N = n*m

    #Hopping integral in eV puts matrix elements in with sensible numbers
    HI = -2.7#should really define two hopping integrals ... one in the x direction and another in the y ...

    #it seems to be insensitive to the numerical value here ... this is concerning!
    H_1 = (HI*1j*constants.e*dt*0.25)/constants.hbar
    H_2 = (HI*1j*constants.e*dt*0.5)/constants.hbar

    H_V = 1j*constants.e*dt*0.5/constants.hbar#need to check this ...

    #Constructing the set of tridiagonal matrices for H_m
    TH1P = np.zeros((3,N),dtype=complex)
    TH1N = np.zeros((3,N),dtype=complex)

    #Setting diagonal elements

    #Need to specify if there is an external potential
    if V:

        #call the potential function ...

#        DP = oneDdisorderpotential(m,n)

        #Initial counters are required to populate Hamiltonian
        ip = 0
        fp = n

        #Need to loop over each column of carbon atoms to set the disorder potential
        for j in range(m):

            TH1P[1,ip:fp] = 1 + 0.25*H_V*DP[j,0]
            TH1N[1,ip:fp] = 1 - 0.25*H_V*DP[j,0]

            ip += n
            fp += n
    else:

        TH1P[1] = 1
        TH1N[1] = 1

    #Setting off-diagonal elements
    TH1P[0] = H_1
    TH1N[0] = -H_1
    TH1P[2] = H_1
    TH1N[2] = -H_1

    TH1P[0,0:N:n] = 0
    TH1N[0,0:N:n] = 0
    TH1P[2,n-1:N:n] = 0
    TH1N[2,n-1:N:n] = 0

    #Construct the set of tridiagonal matrices for H_n
    TH2P = np.zeros((3,N),dtype=complex)
    TH2N = np.zeros((3,N),dtype=complex)

    #Setting diagonal elements

    #Need to specify if there is an external potential
    if V:

        #Initial counters are required to populate Hamiltonian
        ip = 0
        fp = m

        #Need to loop over each column of carbon atoms to set the disorder potential
        for j in range(n):

            TH2P[1,ip:fp] = 1 + H_V*DP[0:m,0]
            TH2N[1,ip:fp] = 1 - H_V*DP[0:m,0]

            ip += m
            fp += m
    else:

        TH2P[1] = 1
        TH2N[1] = 1


    #Setting off-diagonal elements of m matrix
    TH2P[0] = H_2
    TH2N[0] = -H_2
    TH2P[2] = H_2
    TH2N[2] = -H_2

    TH2P[0,0:N:m] = 0
    TH2N[0,0:N:m] = 0
    TH2P[2,m-1:N:m] = 0
    TH2N[2,m-1:N:m] = 0

    #also need to work out if we can do periodic boundary conditions ...

    #return matrix forms of hamiltonians for split operator technique

    return TH1P, TH1N, TH2P, TH2N

if __name__ == "__main__":
    TBH(5,5,dt=0.1e-15,V=True)
