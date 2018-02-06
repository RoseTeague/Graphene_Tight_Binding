"""Module for 1 D problem
"""

import numpy as np
from scipy import constants
from scipy import sparse

def TBH(DP, n=10,m=10,dt=0.1e-15, V=False):
    """
    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"

    N = n*m
    #Hopping integral in eV puts matrix elements in with sensible numbers
    HI = -2.7

    #Constants that populate the off diagonal elements of the hamiltonian
    H_1 = (HI*1j*constants.e*dt)/constants.hbar
    H_V = 1j*constants.e*dt/constants.hbar

    #Constructing the set of tridiagonal matrices for H_m
    #H = np.zeros((3,n),dtype=complex)

    H_d = np.full(N,0,dtype=complex)

    if V:

        #call the potential function ...
#        DP = oneDdisorderpotential(m,n)

        #Initial counters are required to populate Hamiltonian
        ip = 0
        fp = n

        #Need to loop over each column of carbon atoms to set the disorder potential
        for j in range(m):

            #The even and odd column numbers are populated from the disorder potential
            #differently, because the atoms are in a zig zag.
            if j % 2 == 0:

                H_d[ip:fp:2] = H_V*DP[2*j+1,0]
                H_d[ip:fp:2] = H_V*DP[2*j+1,0]

                H_d[ip+1:fp+1:2] = H_V*DP[2*j,0]
                H_d[ip+1:fp+1:2] = H_V*DP[2*j,0]

            else:

                H_d[ip:fp:2] = H_V*DP[2*j,0]
                H_d[ip:fp:2] = H_V*DP[2*j,0]

                H_d[ip+1:fp+1:2] = H_V*DP[2*j+1,0]
                H_d[ip+1:fp+1:2] = H_V*DP[2*j+1,0]

            #conting to the next set of atoms
            ip += n
            fp += n


    H_cd_u = np.full(N-1,-H_1,dtype=complex)#.reshape((n*m-1,1))
    H_cd_l = np.full(N-1,-H_1,dtype=complex)

    #need to check these ... might not actually need two of these ...
    #Should this be n or m ... ?
    H_cd_u[n-1:N:n] = 0
    H_cd_l[n-1:N:n] = 0
    #print(H_cd_u)
    #need to set every other one of these to zero ...

    H_rd_u = np.full(n*m-n,-H_1)
    H_rd_l = np.full(n*m-n,-H_1)

    H_cd_u[1:N:2] = 0
    H_cd_l[1:N:2] = 0

    H = sparse.csr_matrix(np.diag(H_d, 0) + np.diag(H_rd_l, -n) + np.diag(H_rd_u, n) + np.diag(H_cd_l, -1) + np.diag(H_cd_u, 1))

    return H
