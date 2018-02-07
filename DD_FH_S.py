"""Module for 2 D problem
"""

import numpy as np
from scipy import constants
from scipy import sparse

def FTBH(DP,n=10,m=10,dt=0.1e-15,V='no'):
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

    if V == 'one dimensional':

        #Initial counters are required to populate Hamiltonian
        ip = 0
        fp = n

        #Need to loop over each column of carbon atoms to set the disorder potential
        for j in range(m):

            H_d[ip:fp] = H_V*DP[j,0]
            H_d[ip:fp] = H_V*DP[j,0]

            ip += n
            fp += n
            
    elif V == 'two dimensional':
        H_d = H_V*DP

    #Need to add a part in for the two dimensional part ... file needs to be written first!

    H_cd_u = np.full(N-1,-H_1,dtype=complex)#.reshape((n*m-1,1))
    H_cd_l = np.full(N-1,-H_1,dtype=complex)

    #need to check these ... might not actually need two of these ...
    H_cd_u[n-1:N:n] = 0
    H_cd_l[n-1:N:n] = 0

    #need to set every other one of these to zero ...

    H_rd = np.full(n*m-n,-H_1)

    H = sparse.csr_matrix(np.diag(H_d, 0) + np.diag(H_rd, -n) + np.diag(H_rd, n) + np.diag(H_cd_l, -1) + np.diag(H_cd_u, 1))

    return H
