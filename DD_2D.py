"""Module for 2 D problem
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import constants
from scipy import sparse
import scipy.linalg
from scipy.sparse.linalg import expm_multiply
from scipy.sparse import dia_matrix

def TBH(DP,n=10,m=10,dt=0.1e-15, V=False):
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

            H_d[ip:fp] = H_V*DP[j,0]
            H_d[ip:fp] = H_V*DP[j,0]

            ip += n
            fp += n
    
    H_cd_u = np.full(N-1,-H_1,dtype=complex)#.reshape((n*m-1,1))
    H_cd_l = np.full(N-1,-H_1,dtype=complex)

    #need to check these ... might not actually need two of these ...
    H_cd_u[n-1:N:n] = 0
    H_cd_l[n-1:N:n] = 0
    #print(H_cd_u)
    #need to set every other one of these to zero ...

    H_rd = np.full(n*m-n,-H_1)

    H = sparse.csr_matrix(np.diag(H_d, 0) + np.diag(H_rd, -n) + np.diag(H_rd, n) + np.diag(H_cd_l, -1) + np.diag(H_cd_u, 1))

    return H

def TB_solver_2D(DP, n, m, pos, wfc, dt, DT, V=False):
    """Split operator technique for propogation of wave packets with square
        lattice tight-binding hamiltonian.

    """

    Ns = round(DT/dt)+1

    wvf2D = wfc

    H  = TBH(DP,n,m,dt, V)

    #points = Crystal(m, n)

    for i in range(Ns):

        wvf2D = expm_multiply(H, wvf2D)

    wvf_c = np.conjugate(wvf2D)

    pd = np.multiply(wvf_c, wvf2D)
    pd = np.reshape(pd,(n,m))

    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)),pd)#,cmap='RdGy'
    plt.show()

    return pd

if __name__ == "__main__":
    #Crystal(5,5)
    #TBH(5,5,dt=0.1e-15)
    n = 100
    m = 100
    dt=0.1e-15
    DT = 1e-15
    wvf2=TB_solver_2D(n,m,pos, wfc,dt,DT)
