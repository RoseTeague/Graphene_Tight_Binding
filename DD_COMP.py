"""
Module for comparing the two different solvers
"""

import math
import matplotlib.pyplot as plt
import numpy as np

def Comparison(lattice,potential):
    """
    """

    n = 100
    m = 100
    lc = 1
    s = 5*lc
    kx = math.pi/(5*lc)
    ky = math.pi/(5*lc)
    dt = 0.1e-15
    T = 0#0.1e-15
    V = True

    if lattice == 'square':
        from DD_WP_S import Crystal                         # Wavefunction module
        from DD_WP_S import Psi                             # Wavefunction module
        from DD_DP_S import oneDdisorderpotential           # Potential module
        from DD_SH import TBH                               # Hamiltonian module

    if lattice = 'graphene':
        from DD_WP_G import Crystal                         # Wavefunction module
        from DD_WP_G import Psi                             # Wavefunction module
        from DD_DP_G import oneDdisorderpotential           # Potential module
        from DD_GH import TBH                               # Hamiltonian module


    pos = Crystal(n,m)
    wfc = Psi(s,kx,ky,n,m,pos)
    DP = oneDdisorderpotential(m,n,lc,pos)

    from DD_2D import TB_solver_2D                     # Solver module

    wvf2D = TB_solver_2D(DP,n, m, pos, wfc, dt, T, V)

    # Import modules for Square lattice

    from DD_TB_S import TB_solver_S                     # Solver module

    H = TBH(DP,n,m,dt,V)

    wvf_S = TB_solver_S(n,m,pos,wfc,H,T,dt, False)

    difference = wvf2D - wvf_S

    #print(np.amax(difference)/np.amax(wvf2D)*100)

    #Plotting
    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), difference)
    plt.show()

    #return wvf2D, wvf_S, differnce

if __name__ == '__main__':
    #TB('Square', 'no')
    Comparison('square', 'no')
