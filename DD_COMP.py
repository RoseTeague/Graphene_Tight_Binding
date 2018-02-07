"""
Module for comparing the two different solvers
"""

import math
import matplotlib.pyplot as plt

def Comparison(lattice,V):
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
    V = 'no'#'one dimensional'
    Ns = round(T/dt) + 1

    if lattice == 'square':
        from DD_WP_S import Crystal                         # Wavefunction module
        from DD_WP_S import Psi                             # Wavefunction module
        from DD_DP_S import oneDdisorderpotential           # Potential module
        from DD_DP_S2 import twoDdisorderpotential          # Potential module
        from DD_SH import TBH                               # Hamiltonian module        
        from DD_FH_S import FTBH                            # Hamiltonian module

    if lattice == 'graphene':
        from DD_WP_G import Crystal                         # Wavefunction module
        from DD_WP_G import Psi                             # Wavefunction module
        from DD_DP_G import oneDdisorderpotential           # Potential module
        from DD_DP_G2 import twoDdisorderpotential          # Potential module 
        from DD_GH import TBH                               # Hamiltonian module
        from DD_FH_G import FTBH                            # Hamiltonian module

    pos = Crystal(n,m)
    wfc = Psi(s,kx,ky,n,m,pos)
    DP = oneDdisorderpotential(m,n,lc,pos)
    
    if V == 'two dimensional':
        DP = twoDdisorderpotential(m,n,lc,pos)


    from DD_SS import TB_ss                                 # Solver module
    FH = FTBH(DP,n,m,dt,V)

    wvf_ss = TB_ss(n, m, pos, wfc, FH, T, dt)

    # Import modules for Square lattice

    from DD_TB_S import TB_solver_S                         # Solver module

    H = TBH(DP,n,m,dt,V)

    wvf_S = TB_solver_S(n,m,pos,wfc,H,T,dt,False)

    difference = wvf_ss - wvf_S

    #print(np.amax(difference)/np.amax(wvf2D)*100)

    #Plotting ... what does the 100 do?
    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), difference, 100, cmap = 'gnuplot')
    plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
    plt.show()

if __name__ == '__main__':
    #TB('Square', 'no')
    Comparison('square', 'no')
