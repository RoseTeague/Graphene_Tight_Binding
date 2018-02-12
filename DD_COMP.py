"""
Module for comparing the two different solvers
"""

import math
import matplotlib.pyplot as plt

def Comparison(lattice,V):
    """
    ==========================================================================
    Comparison of simple and fast solvers
    ==========================================================================

    This function calculates the difference in the final wavefunctions when
    either the simple solver or the fast solver are used to propagate the
    wavepacket. This is used to validate the different approaches.

    The system parameters on lines ___ to __ can be modified to change the
    conditions on the system.

    Inputs
    ---------
    lattice - string
            - description of the lattice to be modelled. This can be either
                'square', 'graphene', or '1D square'.

    V       - string
            - description of the potential to be used on the system. This can be
                either 'None', 'one dimensional' or 'two dimensional'.
                If the lattice is '1D square', the potential MUST be 'None'

    Returns
    ---------
    N.A

    Displays a plot of the differences between the two final wavefunctions,
    found by propagating the system with either the simple or the fast solver.

    """

    # System parameters
    m = 100                         # Number of atoms along x
    n = 100                         # Number of atoms along y
    lc = 1                          # Correlation length of disorder potential
    s = 5*lc                        # Width of initial gaussian packet
    kx = math.pi/(5*lc)             # Momentum eigenvalue (wavenumber) along x
    ky = math.pi/(5*lc)             # Momentum eigenvalue (wavenumber) along y
    dt = 0.1e-15                    # Time interval to be sampled
    T = 0#0.1e-15                   # Total time of wavepacket propagation
    Ns = round(T/dt) + 1            # Integer number of time steps

    if lattice == 'square':
        # Import modules for square lattice
        from DD_WP_S import Crystal                         # Wavefunction module
        from DD_WP_S import Psi                             # Wavefunction module
        from DD_DP_S import oneDdisorderpotential           # Potential module
        from DD_DP_S2 import twoDdisorderpotential          # Potential module
        from DD_SH import TBH                               # Hamiltonian module
        from DD_FH_S import FTBH                            # Hamiltonian module

    if lattice == 'graphene':
        # Import modules for graphene lattice
        from DD_WP_G import Crystal                         # Wavefunction module
        from DD_WP_G import Psi                             # Wavefunction module
        from DD_DP_G import oneDdisorderpotential           # Potential module
        from DD_DP_G2 import twoDdisorderpotential          # Potential module
        from DD_GH import TBH                               # Hamiltonian module
        from DD_FH_G import FTBH                            # Hamiltonian module

    # Run the system modules to descibe the lattice and initial wavepacket
    pos = Crystal(n,m)
    wfc = Psi(s,kx,ky,n,m,pos)

    # Run the potential module for the system to be simulated
    if V == 'None':
        # If the system is to be run with no external potential, assign an
        # arbitrary value of 1 to the potential variable.
        DP = 1
        V = 'no'

    if V == 'one dimensional':
        DP = oneDdisorderpotential(m,n,lc,pos)

    if V == 'two dimensional':
        DP = twoDdisorderpotential(m,n,lc,pos)

    # Import and Run the simple solver module
    from DD_SS import TB_ss                                 # Solver module
    FH = FTBH(DP,n,m,dt,V)

    wvf_ss = TB_ss(n, m, pos, wfc, FH, T, dt)

    # Import and Run the full (fast) solver module
    from DD_TB_S import TB_solver_S                         # Solver module
    H = TBH(DP,n,m,dt,V)

    wvf_S = TB_solver_S(n,m,pos,wfc,H,T,dt,False)


    # Calculates the difference between the final wavefunctions from each method.
    difference = wvf_ss - wvf_S

    #Plotting
    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), difference, 100, cmap = 'gnuplot')
    plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
    plt.show()

if __name__ == '__main__':
    #TB('Square', 'None')
    Comparison('square', 'None')
