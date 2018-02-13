#!/usr/bin/env python3
"""
Module for comparing the two different solvers
"""

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#Defining the plotting stuff ... this was taken from the following websites. It just puts figures in a nice LaTeX format.

#http://bkanuka.com/articles/native-latex-plots/
#http://sbillaudelle.de/2015/02/23/seamlessly-embedding-matplotlib-output-into-latex.html

mpl.use('pgf')

def figsize(scale):
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1)/2.0              # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "text.fontsize": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(1),       # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

#End of plotting stuff


def Comparison(lattice,V):
    """
    ==========================================================================
                    Comparison of Simple and Fast Solvers
    ==========================================================================

    This function calculates the difference in the final wavefunctions when
    either the simple solver (DD_SS.py) or the fast solver (DD_TB_S.py) are
    used to propagate the wavepacket. This is used to validate the fast approach.

    The system parameters on lines __ to __: can be modified to change the
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
    T = 0.1e-15                     # Total time of wavepacket propagation
    Ns = round(T/dt) + 1            # Integer number of time steps

    if lattice == 'square':
        # Import modules for square lattice
        from wavepackets.DD_WP_S import Crystal                                 # Wavefunction module
        from wavepackets.DD_WP_S import Psi                                     # Wavefunction module
        from disorder_potentials.DD_DP_S import oneDdisorderpotential           # Potential module
        from disorder_potentials.DD_DP_S2 import twoDdisorderpotential          # Potential module
        from hamiltonians.DD_SH import TBH                                      # Hamiltonian module
        from hamiltonians.DD_FH_S import FTBH                                   # Hamiltonian module

    if lattice == 'graphene':
        # Import modules for graphene lattice
        from wavepackets.DD_WP_G import Crystal                                 # Wavefunction module
        from wavepackets.DD_WP_G import Psi                                     # Wavefunction module
        from disorder_potentials.DD_DP_G import oneDdisorderpotential           # Potential module
        from disorder_potentials.DD_DP_G2 import twoDdisorderpotential          # Potential module
        from hamiltonians.DD_GH import TBH                                      # Hamiltonian module
        from hamiltonians.DD_FH_G import FTBH                                   # Hamiltonian module

    # Run the system modules to descibe the lattice and initial wavepacket
    pos = Crystal(n,m)
    wfc = Psi(s,kx,ky,n,m,pos)

    # Run the potential module for the system to be simulated
    if V == 'None':
        # If the system is to be run with no external potential, assign an
        # arbitrary value of 1 to the potential variable.
        DP = 1

    if V == 'one dimensional':
        DP = oneDdisorderpotential(m,n,lc,pos)

    if V == 'two dimensional':
        DP = twoDdisorderpotential(m,n,lc,pos)

    # Import and Run the simple (VERY SLOW!) solver module
    from solvers.DD_SS import TB_ss                                             # Solver module
    FH = FTBH(DP,n,m,dt,V)

    wvf_ss = TB_ss(n, m, pos, wfc, FH, T, dt)

    # Import and Run the full (FAST) solver module
    from solvers.DD_TB_S import TB_solver_S                                     # Solver module
    H = TBH(DP,n,m,dt,V)

    wvf_S = TB_solver_S(n,m,pos,wfc,H,T,dt,False)


    # Calculates the difference between the final wavefunctions from each method.
    difference = wvf_ss - wvf_S

    #Plotting
    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), difference, 100, cmap = 'gnuplot')
    plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
    plt.show()

if __name__ == '__main__':
    Comparison('square', 'None')
