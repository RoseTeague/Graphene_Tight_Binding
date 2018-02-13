#!/usr/bin/env python3
"""
Main program for studying the propagation of a Gaussian wavepacket through 2D
crystals using the Tight Binding Method.
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
    golden_mean = (np.sqrt(5.0)-1)/2.0            # Aesthetic ratio (you could change this)
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
    "figure.figsize": figsize(1),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

#End of plotting stuff


def TB(lattice,V,animate=True):
    """
    ===========================================================================
              Main file for the split operator tight binding solver
    ===========================================================================

    This function collects all modules which define the lattice, the initial
    wavepacket, the Tight Binding Hamiltonian and the full split operator
    solver. It applies the same parameters to each module and allows the program
    to be run from one place.

    The following variables are defined on lines ___ to ___:
    the size of the system, the correlation length, the
    width of the Gaussian packet, the momentum eigenvalues in x and y and the
    time intervals and length to be simulated can be adjusted. These can be
    changed to simulate different conditions.

    Inputs
    -----------
    lattice - string
            - description of the lattice to be modelled. This can be either
                'square' or 'graphene'.

    V       - string
            - description of the potential to be used on the system. This can be
                either 'None', 'one dimensional' or 'two dimensional'.

    animate - boolean
            - Determines whether or not a .mp4 file of the propagation will be
                produced. N.B this takes more time so should only be set if it
                is necessary to study the propagation.

    Returns
    ----------
    N.A

    If 'animate' is set to be true, an .mp4 file will be saved.
    A plot of the final state of the wavepacket will be displayed as a contour
    plot where the colour represents the amplitude.
    """

    # System parameters: can be modified
    m = 500                     # Number of atoms along x
    n = 500                     # Number of atoms along y
    lc = 1                      # Correlation length for disorder potential in Angstroms
    s = 5*lc                    # Width of initial Gaussian wavepacket in Angstroms
    kx = math.pi/(5*lc)         # Momentum eigenvalue (wavenumber) along x in reciprocal Angstroms
    ky = math.pi/(5*lc)         # Momentum eigenvalue (wavenumber) along y in reciprocal Angstroms
    dt = 0.1e-15                # Time step in fs
    T = 5e-15                   # Total time of wavepacket propagation in fs
    Ns = round(T/dt) + 1        # Integer number of time steps


    # Import the full (FAST) Hamiltonian solver
    from solvers.DD_TB_S import TB_solver_S                                     # Solver module

    if lattice == 'square':
        # Import system modules for Square lattice
        from wavepackets.DD_WP_S import Crystal                                 # Wavefunction module
        from wavepackets.DD_WP_S import Psi                                     # Wavefunction module
        from disorder_potentials.DD_DP_S import oneDdisorderpotential           # Potential module
        from disorder_potentials.DD_DP_S2 import twoDdisorderpotential          # Potential module
        from hamiltonians.DD_SH import TBH                                      # Hamiltonian module

    if lattice == 'graphene':
        # Import system modules for Graphene (hexagonal) lattice
        from wavepackets.DD_WP_G import Crystal                                 # Wavefunction module
        from wavepackets.DD_WP_G import Psi                                     # Wavefunction module
        from disorder_potentials.DD_DP_G import oneDdisorderpotential           # Potential module
        from disorder_potentials.DD_DP_G2 import twoDdisorderpotential          # Potential module
        from hamiltonians.DD_GH import TBH                                      # Hamiltonian module


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

    # Run the Hamiltonian module for the system
    H = TBH(DP,n,m,dt,V)

    if animate:

        # Solves the tight binding hamiltonian using a split operator method
        # described in the file DD_TB_S. Records an image after each time step and
        # compiles into an mp4 movie file.
        from animate import MakeMovie

        pd = TB_solver_S(n,m,pos,wfc,H,T,dt,animate)
        MakeMovie('Tight Binding in ' + lattice + ' with ' + V + ' potential')

    else:

        # Solves the tight binding hamiltonian using a split operator method
        # describes in the file DD_TB_S
        pd = TB_solver_S(n,m,pos,wfc,H,T,dt,animate)

    #Plotting
    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), pd, 100, cmap = 'gnuplot')
    plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
    plt.show()

def TBS(lattice,V, animate = False):
    """
    ===========================================================================
           Main file for the simple split operator tight binding solver
    ===========================================================================

    This function collects all modules which define the lattice, the initial
    wavepacket, the Tight Binding Hamiltonian and the full split operator
    solver. It applies the same parameters to each module and allows the program
    to be run from one place.

    The following variables are defined on lines ___ to ___:
    the size of the system, the correlation length, the
    width of the Gaussian packet, the momentum eigenvalues in x and y and the
    time intervals and length to be simulated can be adjusted. These can be
    changed to simulate different conditions.

    Inputs
    -----------
    lattice - string
            - description of the lattice to be modelled. This can be either
                'square', 'graphene', or '1D square'.

    V       - string
            - description of the potential to be used on the system. This can be
                either 'None', 'one dimensional' or 'two dimensional'.
                If the lattice is '1D square', the potential MUST be 'None'

    animate - boolean
            - Determines whether or not a .mp4 file of the propagation will be
                produced. N.B this takes more time so should only be set if it
                is necessary to study the propagation.

    Returns
    ----------
    N.A

    If 'animate' is set to be true, an .mp4 file will be saved.
    A plot of the final state of the wavepacket will be displayed as a contour
    plot where the colour represents the amplitude.

    """

    # System Parameters
    m = 100                     # Number of atoms along x
    n = 100                     # Number of atoms along y
    lc = 2                      # Correlation length for disorder potential in Angstroms
    s = 5*lc                    # Width of initial Gaussian wavepacket in Angstroms
    kx = math.pi/(5*lc)         # Momentum eigenvalue (wavenumber) along x in reciprocal Angstroms
    ky = math.pi/(5*lc)         # Momentum eigenvalue (wavenumber) along y in reciprocal Angstroms
    dt = 0.1e-15                # Time step in fs
    T = 4e-15                   # Total time of wavepacket propagation in fs
    Ns = round(T/dt) + 1        # Integer number of time steps


    # Import the simple Hamiltonian solver
    from solvers.DD_SS import TB_ss                                             # Solver module

    if lattice == 'square':
        # Import modules for Square lattice
        from wavepackets.DD_WP_S import Crystal                                 # Wavefunction module
        from wavepackets.DD_WP_S import Psi                                     # Wavefunction module
        from disorder_potentials.DD_DP_S import oneDdisorderpotential           # Potential module
        from disorder_potentials.DD_DP_S2 import twoDdisorderpotential          # Potential module
        from hamiltonians.DD_FH_S import FTBH                                   # Hamiltonian module

    if lattice == 'graphene':
        # Import modules for Graphene lattice
        from wavepackets.DD_WP_G import Crystal                                 # Wavefunction module
        from wavepackets.DD_WP_G import Psi                                     # Wavefunction module
        from disorder_potentials.DD_DP_G import oneDdisorderpotential           # Potential module
        from disorder_potentials.DD_DP_G2 import twoDdisorderpotential          # Potential module
        from hamiltonians.DD_FH_G import FTBH                                   # Hamiltonian module

    if lattice == '1D square':
        # Import modules for 1D square lattice
        # Must have no potential for this 1D solver
        from other.DD_1D import Crystal                                         # Crystal Module
        from other.DD_1D import Psi                                             # Wavefunction Module
        from other.DD_1D import FTBH                                            # Hamiltonian Module
        from other.DD_1D import TB_ss                                           # Solver Module

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

    # Run the Hamiltonian module for the system
    FH = FTBH(DP,n,m,dt,V)

    if animate:

        # Solves the tight binding hamiltonian using time propogation operator.
        # Records an image after each time step and compiles into an mp4 movie file.
        from animate import MakeMovie

        pd = TB_ss(n,m,pos,wfc,FH,T,dt,animate)
        MakeMovie('Tight Binding in ' + lattice + ' with ' + V + ' potential')

    else:

        # Solves the tight binding hamiltonian using time propogation operator
        pd = TB_ss(n,m,pos,wfc,FH,T,dt)

    #Plotting
    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), pd, 100, cmap = 'gnuplot')
    plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
    plt.show()


if __name__ == '__main__':
    TB('graphene', 'one dimensional')
