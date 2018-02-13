#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module that tests other modules in subdirectories
"""

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

def disorder_potential(lattice, V):
    """
    ===========================================================================
                       Checking the disorder potential
    ===========================================================================

    Function that tests the different disorder potentials and returns it

    Inputs
    -----------
    lattice - string,
            description of the lattice to be modelled. This can be either
                'square' or 'graphene'.

    V       - string,
            description of the potential to be used on the system. This can be
                either 'None', 'one dimensional' or 'two dimensional'.

    Prints
    ----------
    DP - array of the disorder potential

    """

    # System parameters: can be modified
    n = 10                                  # Number of atoms along x
    m = 10                                  # Number of atoms along x
    lc = 10                                 # Correlation length for disorder potential in Angstroms

    if lattice == 'square':
        # Import system modules for Square lattice
        from wavepackets.DD_WP_S import Crystal                                     # Wavefunction module
        pos = Crystal(m,n)

        #  Check for the potential
        if V == 'one dimensional':
            from disorder_potentials.DD_DP_S import oneDdisorderpotential           # Potential module
            DP = oneDdisorderpotential(m,n,lc,pos)

        else:
            from disorder_potentials.DD_DP_S2 import twoDdisorderpotential          # Potential module
            DP = twoDdisorderpotential(m,n,lc,pos)


    if lattice == 'graphene':
        # Import system modules for Graphene (hexagonal) lattice
        from wavepackets.DD_WP_G import Crystal                                     # Wavefunction module
        pos = Crystal(m,n)

        # Importing the potential
        if V == 'one dimensional':
            from disorder_potentials.DD_DP_G import oneDdisorderpotential           # Potential module
            DP = oneDdisorderpotential(m,n,lc,pos)

        else:
            from disorder_potentials.DD_DP_G2 import twoDdisorderpotential          # Potential module
            DP = twoDdisorderpotential(m,n,lc,pos)

    # Prinf the disorder potential as an array
    print(DP)


def hamiltonians (lattice, V, form):
    """
    ===========================================================================
                           Checking of the Hamiltonian
    ===========================================================================

    Function that tests the different hamiltonians and returns it

    Inputs
    -----------
    lattice - string,
            description of the lattice to be modelled. This can be either
                'square' or 'graphene'.

    V       - string,
            description of the potential to be used on the system. This can be
                either 'None', 'one dimensional' or 'two dimensional'.

    form    - string,
            description of the form of the hamiltonian to be used on the system.
            This can be either 'full', 'tridiagonal', by default it runs as tridiagonal.

    Prints
    ------
    H - array of the hamiltonian

    """

    # System parameters: can be modified
    n = 4                                  # Number of atoms along x
    m = 4                                  # Number of atoms along y
    lc = 5                                 # Correlation length for disorder potential in Angstroms
    dt = 0.1e-15                           # Time step in seconds

    if lattice == 'square':
        # Import system modules for Square lattice
        from wavepackets.DD_WP_S import Crystal                                 # Wavefunction module

        pos = Crystal(m,n)

        # Import the potential
        if V == 'one dimensional':
            from disorder_potentials.DD_DP_S import oneDdisorderpotential       # Potential module
            DP = oneDdisorderpotential(m,n,lc,pos)

        else:
            from disorder_potentials.DD_DP_S2 import twoDdisorderpotential      # Potential module
            DP = twoDdisorderpotential(m,n,lc,pos)

        # Import the hamiltonian for printing
        if form == 'full':
            from hamiltonians.DD_FH_S import FTBH
            H = FTBH(DP, n, m, dt, V)

        else:
            from hamiltonians.DD_SH import TBH
            H = TBH(DP, n, m, dt, V)


    if lattice == 'graphene':
        # Import system modules for Graphene (hexagonal) lattice
        from wavepackets.DD_WP_G import Crystal                                 # Wavefunction module
        pos = Crystal(m,n)

        #  Import the potential
        if V == 'one dimensional':
            from disorder_potentials.DD_DP_G import oneDdisorderpotential       # Potential module
            DP = oneDdisorderpotential(m,n,lc,pos)


        else:
            from disorder_potentials.DD_DP_G2 import twoDdisorderpotential          # Potential module
            DP = twoDdisorderpotential(m,n,lc,pos)

        # Produce the hamiltonian for printing
        if form == 'full':
            from hamiltonians.DD_FH_G import FTBH
            H = FTBH(DP, n, m, dt, V)

        else:
            from hamiltonians.DD_GH import TBH
            H = TBH(DP, n, m, dt, V)

    # Print of the Hamiltonian as an array
    print(H)


def wavepackets(lattice):
    """
    ===========================================================================
                            Checking of the Wavepackets
    ===========================================================================

    Function that tests the wavepackets and plots it

    Inputs
    -----------
    lattice - string,
            description of the lattice to be modelled. This can be either
                'square' or 'graphene'.

    Plots
    -----
    Contour plot of the wavepacket on a lattice

    """

    # System parameters: can be modified
    m = 100                     # Number of atoms along x
    n = 100                     # Number of atoms along y
    lc = 1                      # Correlation length for disorder potential in Angstroms
    s = 5*lc                    # Width of initial Gaussian wavepacket in Angstroms
    kx = 1/(5*lc)               # Momentum eigenvalue (wavenumber) along x in reciprocal Angstroms
    ky = 1/(5*lc)               # Momentum eigenvalue (wavenumber) along y in reciprocal Angstroms

    if lattice == 'square':
        # Import system modules for Square lattice
        from wavepackets.DD_WP_S import Crystal                         # Wavefunction module
        from wavepackets.DD_WP_S import Psi                             # Wavefunction module
        pos = Crystal(m,n)
        Psi = Psi(s, kx, ky, m, n, pos)

    if lattice == 'graphene':
        # Import system modules for Graphene (hexagonal) lattice
        from wavepackets.DD_WP_G import Crystal                         # Wavefunction module
        from wavepackets.DD_WP_G import Psi                             # Wavefunction module
        pos = Crystal(m,n)
        Psi = Psi(s, kx, ky, m, n, pos)

    X, Y = pos
    pd = np.abs(Psi)**2

    # Contour plot
    plt.contourf(X.reshape((n,m)),Y.reshape((n,m)),pd.reshape((n,m)), 100, cmap = 'gnuplot')
    plt.show()

if __name__ == '__main__':
    wavepackets('square')
    hamiltonians('square', 'one dimensional', 'full')
    disorder_potential('square', 'one dimensional')
