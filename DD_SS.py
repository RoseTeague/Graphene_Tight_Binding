#!/usr/bin/env python3
"""Module for simple 2D solver
"""

import numpy as np
from scipy.sparse.linalg import expm_multiply
import matplotlib.pyplot as plt

def TB_ss(n, m, pos, wfc, H, T, dt,video=False):
    """
    ============================================================================
                Time Propogation Operator Acting on Wave Packet
    ============================================================================

    Applying time popogation operator to wave packet. The function takes the
    full hamiltonian of the system (imported from DD_FH_G or DD_FH_S where it
    is saved in a sparse form) and the initial wave packet (from DD_WP_S or
    DD_WP_G). The solver takes the dot product of the exponential of the
    hamiltonian and the wave packet to find the wave packet at the next time
    step.

    It is advised that the following paper is consulted for guidance: H. J.
    Korsch and K. Rapedius, Computations in quantum mechanics made easy,
    Eur. Phys. J., 2016, 37, 055410.

    Inputs
    ------
    n - int,
        number of rows of carbon atoms

    m - int,
        number of columns of carbon atoms

    pos - array (n*m,1),
        positions of all atoms for plotting

    wvf - arrary (n*m,1),
        wavefunction at each atom

    H - four arrays from DD_GH or DD_SH,
        tridiagonal matrices (3,n*m) for hopping in columns and rows for split
        operator technique

    T - float or int,
        duration of calculation in seconds

    dt - float or int,
        time step in seconds

    video - boolean,
        determines if a video is produced


    Returns
    -------
    pd - array,
        vector of probability density on each atomic site

    video - mp4,
        movie of wave packet propogation

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"
    assert n % 2 == 0, "The Hamiltonian can only be constructed for an even number of rows"
    assert m % 2 == 0, "The Hamiltonian can only be constructed for an even number of columns"
    assert type(dt) is float or int, "The time step must be numeric"
    assert type(T) is float or int, "The duration of the calculation must be numeric"
    assert type(video) is bool, "video must be a boolean"

    #Determine the number of time steps.
    Ns = round(T/dt) + 1

    #Perform time propogation operator to wave packet Ns times
    for i in range(Ns):

        #Do the dot product of the exponential of the hamiltonian with the
        #wave packet
        wfc = expm_multiply(H, wfc)

        if video:
            #Calculate conjugate wave function
            wfc_c = np.conjugate(wfc)

            #Probability density function by element wise multiplication
            #Neglect imaginary part
            pd = np.multiply(wfc_c, wfc)
            pd = np.reshape(pd,(n,m))

            #Plotting
            plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)),pd, 100, cmap = 'gnuplot')
            plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
            plt.savefig('Images/'+str(i))

    #Calculate conjugate wave function
    wfc_c = np.conjugate(wfc)

    #Probability density function by element wise multiplication
    #Neglect imaginary part
    pd = np.multiply(wfc_c, wfc)
    pd = np.reshape(pd,(n,m))

    return pd
