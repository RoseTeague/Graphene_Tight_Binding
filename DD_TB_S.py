#!/usr/bin/env python3
"""Module for split operator solver of with tight-binding hamiltonian
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg

def TB_solver_S(n, m, pos, wvf, H, T, dt, video=False):
    """
    ============================================================================
        Split Operator Technique for Propogation of a Wave Packet with
                            Tight-Binding Hamiltonian.
    ============================================================================

    Function that takes the wavef packet and hamiltonian (on a square or graphene
    lattice - the solver is generic) and uses the split operator technique to
    efficiently propogate the wave packet.

    It is advised that the following paper is consulted for guidance: A. Chaves,
    L. Covaci, Kh. Yu. Rakhimov, G. A. Farias and F. M. Peeters, Wave packet
    dynamics and valley filter in strained graphene, Phys. Rev. B, 82, 2010,
    205430. DOI:https://doi.org/10.1103/PhysRevB.82.205430

    Initially, the wave packet and trigdiagonal matrices that form the hamiltonain
    are imported and arrays are created for populating. For each time step the
    wave packet must be multiplied by the hamiltonian matrix (column form to start)
    and then the linear equation is solved. The resulting vector is then reshaped
    (into a row form) to permit the second multiplication and solution to the
    linear equation. Fianlly, the vector is reshaped (back to the column form),
    multiplied by a matrix again and the linear equation is solved to yield the
    wave packet at the next time step.

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

    #Number of time steps to be taken.
    Ns = int(T/dt) + 1

    #Total number of atoms in calculation.
    N = n*m

    #Importing hamiltonain from module. Exctract the different matrices.
    TH1P = H[0]
    TH1N = H[1]
    TH2P = H[2]
    TH2N = H[3]

    #Vectors of the matrices are defined so a tridiagonal matrix multiplication
    #can be performed.
    TH1N0 = TH1N[0,2:N].reshape((N-2,1))
    TH1N1 = TH1N[1,1:N-1].reshape((N-2,1))
    TH1N2 = TH1N[2,0:N-2].reshape((N-2,1))

    TH2N0 = TH2N[0,2:N].reshape((N-2,1))
    TH2N1 = TH2N[1,1:N-1].reshape((N-2,1))
    TH2N2 = TH2N[2,0:N-2].reshape((N-2,1))

    #Empty vectors for populating the intermediate `wavefunction' of the split
    #operator technique
    psi_p = np.zeros((N,1),dtype=complex)
    eta_p = np.zeros((N,1),dtype=complex)
    xi_p = np.zeros((N,1),dtype=complex)

    for i in range(Ns):

        #Matrix multiplication with wavefunction and TH1N
        psi_p[0] = TH1N[1,0]*wvf[0] + TH1N[0,1]*wvf[1]
        psi_p[-1] = TH1N[1,-1]*wvf[-1] + TH1N[2,-2]*wvf[-2]
        psi_p[1:N-1] = np.multiply(TH1N0, wvf[2:N]) + np.multiply(TH1N1, wvf[1:N-1]) + np.multiply(TH1N2, wvf[0:N-2])

        #Solving linear system of equations
        eta_c = scipy.linalg.solve_banded((1,1), TH1P, psi_p)

        #Reshape vector for next calculation
        eta_r = eta_c.reshape((n,m)).T.reshape((N,1))

        #Matrix multiplication with eta and TH2N
        eta_p[0] = TH2N[1,0]*eta_r[0] + TH2N[0,1]*eta_r[1]
        eta_p[-1] = TH2N[1,-1]*eta_r[-1] + TH2N[2,-2]*eta_r[-2]
        eta_p[1:N-1] = np.multiply(TH2N0, eta_r[2:N]) + np.multiply(TH2N1, eta_r[1:N-1]) + np.multiply(TH2N2, eta_r[0:N-2])

        #Solving linear system of equations
        xi_r = scipy.linalg.solve_banded((1,1), TH2P, eta_p)

        #Reshape vector for next calculation
        xi_c = xi_r.reshape((m,n)).T.reshape((N,1))

        #Matrix multiplication with xi and TH1N
        xi_p[0] = TH1N[1,0]*xi_c[0] + TH1N[0,1]*xi_c[1]
        xi_p[-1] = TH1N[1,-1]*xi_c[-1] + TH1N[2,-2]*xi_c[-2]
        xi_p[1:N-1] = np.multiply(TH1N0, xi_c[2:N]) + np.multiply(TH1N1, xi_c[1:N-1]) + np.multiply(TH1N2, xi_c[0:N-2])

        #Solving for wavefunction at next time step
        wvf = scipy.linalg.solve_banded((1,1), TH1P, xi_p)

        if video:
            #Calculate conjugate wave function
            wvf_conj = np.conjugate(wvf)

            #Probability density function by element wise multiplication
            #Neglect imaginary part
            pd = np.multiply(wvf_conj,wvf)
            pd = np.reshape(pd,(n,m))

            #Plotting
            plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)),pd, 100, cmap = 'gnuplot')
            plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
            plt.savefig('Images/'+str(i))


    #Calculate conjugate wave function
    wvf_conj = np.conjugate(wvf)

    #Probability density function by element wise multiplication.
    #Neglect imaginary part.
    pd = np.multiply(wvf_conj,wvf)
    pd = np.reshape(pd,(n,m))

    return pd
