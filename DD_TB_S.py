"""Module for solver of with TB split operator
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg

def TB_solver_S(n, m, pos, wvf, H, T, dt, video=False):
    """
    ============================================================================
        Split operator technique for propogation of a wave packet with
                            tight-binding hamiltonian.
    ============================================================================

        ... it is suggested that you read ... section ... before using this
        function so that it becomes clear how the calculation is being performed
        ...

    Inputs
    ------
    n - int, number of rows of carbon atoms

    m - int, number of columns of carbon atoms

    pos -

    wvf -

    DP - ... does this actully need to take DP .. ? I do not think so ...

    H -

    T -

    dt - float, time step in seconds

    video - boolean,
        determines if there is an external potential


    Returns
    -------
    pd - array,
        vector of probability density on each atomic site

    """

    #need to add some assert statements ...
    #T needs to be a number
    #video needs to be a boolean

    #might want to check the sizes of each of the inputs ... ?

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

        #Solving for wavefunction
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

    #Probability density function by element wise multiplication
    #Neglect imaginary part
    pd = np.multiply(wvf_conj,wvf)
    pd = np.reshape(pd,(n,m))

    #Plotting
    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), pd, 100, cmap = 'gnuplot')
    plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
    plt.show()

    return pd
