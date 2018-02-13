#!/usr/bin/env python3
"""
Module to create a 2-Dimensional Gaussian Wavepacket
"""
import numpy as np
import matplotlib.pyplot as plt


def Crystal(m, n):
    """
    ============================================================================
                    Function to create a square/rectangular crystal
    ============================================================================

    A square/rectangular lattice. Each atom has a unique label {n,m}, where n
    refers to the row number and m to the column number.

    In later calculations, we require a vector form of all the coordinates. This
    vector is arrange as a column vector composed of each column of atoms.

    Inputs
    -----------
    m - integer,
        Number of atoms along the x-direction

    n - integer,
        Numver of atoms along the y-direction

    Parameters
    -----------
    a - float,
        Lattice parameter

    b - float,
        Lattice parameter. By default this is set to a to study a square lattice

    Returns
    ----------
    X - arrary (n*m,1),
        vector with all of the x positions

    Y - arrary (n*m,1),
        vector with all of the y positions
    X_0  - float,
        initial x position for gaussian wave packet

    Y_0 - float,
        initial y position for gaussian wave packet

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"

    #Lattice parameters in the y and x direction, respectively.
    a = 1.42
    b = 1.42#Can make rectangular lattice by changing b

    #Creating initial arrays for populating the x and y values of each carbon atom.
    X = np.zeros((m,n))
    Y = np.zeros((n,m))

    #Creating a list from 0 to n for each atom in a column. This is reshaped
    #such that further operations can be performed.
    y_1 = np.arange(n)
    y_1 = y_1.reshape((n,1))

    #Each column has regularly spaced atoms.
    Y[:,0:m] = - a*y_1

    #Reshaping for further calculations
    Y_T = Y.T
    Y = Y_T.reshape((n*m,1))

    #Creating a list from 0 to n for each atom in a column. This is reshaped
    #such that further operations can be performed.
    x_1 = np.arange(m)
    x_1 = x_1.reshape((m,1))
    X[:,] = b*x_1

    #Reshaping for further calculations
    X = X.reshape((n*m,1))

    #Defining initial position of wave packet.
    X_0 = X[int(n*m/2 - m/2)-1]
    Y_0 = Y[int(n*m/2 - m/2)-1]

    return X, Y, X_0, Y_0


def Psi(s, kx, ky, m, n,pos):
    """
    ===========================================================================
                        Creation of a 2D Gaussian wavepacket
    ===========================================================================
    Initial Gaussian wave packet distributed on the sites of each atom.

    Takes the positions from Crystal and caculations the wavefunction on each
    atom.

    Inputs
    -----------
    s - float,
        Width of gaussian Wavepacket

    kx, ky - float,
        wavenumbers along x and y directions

    m  - int,
        Number of atoms along x

    n  - int,
        Number of atoms along y

    pos - array,
        vector of all atomic positions in square lattice

    Returns
    -----------
    Psi : ndarray (mxn,1), complex
        a matrix with the value of the wavefunction defined at each atom

    """

    assert type(n) is int, "Initial number of rows of carbon atoms must be an integer"
    assert type(m) is int, "Initial number of columns of carbon atoms must be an integer"

    #Calling Crystal function to import all of the positions of atoms and
    #the initial wave packet position
    X, Y, X_0, Y_0 = pos

    #Defining and empty array of complex type for populating
    Psi = np.zeros((n*m,1),dtype=complex)

    #Calculating normalised wave packet
    Psi = (np.exp(-0.5*(((X - X_0)/s)**2 + ((Y - Y_0)/s)**2))*np.exp((kx*X + ky*Y)*1j))/np.sqrt(4*np.pi*s)

    return Psi
