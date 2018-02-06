"""
Main program for studying the propagation of a Gaussian wavepacket through 2D
crystals using the Tight Binding Method.


For animations to work, create an empty folder called Images and download cv2
with the following command:
                pip install opencv-python

Functions required for this method and their necessary output formats are:
    Crystal                 :  Purpose - describes the 2D lattice
                                Output  - X - list of coordinate x-values going down each column
                                                of the Lattice
                                        - Y - list of coordinate y-values going down each column
                                                of the Lattice
                                        - X0 - initial x value
                                        - Y0 - initial y value
    Psi                     :  Purpose - Describes the initial wavepacket
                                Output  - Psi - complex value of the wavepacket at each
                                            coordinate described by (X,Y)
    oneDdisorderpotential   :  Purpose - Define the potential
                                Output  -
    TBH                     :  Purpose - Define the Tight Binding Hamiltonian
                                Output  - TH1P -
                                        - TH2P -
                                        - TH1N -
                                        - TH2N -
"""

import math
import matplotlib.pyplot as plt
import numpy as np

def Comparison_S():

    n = 100
    m = 100
    lc = 1
    s = 5*lc
    kx = math.pi/(5*lc)
    ky = math.pi/(5*lc)
    dt = 0.1e-15
    T = 0#0.1e-15
    V = False

    from DD_WP_S import Crystal                         # Wavefunction module
    from DD_WP_S import Psi                             # Wavefunction module

    pos = Crystal(n,m)
    wfc = Psi(s,kx,ky,n,m,pos)

    from DD_2D import TB_solver_2D

    wvf2D = TB_solver_2D(n, m, pos, wfc, dt, T)

    # Import modules for Square lattice

    from DD_DP_S import oneDdisorderpotential           # Potential module
    from DD_SH import TBH                               # Hamiltonian module
    from DD_TB_S import TB_solver_S                     # Solver module

    DP = oneDdisorderpotential(m,n,lc,pos)
    H = TBH(DP,n,m,dt,V)

    wvf_S = TB_solver_S(n,m,pos,wfc,DP,H,T,dt, False)

    difference = wvf2D - wvf_S

    print(np.amax(difference)/np.amax(wvf2D)*100)

    #Plotting
    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), difference)
    plt.show()

    #return wvf2D, wvf_S, differnce

if __name__ == '__main__':
    #TB('Square', 'no')
    Comparison_S()
