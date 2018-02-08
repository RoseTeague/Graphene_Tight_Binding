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

# Describe system for use in image/video titles
def TB(lattice,V):

    n = 50#500
    m = 0#500
    lc = 1
    s = 5*lc
    kx = 1#math.pi/(5*lc)
    ky = 0#math.pi/(5*lc)
    dt = 0.1e-15
    T = 5e-15
    Ns = round(T/dt) + 1

    #We should just be able to import the solver ... there does not need to be
    #Any specification on which one it is, since both are the same ...


    if lattice == 'square':
        # Import modules for Square lattice
        from DD_WP_S import Crystal                         # Wavefunction module
        from DD_WP_S import Psi                             # Wavefunction module
        from DD_DP_S import oneDdisorderpotential           # Potential module
        from DD_DP_S2 import twoDdisorderpotential          # Potential module
        from DD_SH import TBH                               # Hamiltonian module
        from DD_TB_S import TB_solver_S                         # Solver module


    if lattice == 'graphene':
        # Import modules for Square lattice
        from DD_WP_G import Crystal                         # Wavefunction module
        from DD_WP_G import Psi                             # Wavefunction module
        from DD_DP_G import oneDdisorderpotential           # Potential module
        from DD_DP_G2 import twoDdisorderpotential          # Potential module
        from DD_GH import TBH                               # Hamiltonian module
        from DD_TB_S import TB_solver_S                         # Solver module


    if lattice == '1D square':
        # Import modules for 1D square lattice
        from DD_1D_modified import Crystal
        from DD_1D_modified import Psi
        from DD_1D_modified import TBH
        from DD_1D_modified import TB_solver_S

        def oneDdisorderpotential(m,n,lc,pos):
            return 1

        #Need to have some assert statements here ...

    pos = Crystal(n,m)
    wfc = Psi(s,kx,ky,n,m,pos)
    #this means that we calcualte it even if we do not need it ... I don't like that
    #Need to add it so that the two dimensional disorder potential can be chosen here ...

    #Choose the 1D potential first ... if it is not correct, then write over it ...
    DP = oneDdisorderpotential(m,n,lc,pos)
    if V == 'two dimensional':
        DP = twoDdisorderpotential(m,n,lc,pos)

    H = TBH(DP,n,m,dt,V)

    #Change this to True to produce an mp4 video
    animate = True

    if animate:

        from animate import MakeMovie

        pd = TB_solver_S(n,m,pos,wfc,H,T,dt,animate)
        MakeMovie('Tight Binding in ' + lattice + ' with ' + V + ' potential')
    else:

        pd = TB_solver_S(n,m,pos,wfc,H,T,dt,animate)

    #Plotting
    # plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), pd, 100, cmap = 'gnuplot')
    # plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
    # plt.show()

def TBS(lattice,V):
    """
    """

    n = 100
    m = 100
    lc = 2
    s = 5*lc
    kx = math.pi/(5*lc)
    ky = math.pi/(5*lc)
    dt = 0.1e-15
    T = 4e-15
    Ns = round(T/dt) + 1

    from DD_SS import TB_ss                                 # Solver module

    if lattice == 'square':
        # Import modules for Square lattice
        from DD_WP_S import Crystal                         # Wavefunction module
        from DD_WP_S import Psi                             # Wavefunction module
        from DD_DP_S import oneDdisorderpotential           # Potential module
        from DD_DP_S2 import twoDdisorderpotential          # Potential module
        from DD_FH_S import FTBH                             # Hamiltonian module

    if lattice == 'graphene':
        # Import modules for Square lattice
        from DD_WP_G import Crystal                         # Wavefunction module
        from DD_WP_G import Psi                             # Wavefunction module
        from DD_DP_G import oneDdisorderpotential           # Potential module
        from DD_DP_G2 import twoDdisorderpotential          # Potential module
        from DD_FH_G import FTBH                             # Hamiltonian module

    pos = Crystal(n,m)
    wfc = Psi(s,kx,ky,n,m,pos)
    DP = oneDdisorderpotential(m,n,lc,pos)

    if V == 'two dimensional':
        DP = twoDdisorderpotential(m,n,lc,pos)

    FH = FTBH(DP,n,m,dt,V)

    #Change this to True to produce an mp4 video
    animate = True

    if animate:

        from animate import MakeMovie

        pd = TB_ss(n,m,pos,wfc,FH,T,dt,animate)
        MakeMovie('Tight Binding in ' + lattice + ' with ' + V + ' potential')
    else:


    #why does the solver take the DP argument?
        pd = TB_ss(n,m,pos,wfc,FH,T,dt)

    plt.contourf(pos[0].reshape((n,m)),pos[1].reshape((n,m)), pd, 100, cmap = 'gnuplot')
    plt.title('n='+str(n)+' m='+str(m)+' t='+str(Ns*0.1)+'fs')
    plt.show()

    # #Change this to True to produce an mp4 video
    # animate = False
    #
    # if animate:
    #
    #     from animate import MakeMovie
    #
    #     TB_ss(DP,n,m,pos,wfc,H,T,dt)
    #     MakeMovie('Tight Binding in ' + lattice + ' with ' + potential + ' potential')
    # else:
    #     TB_ss(DP,n,m,pos,wfc,H,T,dt)

if __name__ == '__main__':
    TBS('graphene', 'one dimensional')
