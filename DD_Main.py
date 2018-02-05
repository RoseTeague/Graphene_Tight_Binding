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

# Import modules for graphene
# from DD_WP_G import *           # Wavefunction module
# from DD_DP_G import *           # Potential module
# from DD_GH import *             # Hamiltonian module
# from DD_TB_G import *           # Solver module

# Import modules for Square lattice
from DD_WP_S import *           # Wavefunction module
from DD_DP_S import *           # Potential module
from DD_SH import *             # Hamiltonian module
from DD_TB_S import *           # Solver module

# Import video creation module
from animate import *

# Describe system for use in image/video titles
def TB(lattice,potential):

    n=100
    m=100
    lc = 1
    s = 5*lc
    kx = math.pi/(5*lc)
    ky = math.pi/(5*lc)
    dt = 0.1e-15
    T = 7e-15

    pos = Crystal(n,m)
    wfc = Psi(s,kx,ky,n,m,pos)
    DP = oneDdisorderpotential(m,n,lc,pos)
    H = TBH(DP,n,m,dt,V=False)

    #Change this to True to produce an mp4 video
    animate = True

    if animate:
        TB_solver(n,m,pos,wfc,DP,H,T,dt,video=True)
        MakeMovie('Tight Binding in ' + lattice + ' with ' + potential + ' potential')
    else:
        TB_solver(n,m,pos,wfc,DP,H,T,dt,False)



if __name__ == '__main__':
    TB('Square', 'no')
