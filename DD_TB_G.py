"""Module for solver of with TB split operator
"""
import math #can not remeber if this is how we do this ...
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import constants
from scipy import sparse
import scipy.linalg
from DD_GH import TBH
import DD_WP_G
from DD_WP_G import *

#Defining the plotting stuff ... this was taken from the following websites. It just puts figures in a nice LaTeX format.

#http://bkanuka.com/articles/native-latex-plots/
#http://sbillaudelle.de/2015/02/23/seamlessly-embedding-matplotlib-output-into-latex.html

# #contour plot
# def figsize(scale):
#     fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
#     inches_per_pt = 1.0/72.27                       # Convert pt to inch
#     golden_mean = (np.sqrt(5.0)-1)/2.0            # Aesthetic ratio (you could change this)
#     fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
#     fig_height = fig_width*golden_mean              # height in inches
#     fig_size = [fig_width,fig_height]
#     return fig_size
#
# pgf_with_latex = {                      # setup matplotlib to use latex for output
#     "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
#     "text.usetex": True,                # use LaTeX to write all text
#     "font.family": "serif",
#     "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
#     "font.sans-serif": [],
#     "font.monospace": [],
#     "axes.labelsize": 10,               # LaTeX default is 10pt font.
#     "text.fontsize": 10,
#     "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
#     "xtick.labelsize": 8,
#     "ytick.labelsize": 8,
#     "figure.figsize": figsize(1),     # default fig size of 0.9 textwidth
#     "pgf.preamble": [
#         r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
#         r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
#         ]
#     }
# mpl.rcParams.update(pgf_with_latex)

#End of plotting stuff


#in the function header we need to have the timestep and duration of simulation specified

def TB_solver():
    """Split operator technique for propogation of wave packets with square
        lattice tight-binding hamiltonian.

    """

    #need to add some assert statements ...

    #probably want to have a file with all of the details ...
    Ns = 10#DT/dt #need to then round this off to the nearest integer ...
    #might need to do up to 1000 iterations ...
    #should be able to do this with 400 iteraions tops ...
    n = 4000
    m = 4000
    N = n*m

    l_c = 200#this is probably too short ...

    #correlation length should be 20 nm ...
    points = Crystal(m, n)
    wvf = Psi(4*l_c, math.pi/(5*l_c), math.pi/(5*l_c), m, n)#wavefunction ... just put some matrix there for now ...
    #math.pi/(5*l_c)
    #the wavevector should be like 0.07 ish ... in total ..

    #the value in the Gaussian must be much less than this!

    #Importing hamiltonain from module
    H = TBH(n,m,dt=0.1e-15,V=True)
    TH1P = H[0]
    TH1N = H[1]
    TH2P = H[2]
    TH2N = H[3]

    #these are so we can vectorise the tridiagonal matrix multiplication

    TH1N0 = TH1N[0,2:N].reshape((N-2,1))
    TH1N1 = TH1N[1,1:N-1].reshape((N-2,1))
    TH1N2 = TH1N[2,0:N-2].reshape((N-2,1))

    TH2N0 = TH2N[0,2:N].reshape((N-2,1))
    TH2N1 = TH2N[1,1:N-1].reshape((N-2,1))
    TH2N2 = TH2N[2,0:N-2].reshape((N-2,1))

    #need to specify how it is constructed so we can make the tridagonly matrices

    #Empty vectors for populating
    psi_p = np.zeros((N,1),dtype=complex)
    eta_p = np.zeros((N,1),dtype=complex)
    xi_p = np.zeros((N,1),dtype=complex)

    for i in range(Ns):

        #Matrix multiplication with wavefunction and TH1N
        psi_p[0] = TH1N[1,0]*wvf[0] + TH1N[0,1]*wvf[1]
        psi_p[-1] = TH1N[1,-1]*wvf[-1] + TH1N[2,-2]*wvf[-2]
        psi_0 = np.multiply(TH1N0, wvf[2:N])
        psi_1 = np.multiply(TH1N1, wvf[1:N-1])
        psi_2 = np.multiply(TH1N2, wvf[0:N-2])
        psi_p[1:N-1] = psi_0 + psi_1 + psi_2

        #Solving linear system of equations
        eta_c = scipy.linalg.solve_banded((1,1), TH1P, psi_p)

        #Reshape vector for next calculation
        eta_nm = eta_c.reshape((n,m))
        eta_T = eta_nm.T
        eta_r = eta_T.reshape((N,1))

        #Matrix multiplication with eta and TH2N
        eta_p[0] = TH2N[1,0]*eta_r[0] + TH2N[0,1]*eta_r[1]
        eta_p[-1] = TH2N[1,-1]*eta_r[-1] + TH2N[2,-2]*eta_r[-2]
        eta_0 = np.multiply(TH2N0, eta_r[2:N])
        eta_1 = np.multiply(TH2N1, eta_r[1:N-1])
        eta_2 = np.multiply(TH2N2, eta_r[0:N-2])
        eta_p[1:N-1] = eta_0 + eta_1 + eta_2

        #Solving linear system of equations
        xi_r = scipy.linalg.solve_banded((1,1), TH2P, eta_p)

        #Reshape vector for next calculation
        xi_nm = xi_r.reshape((m,n))
        xi_T = xi_nm.T
        xi_c = xi_T.reshape((N,1))

        #so commenting out the following shifted the wavepacket in the opposite direction, but less so ...

        #Matrix multiplication with xi and TH1N
        xi_p[0] = TH1N[1,0]*xi_c[0] + TH1N[0,1]*xi_c[1]
        xi_p[-1] = TH1N[1,-1]*xi_c[-1] + TH1N[2,-2]*xi_c[-2]
        xi_0 = np.multiply(TH1N0, xi_c[2:N])
        xi_1 = np.multiply(TH1N1, xi_c[1:N-1])
        xi_2 = np.multiply(TH1N2, xi_c[0:N-2])
        xi_p[1:N-1] = xi_0 + xi_1 + xi_2

        #Solving for wavefunction
        wvf = scipy.linalg.solve_banded((1,1), TH1P, xi_p)


    #Calculate conjugate wave function
    wvf_conj = np.conjugate(wvf)

    #Probability density function by element wise multiplication
    #Neglect imaginary part
    pd = np.multiply(wvf_conj,wvf)
    pd = np.reshape(pd,(n,m))

    #Plotting
    plt.contourf(points[0].reshape((n,m)),points[1].reshape((n,m)),pd)#,cmap='RdGy'
    plt.show()

#     #plotting
#     if display:
#
#         plt.figure()
#         plt.contourf(x1, x2, z, locator=ticker.LogLocator(),cmap='RdGy')
#         cb = plt.colorbar()
#         cb.set_label('Z')
#         plt.scatter(1,1,c='black',marker='o')
#         plt.xlabel(r'$x_{1}$', fontsize = 15)
#         plt.ylabel(r'$x_{2}$', fontsize = 15)
#         plt.locator_params('x',tight=True, nbins =8)
#         plt.locator_params('y',tight=True, nbins =8)
#         plt.title('Zachary Goodwin - visualize')
#         plt.tick_params(axis='both',which='major',labelsize =15,direction='in',right=True,top=True)
#         plt.tight_layout()
#         #plt.savefig('hw212.pdf')
#         #plt.savefig('hw212.png')
#         plt.show()

if __name__ == '__main__':
    TB_solver()
