# Tight Binding Solver
#### Version 1.0  -   13/02/2018


A Tight binding solver to study the evolution of Gaussian wavepackets in 2D lattices.
Uses a split Hamiltonian approach to produce fast results and is capable of creating
.mp4 videos of the propagation.

## 1. Setup

This script is written in python 3 and uses pre-installed modules
(numpy, matplotlib, scipy) and the cv2 module, installation instructions below. 

All modules for the **crystal structure**, **initial wavepacket**,
**Tight Binding Hamiltonian**, and **solvers** have been included for
**square (S)** and **hexagonal (G)** lattices. Additionally, an **animate.py**
file in included to produce .mp4 files and, most importantly, the **DD_Main.py**
file, which will be the only file you need to modify to simulate different
systems.

Once you have downlaoded these, install cv2 using the command below which will be
necessary for producing the video files.

```ruby
pip install opencv-python
```

You will also need to create an empty folder called Images into your work directory.

Finally, to make your code executable, run the following command:
```ruby
chmod +x DD_Main.py
```

## 2. Split Operator Method

The method used to solve for the time evolution of a wavepacket in this program
is the Split Operator Method. In this approach, the system is described at each
lattice site of the crystal to create a 'mesh' of atoms, as shown in the image
below, taken from [Wave packet dynamics and valley filter in strained graphene](https://arxiv.org/abs/1105.1125v1).

![Alt Text](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/GrapheneGrid.png)

The Tight Binding Hamiltonian is then constructed as a sparse matrix of on-site
and hopping elements between these lattice sites. The split operator technique
separates the contributions from lateral and vertical hopping into two
separate hamiltonians, which can both be represented as diagonal matrices and
hence leads to a considerable reduction of the computation time.

These matrices are then solved following the approach described in section 4.

## 3. System Modules
If you wish to use your own input modules for different structures, they must have
the following formats:

### Crystal
Describes the structure of your 2D lattice. This module defines the set of
coordinates in real space which are occupied by atoms. The coordinates should move
down each column, then to the top of the next, until all lattice points have been
described.

*Inputs*

* m - number of atoms along the x direction.
* n - number of atoms along the y direction.

*Returns*

* X   - complete list of atomic x-positions.
* Y   - complete list of atomic y-positions.
* X_0 - initial x-position of the centre of the wavepacket.
* Y_0 - initial y-position of the centre of the wavepacket.
    *N.B - X_0 and Y_0 __must__ be atomic positions*

##### Example: Square Lattice; m=4, n=4
o  o  o  o  
o  o  o  o  
o  o  o  o  
**o**  o  o  o    

X  = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]  
Y  = [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]  
X0 = 0  
Y0 = 3  

### Psi
Describes the initial wavepacket to be propagated through the system. This
module determines the value of a gaussian centred at (X_0, Y_0) at each atomic
position on the lattice, again moving down each column from left to right.

*Inputs*

* s     - width of gaussian Wavepacket
* kx,ky - wavenumbers along x and y directions
* m     - number of atoms along the x direction
* n     - number of atoms along the y direction
* pos   - output from the Crystal function

*Returns*

* Psi   - an (n*m)x1 matrix with the value of the wavefunction defined at each atom
          elements moving from top left to bottom right

##### Example : Square Lattice; m=4, n=4
o  o  o  o    
o  o  o  o  
o  o  o  o  
**o**  o  o  o     

Psi = [[0,0,1,5,  0,0,0,1,  0,0,0,0,  0,0,0,0]]

### Potential
Desribes the external potential that the wavepacket is within. This module can
be defined as a **1D** or **2D** potential, where the value of V is determined
at each x-position **or** at each atomic location respectively.

*Inputs*

* m      - number of atoms along the x direction
* n      - number of atoms along the y direction
* pos    - output from the Crystal function
* lc     - (optional), the disorder length

*Returns*

* Wfinal - value of V at each: x-position
                  **OR**   : atomic location

### Hamiltonian
Creates the tight binding Hamiltonian as one or several matrices to be used
to solve for the wavepacket propagation through the lattice. Included are two
modules, one for a 'Full Hamiltonian' to be used in conjunction with a simple
solver, and another which separates the Hamiltonian to be used with the Split
Operator Approach.

#### Full Hamiltonian
Describes the complete tight binding hamiltonian as a sparse matrix. Full
details can be found in the doc strings of [DD_FH_G](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/DD_TB_S.py) and [DD_FH_S](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/DD_TB_S.py).

*Inputs*

* DP - disorder potential
* n  - number of rows of carbon atoms
* m  - number of columns of carbon atoms
* dt - time step in seconds
* V  - determines what type of external potential should be included

*Returns*

* H  - sparse tight binding matrix hamiltonian

#### Split Operator Hamiltonian
Describes the complete tight binding hamiltonian as two tri-diagonal matrices.
Full details can be found in the doc strings of [DD_SH](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/DD_TB_S.py) and [DD_GH](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/DD_TB_S.py).

*Inputs*
* DP - disorder potential
* n  - number of rows of carbon atoms
* m  - number of columns of carbon atoms
* dt - time step in seconds


*Returns*
* TH1P - 3x(n*m) tridiagonal matrix for H_m of split operator technique.
        Positive version in the linear equation.
* TH1N - 3x(n*m) tridiagonal matrix for H_m of split operator technique.
        Negative version in the linear equation.
* TH2P - 3x(n*m) tridiagonal matrix for H_n of split operator technique.
        Positive version in the linear equation.
* TH2N - 3x(n*m) tridiagonal matrix for H_n of split operator technique.
        Negative version in the linear equation.

## 4. Solvers
The modules to solve for the time propagation are described in [DD_SS](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/DD_SS.py) for the
simple solver, and in [DD_TB_S](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/DD_TB_S.py) for the split operator solver. As the time
evolution of the wavepacket is contained within an exponential term involving
the Hamiltonian, it must be expanded in the Cayley form to solve as an
eigenvalue problem.

### Simple solver
This module takes as an input the single, full hamiltonian and takes the dot
product with the wavefunction to find the updated wavefunction at each time step.
More details are given in the doc string for [DD_SS](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/DD_SS.py)

### Split Operator Solver
This module takes all 4 matrices created in the Split Operator Hamiltonian
modules. The method for solving this at each time step is outlined in
[Wave packet dynamics and valley filter in strained graphene](https://arxiv.org/abs/1105.1125v1). In essence, writing the full
hamiltonian as a sum allows the exponential to be factorised into 3 terms; the
first and last depend on the vertical hopping, while the central term depends on
the lateral hopping. Each exponential can then be written in Cayley form and
solved as an eigenvalue problem. This requires a certain amount of restructuring
of the matrices which is described in more detail in the doc string of [DD_TB_S](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/DD_TB_S.py).

This equation is solved for each time step, updating the input wavefunction each
time until the final form is obtained.

## 5. Tests
To test the program is working correctly, complete the following steps:

1. Follow installation instructions provided in section 1.
2. Open a unix terminal and run the command
  ```
  chmod +x DD_test.py
  ```
  This will make the program executable.
3. In your terminal, run the command
  ```
  ./DD_test.py
  ```
  This will run a test code to produce a gaussian packet on a square lattice. If working correctly 
  two lists will be output to the terminal and an image should appear. The first list is the coordinates and
  corresponding values of the hamiltonian, the second list is the vector representing the disorder potential,
  and the image will display a 2D-Gaussian on a square grid, as shown below: 
  
  ![alt text](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/test.png "Initial guassian on a square lattice")
  
  

## 6. Credits
This program was written as a part of the Computational Methods course for the
Centre for Doctoral Training in Theory and Simulation of Materials at Imperial
College London.

It was written by:
* Zachery Goodwin (https://gitlab.com/zachary.goodwin13)
* Georgios Samaras (https://gitlab.com/Samaras)
* Rosemary Teague (https://gitlab.com/RoseTeague)
* Amy Wang (https://gitlab.com/yiyuan.wang16)

## 7. License
GNU LESSER General Public License. Copyright 2007.

Under this [license](https://gitlab.com/RoseTeague/Group_Programming_Project/blob/master/LICENSE), everyone is permitted to copy and distribute verbatim copies
of this document, but changing it is not allowed.
