# Molecular-Liquids

Program for simulating a molecular liquid of Si, using a MonteCarlo approach. The theoretical background for the development of this program is based on the Manual provided by Jeremy Harvey and the paper published by F. H. Stillinger and T. A. Weber. [1]

[1] F. H. Stillinger and T. A. Weber, Phys.Rev. B, 1985, 31, 5626 - 5271

# Requirements
gfortran compiler

# How to compile

- For only compilation, type in your terminal in the main directory: make

- For compilation and directly running the program, type: make run

- For cleaning all object files, module file and program file, type: make clean

# How to run

If the user did not type make run, the executable file is called: program.out. Then for running the program, type: ./program.out

# Input and output files

The input and output files are separated into their own directory, with the same name.

  ## - Input directory.
        
### The program needs the following files:
          
· parameters.inp: In this file all parameters needed for the simulation of the program must be specified by the user.
              
· init-geom.xyz: This file contains the initial structure of the Si molecular liquid in an XYZ file format.
              
   ## - Output directory.
   
### The program will generate the following files:
          
· prod_traj.xyz: This file contains the positions of the atoms in each sweep of the production phase in XYZ file format.
              
· gfunction.dat: This file shows two columns. The first one is the considered distances for the calculation of the RDF and the second one its corresponding values.
              
· SWM-data.out: In this file the user can find useful information about the input parameters, some steps of the simulation in the production phase, the proportion of accepted trial moves and the CPU time for the simulation.

# Program files

### The program consists on the following files:

· main.f90
· module.f90

Additionally, there is a Makefile.

  
