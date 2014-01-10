
    --------------------- Open Numerical Library 3.2 -------------------------

General information
===================


This is OpenNL, a library to easily construct and solve sparse linear systems.
* OpenNL is supplied with a set of built-in iterative solvers (Conjugate
  gradient, BICGSTAB, GMRes) and preconditioners (Jacobi, SSOR). 
* OpenNL can also use other solvers (SuperLU 4.0 is supported as an OpenNL
  extension, CNC can be enabled on compatible hardware, MUMPS will be supported
  in a future version)
  
To install OpenNL please look at the INSTALL.txt 

The OpenNL API is described in User_Guide.txt

    See 
    http://alice.loria.fr/software/OpenNL/
    for information and additional example data sets.



Changelog:

3.2 (Sept. 2010) (fix)
- Added some missing documentation

3.1 (Sept. 2010)
- CNC: support for CUDA 2.0 capabilities, compatibility with Fermi-specific features

3.0 (Feb. 2010)
- added CNC, with CRS, BCRS2, ELL, Hybrid(ELL+COO) kernels, in simple and double
    precision versions
- uses CMake


Contributors:

Bruno LEVY
Luc BUATOIS
Marc FUENTES
Thomas JOST
Nicolas SAUGNIER


