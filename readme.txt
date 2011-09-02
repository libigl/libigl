igl Library - A simple c++ mesh library
Copyright 2011 - Daniele Panozzo, Alec Jacobson, Olga Diamanti
Interactive Geometry Lab - ETH Zurich

Naming standards:

- Every function must be written in a .h file with the same name of the function
- cpp files are NOT allowed
- A function can return a value only if it is a single scalar, elsewhere
  the output parameters must be passed as references. 
- Pointers are not allowed, if you need to make optional parameters 
  you should prepare a wrapper for any possible combination of them
- If an external dependency is needed it must be clearly stated at the
  top of the file. If the dependency is header only it must be placed in the "external"
  folder
- Do not use the using namespace directive anywhere. The only exception is for
  the igl namespace
  
Allowed types:

- Eigen Matrices
- Eigen Sparse Matrices
- bool
- int
- unsigned int
- double (float is NOT allowed)
- string

