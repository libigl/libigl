// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

LIM is a small c++ library to compute local injective mappings of triangle and tet-meshes minimizing an arbitrary deformation energy subject to some linear positional constraints.

The easiest way is to use the functions provided in LIMSolverInterface.h. Function ComputeLIM() takes all the arguments needed and returns the resulting mesh. If you want to visualize the intermediate iterations use the functions InitLIM() and ComputeLIM_Step().
If you want to add your own deformation energy check the implementation of the existing energies which all inherit from the class LIMSolver2D or LIMSolver3D.
