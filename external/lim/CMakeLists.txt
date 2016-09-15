cmake_minimum_required(VERSION 2.8.12)
project(lim)

add_library(lim STATIC
    Dirichlet_LIMSolver2D.cpp
    Dirichlet_LIMSolver3D.cpp
    GreenStrain_LIMSolver2D.cpp
    GreenStrain_LIMSolver3D.cpp
    Identity_LIMSolver2D.cpp
    Identity_LIMSolver3D.cpp
    Laplacian_LIMSolver2D.cpp
    Laplacian_LIMSolver3D.cpp
    LGARAP_LIMSolver2D.cpp
    LGARAP_LIMSolver3D.cpp
    LIMSolver.cpp
    LIMSolver2D.cpp
    LIMSolver3D.cpp
    LSConformal_LIMSolver2D.cpp
    NMSolver.cpp
    Poisson_LIMSolver2D.cpp
    TetrahedronMesh.cpp
    TriangleMesh.cpp
    UniformLaplacian_LIMSolver2D.cpp
    UniformLaplacian_LIMSolver3D.cpp)

target_include_directories(lim PRIVATE
  "../nanogui/ext/eigen"
  "../../include")
