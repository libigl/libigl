cmake_minimum_required(VERSION 2.6)
project(libigl_tutorials)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)

# Chapter 1
add_subdirectory("101_FileIO")
add_subdirectory("102_DrawMesh")
add_subdirectory("103_Events")
add_subdirectory("104_Colors")
add_subdirectory("105_Overlays")
add_subdirectory("106_Picking")

# Chapter 2
add_subdirectory("201_Normals")
add_subdirectory("202_GaussianCurvature")
add_subdirectory("203_CurvatureDirections")
add_subdirectory("204_Gradient")
add_subdirectory("205_Laplacian")

# Chapter 3
add_subdirectory("301_Slice")
add_subdirectory("302_Sort")
add_subdirectory("303_LaplaceEquation")
add_subdirectory("304_LinearEqualityConstraints")
add_subdirectory("305_QuadraticProgramming")

# Chapter 5
add_subdirectory("501_HarmonicParam")
add_subdirectory("502_LSCMParam")
add_subdirectory("503_ARAPParam")
add_subdirectory("504_NRosyDesign")
add_subdirectory("505_MIQ")
add_subdirectory("506_FrameField")

# Chapter 6
add_subdirectory("601_Serialization")
add_subdirectory("602_Matlab")
add_subdirectory("604_Triangle")
add_subdirectory("605_Tetgen")
add_subdirectory("606_AmbientOcclusion")
