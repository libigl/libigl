cmake_minimum_required(VERSION 2.6)
project(libigl_tutorials)

message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")

option(LIBIGL_USE_STATIC_LIBRARY "Use static library" OFF)

if(LIBIGL_USE_STATIC_LIBRARY)
  add_definitions(-DIGL_STATIC_LIBRARY)
endif(LIBIGL_USE_STATIC_LIBRARY)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

IF(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") ### Enable parallel compilation for Visual Studio
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
ELSE(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../" )
ENDIF(MSVC)

find_package(CORK QUIET)
find_package(CGAL QUIET)
find_package(EMBREE QUIET)
find_package(LIBCOMISO QUIET)
find_package(MATLAB QUIET)
find_package(MOSEK QUIET)
find_package(TETGEN QUIET)
find_package(TINYXML2 QUIET)
find_package(TRIANGLE QUIET)

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
include("CMakeLists.shared")

#Compile libigl
add_definitions(-DIGL_NO_OPENGL)

if(LIBIGL_USE_STATIC_LIBRARY)
  add_subdirectory("../optional" "libigl")
endif(LIBIGL_USE_STATIC_LIBRARY)
include_directories("../include")

#Compile nanogui
add_subdirectory("../external/nanogui/" "nanogui")
include_directories("../external/nanogui/include")
include_directories("../external/nanogui/ext/nanovg/src")
list(APPEND SHARED_LIBRARIES "nanogui" "glfw")

message(STATUS "Linking with: ${LIBIGL_LIBRARIES}")
list(APPEND SHARED_LIBRARIES ${LIBIGL_LIBRARIES})

# This is necessary to work around the static linking order dependencies on linux
if(UNIX AND NOT APPLE)
  set(SHARED_LIBRARIES ${SHARED_LIBRARIES} ${SHARED_LIBRARIES})
endif(UNIX AND NOT APPLE)

#message(FATAL_ERROR ${SHARED_LIBRARIES})

# Chapter 1
add_subdirectory("101_FileIO")
add_subdirectory("102_DrawMesh")
add_subdirectory("103_Events")
add_subdirectory("104_Colors")
add_subdirectory("105_Overlays")

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

# Chapter 4
add_subdirectory("401_BiharmonicDeformation")
add_subdirectory("402_PolyharmonicDeformation")
add_subdirectory("403_BoundedBiharmonicWeights")
add_subdirectory("404_DualQuaternionSkinning")
add_subdirectory("405_AsRigidAsPossible")
add_subdirectory("406_FastAutomaticSkinningTransformations")

# Chapter 5
add_subdirectory("501_HarmonicParam")
add_subdirectory("502_LSCMParam")
add_subdirectory("503_ARAPParam")
if(LIBCOMISO_FOUND)
  add_subdirectory("504_NRosyDesign")
  add_subdirectory("505_MIQ")
  add_subdirectory("506_FrameField")
endif(LIBCOMISO_FOUND)
add_subdirectory("507_PolyVectorField")
add_subdirectory("508_ConjugateField")
add_subdirectory("509_Planarization")

# Chapter 6
if(TINYXML2_FOUND)
  add_subdirectory("601_Serialization")
endif(TINYXML2_FOUND)
if(MATLAB_FOUND)
  add_subdirectory("602_Matlab")
endif(MATLAB_FOUND)
add_subdirectory("604_Triangle")
add_subdirectory("605_Tetgen")
if(EMBREE_FOUND)
  add_subdirectory("606_AmbientOcclusion")
  add_subdirectory("607_Picking")
endif(EMBREE_FOUND)
if(LIM_FOUND)
  add_subdirectory("608_LIM")
endif(LIM_FOUND)
if(CGAL_FOUND)
  add_subdirectory("609_Boolean")
endif()
add_subdirectory("701_Statistics")
add_subdirectory("702_WindingNumber")
add_subdirectory("703_Decimation")
add_subdirectory("704_SignedDistance")
