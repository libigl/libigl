cmake_minimum_required(VERSION 2.8.12)
project(libigl_tutorials)
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   OFF)

### libIGL options: choose your dependencies (by default everything is OFF, in this example we need the viewer) ###
find_package(CGAL QUIET COMPONENTS Core)
option(LIBIGL_WITH_CGAL             "Use CGAL"           "${CGAL_FOUND}")
option(LIBIGL_WITH_COMISO           "Use CoMiso"         ON)
option(LIBIGL_WITH_CORK             "Use CORK"           OFF)
option(LIBIGL_WITH_EMBREE           "Use Embree"         ON)
option(LIBIGL_WITH_LIM              "Use LIM"            ON)
find_package(MATLAB QUIET)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         "${MATLAB_FOUND}")
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          "${MOSEK_FOUND}")
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW      "Use OpenGL GLFW"    ON)
option(LIBIGL_WITH_PNG              "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       ON)
option(LIBIGL_WITH_XML              "Use XML"            ON)
### End   to be tested ----

### libIGL options: decide if you want to use the functionalities that depends on cgal
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, cgal strange build system fails without it
  find_package(CGAL REQUIRED COMPONENTS Core)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
  include(${CGAL_USE_FILE})
endif()

### Adding libIGL: choose the path to your local copy libIGL ###
### This is going to compile everything you requested ###
add_subdirectory("${PROJECT_SOURCE_DIR}/../shared/cmake" "libigl")


### Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

### Compilation flags: adapt to your needs ###
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") ### Enable parallel compilation
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /w") # disable all warnings (not ideal but...)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w") # disable all warnings (not ideal but...)
endif()

# Enable openMP if possible
#find_package(OpenMP)
#if (OPENMP_FOUND AND NOT WIN32)
#  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()


### Prepare the build environment

include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})

### Choose which chapters to compile ###
option(TUTORIALS_CHAPTER1 "Compile chapter 1" ON)
option(TUTORIALS_CHAPTER2 "Compile chapter 2" ON)
option(TUTORIALS_CHAPTER3 "Compile chapter 3" ON)
option(TUTORIALS_CHAPTER4 "Compile chapter 4" ON)
option(TUTORIALS_CHAPTER5 "Compile chapter 5" ON)
option(TUTORIALS_CHAPTER6 "Compile chapter 6" ON)
option(TUTORIALS_CHAPTER7 "Compile chapter 7" ON)

# Store location of tutorial/shared directory
set(TUTORIAL_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/shared CACHE PATH "location of shared tutorial resources")
add_definitions("-DTUTORIAL_SHARED_PATH=\"${TUTORIAL_SHARED_PATH}\"")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Chapter 1
if(TUTORIALS_CHAPTER1)
  add_subdirectory("101_FileIO")
  add_subdirectory("102_DrawMesh")
  add_subdirectory("103_Events")
  add_subdirectory("104_Colors")
  add_subdirectory("105_Overlays")
  add_subdirectory("106_ViewerMenu")
endif()

# Chapter 2
if(TUTORIALS_CHAPTER2)
  add_subdirectory("201_Normals")
  add_subdirectory("202_GaussianCurvature")
  add_subdirectory("203_CurvatureDirections")
  add_subdirectory("204_Gradient")
  add_subdirectory("205_Laplacian")
endif()

# Chapter 3
if(TUTORIALS_CHAPTER3)
  add_subdirectory("301_Slice")
  add_subdirectory("302_Sort")
  add_subdirectory("303_LaplaceEquation")
  add_subdirectory("304_LinearEqualityConstraints")
  add_subdirectory("305_QuadraticProgramming")
  add_subdirectory("306_EigenDecomposition")
endif()

# Chapter 4
if(TUTORIALS_CHAPTER4)
  add_subdirectory("401_BiharmonicDeformation")
  add_subdirectory("402_PolyharmonicDeformation")
  add_subdirectory("403_BoundedBiharmonicWeights")
  add_subdirectory("404_DualQuaternionSkinning")
  add_subdirectory("405_AsRigidAsPossible")
  add_subdirectory("406_FastAutomaticSkinningTransformations")
  add_subdirectory("407_BiharmonicCoordinates")
endif()

# Chapter 5
if(TUTORIALS_CHAPTER5)
  add_subdirectory("501_HarmonicParam")
  add_subdirectory("502_LSCMParam")
  add_subdirectory("503_ARAPParam")
  if(LIBIGL_WITH_COMISO)
    add_subdirectory("504_NRosyDesign")
    add_subdirectory("505_MIQ")
    add_subdirectory("506_FrameField")
  endif()
  add_subdirectory("507_PolyVectorField")
  add_subdirectory("508_ConjugateField")
  add_subdirectory("509_Planarization")
  add_subdirectory("510_Integrable")
  add_subdirectory("511_PolyVectorFieldGeneral")
endif()

# Chapter 6
if(TUTORIALS_CHAPTER6)
  if(LIBIGL_WITH_XML)
    add_subdirectory("601_Serialization")
  endif()
  if(LIBIGL_WITH_MATLAB)
    add_subdirectory("602_Matlab")
  endif()
  if(LIBIGL_WITH_TRIANGLE)
    add_subdirectory("604_Triangle")
  endif()
  if(LIBIGL_WITH_TETGEN)
    add_subdirectory("605_Tetgen")
  endif()
  if(LIBIGL_WITH_EMBREE)
    add_subdirectory("606_AmbientOcclusion")
  endif()
  if(LIBIGL_WITH_PNG)
    add_subdirectory("607_ScreenCapture")
  endif()
  if(LIBIGL_WITH_LIM)
    add_subdirectory("608_LIM")
  endif()
  if(LIBIGL_WITH_CGAL)
    add_subdirectory("609_Boolean")
    add_subdirectory("610_CSGTree")
  endif()
endif()

# Chapter 7
if(TUTORIALS_CHAPTER7)
  add_subdirectory("701_Statistics")
  add_subdirectory("702_WindingNumber")
  add_subdirectory("703_Decimation")
  add_subdirectory("704_SignedDistance")
  add_subdirectory("705_MarchingCubes")
  if(LIBIGL_WITH_EMBREE)
    add_subdirectory("706_FacetOrientation")
  endif()
  add_subdirectory("707_SweptVolume")
  add_subdirectory("708_Picking")
  if(LIBIGL_WITH_COMISO)
    add_subdirectory("709_VectorFieldVisualizer")
  endif()
  add_subdirectory("710_SLIM")
  add_subdirectory("711_Subdivision")
endif()
