cmake_minimum_required(VERSION 2.6)
project(libigl)

SET(CMAKE_SKIP_RULE_DEPENDENCY 1)

# SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../lib)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/tutorial/cmake)
find_package(EIGEN REQUIRED)

add_definitions(-DIGL_STATIC_LIBRARY)

## Check for GLFW
find_package(GLFW QUIET)
if (GLFW_FOUND)
  include_directories( ${GLFW_INCLUDE_DIR})
endif(GLFW_FOUND)

find_package(OpenGL QUIET)

## Check for Anttweakbar
find_package(ANTTWEAKBAR QUIET)
if (ANTTWEAKBAR_FOUND)
  include_directories( ${ANT_TWEAK_BAR_INCLUDE_DIR})
endif(ANTTWEAKBAR_FOUND)


## Check for CoMiSo, if not available skip the examples that depends on it
find_package(LIBCOMISO QUIET)
if (LIBCOMISO_FOUND)
  include_directories( ${LIBCOMISO_INCLUDE_DIRS})
endif(LIBCOMISO_FOUND)

## Check for MATLAB, if not available skip the examples that depends on it
find_package(MATLAB QUIET)
if (MATLAB_FOUND)
  include_directories( ${MATLAB_INCLUDE_DIR})
endif(MATLAB_FOUND)

## Check for EMBREE, if not available skip the examples that depends on it
find_package(EMBREE QUIET)
if (EMBREE_FOUND)
  include_directories( ${EMBREE_INCLUDE_DIR})
endif(EMBREE_FOUND)

## Check for CGAL, if not available skip the examples that depends on it
find_package(CGAL QUIET)

## Check for mosek
find_package(MOSEK QUIET)
if(NOT MOSEK_FOUND)
  add_definitions(-DIGL_NO_MOSEK)
endif(NOT MOSEK_FOUND)

## Check for CORK
find_package(CORK QUIET)
if (NOT CORK_FOUND)
  add_definitions(-DIGL_NO_CORK)
else(NOT CORK_FOUND)
  include_directories( ${CORK_INCLUDE_DIR})
endif(NOT CORK_FOUND)

## Check for LIM
find_package(LIM QUIET)
if(LIM_FOUND)
  include_directories( ${LIM_INCLUDE_DIR})
endif(LIM_FOUND)

## Check for SVD3X3
find_package(SVD3X3 QUIET)
if(SVD3X3_FOUND)
  include_directories( ${SVD3X3_INCLUDE_DIR})
endif(SVD3X3_FOUND)

## Check for TETGEN
find_package(TETGEN QUIET)
if(TETGEN_FOUND)
  include_directories( ${TETGEN_INCLUDE_DIR})
endif(TETGEN_FOUND)

## Check for TRIANGLE
find_package(TRIANGLE QUIET)
if(TRIANGLE_FOUND)
  include_directories( ${TRIANGLE_INCLUDE_DIR})
endif(TRIANGLE_FOUND)

## Check for TINYXML2
find_package(TINYXML2 QUIET)
if(TINYXML2_FOUND)
  include_directories( ${TINYXML2_INCLUDE_DIR})
endif(TINYXML2_FOUND)

## Check for COMISO
find_package(COMISO QUIET)
if(COMISO_FOUND)
  include_directories( ${COMISO_INCLUDE_DIR})
endif(COMISO_FOUND)

## Use openMP if available
find_package(OpenMP)
if (OPENMP_FOUND AND NOT WIN32)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#### Libigl requires a modern C++ compiler that supports c++11
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#### Compile the core library that depends only on EIGEN ####
include_directories( ${EIGEN_INCLUDE_DIR})
include_directories( ${PROJECT_SOURCE_DIR}/include/)

file(GLOB SOURCES
  "${PROJECT_SOURCE_DIR}/include/igl/*.cpp"
)

### HACKS TO CLEAN
list(REMOVE_ITEM SOURCES ${ui} ${PROJECT_SOURCE_DIR}/include/igl/cocoa_key_to_anttweakbar_key.cpp)
list(REMOVE_ITEM SOURCES ${ui} ${PROJECT_SOURCE_DIR}/include/igl/ReAntTweakBar.cpp)

add_library(igl STATIC ${SOURCES})

#### Compile the BBW part

file(GLOB SOURCES_BBW
  "${PROJECT_SOURCE_DIR}/include/igl/bbw/*.cpp"
)

add_library(iglbbw STATIC ${SOURCES_BBW})

#### Compile the mosek part (untested)
if (MOSEK_FOUND)
  file(GLOB SOURCES_MOSEK
    "${PROJECT_SOURCE_DIR}/include/igl/mosek/*.cpp"
  )

add_library(iglmosek STATIC ${SOURCES_MOSEK})
endif (MOSEK_FOUND)

#### Compile the cgal part
if (CGAL_FOUND)
  file(GLOB SOURCES_CGAL
    "${PROJECT_SOURCE_DIR}/include/igl/cgal/*.cpp"
  )

add_library(iglcgal STATIC ${SOURCES_CGAL})

#### Compile the boolean part
if (NOT CORK_FOUND)
  add_definitions(-DIGL_NO_CORK)
endif(NOT CORK_FOUND)

file(GLOB SOURCES_BOOLEAN
  "${PROJECT_SOURCE_DIR}/include/igl/boolean/*.cpp"
)

add_library(iglboolean STATIC ${SOURCES_BOOLEAN})
endif (CGAL_FOUND)

#### Compile the embree part
if (EMBREE_FOUND)
  file(GLOB SOURCES_EMBREE
    "${PROJECT_SOURCE_DIR}/include/igl/embree/*.cpp"
  )

  add_library(iglembree STATIC ${SOURCES_EMBREE})
endif (EMBREE_FOUND)

#### Compile the lim part
if (LIM_FOUND)
  file(GLOB SOURCES_LIM
    "${PROJECT_SOURCE_DIR}/include/igl/lim/*.cpp"
  )

  add_library(igllim STATIC ${SOURCES_LIM})
endif (LIM_FOUND)

#### Compile the matlab part
if (MATLAB_FOUND)
  file(GLOB SOURCES_MATLAB
    "${PROJECT_SOURCE_DIR}/include/igl/matlab/*.cpp"
  )

  add_library(iglmatlab STATIC ${SOURCES_MATLAB})
endif (MATLAB_FOUND)

#### Compile the svd3x3 part
if (SVD3X3_FOUND)
  file(GLOB SOURCES_SVD3X3
    "${PROJECT_SOURCE_DIR}/include/igl/svd3x3/*.cpp"
  )

  add_library(iglsvd3x3 STATIC ${SOURCES_SVD3X3})
endif (SVD3X3_FOUND)

#### Compile the tetgen part
if (TETGEN_FOUND)
  file(GLOB SOURCES_TETGEN
    "${PROJECT_SOURCE_DIR}/include/igl/tetgen/*.cpp"
  )

  add_library(igltetgen STATIC ${SOURCES_TETGEN})
endif (TETGEN_FOUND)

#### Compile the triangle part
if (TRIANGLE_FOUND)
  file(GLOB SOURCES_TRIANGLE
    "${PROJECT_SOURCE_DIR}/include/igl/triangle/*.cpp"
  )

  add_library(igltriangle STATIC ${SOURCES_TRIANGLE})
endif (TRIANGLE_FOUND)

#### Compile the xml part
if (TINYXML2_FOUND)
  file(GLOB SOURCES_XML
    "${PROJECT_SOURCE_DIR}/include/igl/xml/*.cpp"
  )

  add_library(iglxml STATIC ${SOURCES_XML})
endif (TINYXML2_FOUND)

#### Compile the xml part
if (LIBCOMISO_FOUND)
  file(GLOB SOURCES_COMISO
    "${PROJECT_SOURCE_DIR}/include/igl/comiso/*.cpp"
  )

  add_library(iglcomiso STATIC ${SOURCES_COMISO})
endif (LIBCOMISO_FOUND)

#### Compile the viewer
if (GLFW_FOUND AND ANTTWEAKBAR_FOUND)
  file(GLOB SOURCES_XML
    "${PROJECT_SOURCE_DIR}/include/igl/viewer/*.cpp"
  )

  add_library(iglviewer STATIC ${SOURCES_XML})
endif (GLFW_FOUND AND ANTTWEAKBAR_FOUND)
