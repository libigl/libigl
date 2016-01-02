cmake_minimum_required(VERSION 2.6)
project(libigl_tutorials)

message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")

option(LIBIGL_USE_STATIC_LIBRARY "Use static library" OFF)
option(LIBIGL_VIEWER_WITH_NANOGUI "Enable nanogui " OFF)

if(LIBIGL_USE_STATIC_LIBRARY)
  add_definitions(-DIGL_STATIC_LIBRARY)
endif(LIBIGL_USE_STATIC_LIBRARY)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

IF(MSVC)
  # Enable parallel compilation for Visual Studio
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj /w") 

  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )

ELSE(MSVC)
  # Libigl requires a modern C++ compiler that supports c++11
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../" )
ENDIF(MSVC)

find_package(CORK QUIET)
find_package(CGAL QUIET)
find_package(EMBREE QUIET)
find_package(LIBCOMISOH QUIET)
find_package(MATLAB QUIET)
find_package(MOSEK QUIET)
find_package(TETGEN QUIET)
find_package(TINYXML2 QUIET)
find_package(TRIANGLE QUIET)
find_package(LIM QUIET)

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
include("CMakeLists.shared")

#Compile libigl
add_definitions(-DIGL_NO_OPENGL)

if(LIBIGL_USE_STATIC_LIBRARY)
  add_subdirectory("../optional" "libigl")
endif(LIBIGL_USE_STATIC_LIBRARY)
include_directories("../include")


#Compile embree
option(ENABLE_ISPC_SUPPORT " " OFF)
option(RTCORE_TASKING_SYSTEM " " INTERNAL)

# embree static build is broken on visual studio
IF(NOT MSVC) 
  option(ENABLE_STATIC_LIB " " ON)
  set(ENABLE_STATIC_LIB ON CACHE BOOL "force static build of embree" FORCE)
ENDIF(NOT MSVC)

IF(MSVC)
  option(USE_STATIC_RUNTIME " " OFF)
  set(USE_STATIC_RUNTIME OFF CACHE BOOL "disable changes in the runtime for VS" FORCE)
ENDIF(MSVC)

option(ENABLE_TUTORIALS " " OFF)

IF(EMBREE_FOUND)
  add_subdirectory("../external/embree/" "embree")    
  include_directories("../external/embree/include")    
  list(APPEND SHARED_LIBRARIES "embree")    
ENDIF(EMBREE_FOUND)

#Compile nanogui
if(LIBIGL_VIEWER_WITH_NANOGUI)
  set(NANOGUI_BUILD_PYTHON OFF CACHE BOOL " " FORCE)
  set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
  set(NANOGUI_BUILD_SHARED OFF CACHE BOOL " " FORCE)
  add_definitions(-DIGL_VIEWER_WITH_NANOGUI)
  add_subdirectory("../external/nanogui/" "nanogui")
  include_directories("../external/nanogui/include")
  include_directories("../external/nanogui/ext/nanovg/src")
  list(APPEND SHARED_LIBRARIES "nanogui")
else()
  # still need to build GLFW
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_INSTALL OFF CACHE BOOL " " FORCE)
  add_subdirectory("${GLFW_INCLUDE_DIR}/../" "glfw")
  list(APPEND SHARED_LIBRARIES "glfw")

endif()


#Compile CoMISo (if found)
# NOTE: this cmakefile works only with the
# comiso available here: https://github.com/libigl/CoMISo
IF(LIBCOMISO_FOUND)

  list(GET LIBCOMISO_INCLUDE_DIRS 0 COMISO_ROOT)
  # message( FATAL_ERROR "${COMISO_ROOT}" )

  if(APPLE)
    find_library(accelerate_library Accelerate)
    list(APPEND SHARED_LIBRARIES "CoMISo" ${accelerate_library})
  elseif(UNIX)
    find_package(BLAS REQUIRED)
    list(APPEND SHARED_LIBRARIES "CoMISo" ${BLAS_LIBRARIES})
  endif(APPLE)

  if(MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_SCL_SECURE_NO_DEPRECATE")
    #link_directories("${COMISO_ROOT}/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/lib/")
    list(APPEND SHARED_LIBRARIES "CoMISo" "${COMISO_ROOT}/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/lib/libopenblas.dll.a.lib")
  endif(MSVC)


  include_directories("${COMISO_ROOT}/CoMISo/ext/gmm-4.2/include")
  include_directories("${COMISO_ROOT}/")
  add_subdirectory("${COMISO_ROOT}/CoMISo/" "CoMISo")

  if(MSVC)
    # Copy the dll
    add_custom_target(Copy-CoMISo-DLL        # Adds a post-build event to MyTest
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${COMISO_ROOT}/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libopenblas.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/libopenblas.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${COMISO_ROOT}/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libgcc_s_seh-1.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/libgcc_s_seh-1.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${COMISO_ROOT}/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libgfortran-3.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/libgfortran-3.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${COMISO_ROOT}/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libquadmath-0.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/libquadmath-0.dll")
  endif(MSVC)

ENDIF(LIBCOMISO_FOUND)

message(STATUS "Linking with: ${LIBIGL_LIBRARIES}")
list(APPEND SHARED_LIBRARIES ${LIBIGL_LIBRARIES})

# This is necessary to work around the static linking order dependencies on linux
if(UNIX AND NOT APPLE)
  set(SHARED_LIBRARIES ${SHARED_LIBRARIES} ${SHARED_LIBRARIES})
endif(UNIX AND NOT APPLE)

#message(FATAL_ERROR ${SHARED_LIBRARIES})

# Store location of tutorial/shared directory 
set(TUTORIAL_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/shared CACHE PATH "location of shared tutorial resources")
add_definitions("-DTUTORIAL_SHARED_PATH=\"${TUTORIAL_SHARED_PATH}\"")

# Chapter 1
add_subdirectory("101_FileIO")
add_subdirectory("102_DrawMesh")
add_subdirectory("103_Events")
add_subdirectory("104_Colors")
add_subdirectory("105_Overlays")
add_subdirectory("106_ViewerMenu")

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
add_subdirectory("306_EigenDecomposition")

# Chapter 4
add_subdirectory("401_BiharmonicDeformation")
add_subdirectory("402_PolyharmonicDeformation")
add_subdirectory("403_BoundedBiharmonicWeights")
add_subdirectory("404_DualQuaternionSkinning")
add_subdirectory("405_AsRigidAsPossible")
add_subdirectory("406_FastAutomaticSkinningTransformations")
add_subdirectory("407_BiharmonicCoordinates")

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
add_subdirectory("510_Integrable")

# Chapter 6
if(TINYXML2_FOUND)
  add_subdirectory("601_Serialization")
endif(TINYXML2_FOUND)
if(MATLAB_FOUND)
  add_subdirectory("602_Matlab")
endif(MATLAB_FOUND)
add_subdirectory("604_Triangle")
add_subdirectory("605_Tetgen")

if (EMBREE_FOUND)
  add_subdirectory("606_AmbientOcclusion")
  add_subdirectory("607_Picking")
endif (EMBREE_FOUND)

if(LIM_FOUND)
  add_subdirectory("608_LIM")
endif(LIM_FOUND)
if(CGAL_FOUND)
  add_subdirectory("609_Boolean")
  add_subdirectory("610_CSGTree")
endif()
add_subdirectory("701_Statistics")
add_subdirectory("702_WindingNumber")
add_subdirectory("703_Decimation")
add_subdirectory("704_SignedDistance")
add_subdirectory("705_MarchingCubes")
