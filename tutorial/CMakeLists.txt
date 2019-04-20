cmake_minimum_required(VERSION 3.1)
project(libigl_tutorials)

### Conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

### Compile libigl in header-only mode for Python bindings
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)

### Adding libIGL: choose the path to your local copy libIGL
if(NOT TARGET igl_common)
  include(libigl)
else()
  include(LibiglDownloadExternal)
endif()

### Download data
igl_download_tutorial_data()

### Output directories
if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")
endif()

### Choose which chapters to compile
option(TUTORIALS_CHAPTER1 "Compile chapter 1" ON)
option(TUTORIALS_CHAPTER2 "Compile chapter 2" ON)
option(TUTORIALS_CHAPTER3 "Compile chapter 3" ON)
option(TUTORIALS_CHAPTER4 "Compile chapter 4" ON)
option(TUTORIALS_CHAPTER5 "Compile chapter 5" ON)
option(TUTORIALS_CHAPTER6 "Compile chapter 6" ON)
option(TUTORIALS_CHAPTER7 "Compile chapter 7" ON)

# Store location of tutorial/shared directory
set(TUTORIAL_SHARED_PATH ${LIBIGL_EXTERNAL}/../tutorial/data CACHE PATH "location of shared tutorial resources")
add_library(tutorials INTERFACE)
target_compile_definitions(tutorials INTERFACE "-DTUTORIAL_SHARED_PATH=\"${TUTORIAL_SHARED_PATH}\"")
target_include_directories(tutorials INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Chapter 1
if(TUTORIALS_CHAPTER1)
  add_subdirectory("101_FileIO")
endif()

#######################
if(NOT (LIBIGL_WITH_OPENGL AND LIBIGL_WITH_OPENGL_GLFW) )
  message(WARNING "Most tutorial executables depend on OpenGL and glfw. Use `cmake ../ -DLIBIGL_WITH_OPENGL=ON -DLIBIGL_WITH_OPENGL_GLFW=ON`")
  return()
endif()

if(TUTORIALS_CHAPTER1)
  add_subdirectory("102_DrawMesh")
  add_subdirectory("103_Events")
  add_subdirectory("104_Colors")
  if(LIBIGL_WITH_OPENGL_GLFW_IMGUI)
    add_subdirectory("105_Overlays")
    add_subdirectory("106_ViewerMenu")
  endif()
  add_subdirectory("107_MultipleMeshes")
  add_subdirectory("108_MultipleViews")
endif()

# Chapter 2
if(TUTORIALS_CHAPTER2)
  add_subdirectory("201_Normals")
  add_subdirectory("202_GaussianCurvature")
  add_subdirectory("203_CurvatureDirections")
  add_subdirectory("204_Gradient")
  add_subdirectory("205_Laplacian")
  add_subdirectory("206_GeodesicDistance")
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
  add_subdirectory("507_Planarization")
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
  add_subdirectory("709_SLIM")
  add_subdirectory("710_SCAF")
  add_subdirectory("711_Subdivision")
  add_subdirectory("712_DataSmoothing")
  add_subdirectory("713_ShapeUp")
  if(LIBIGL_WITH_TETGEN)
    add_subdirectory("714_MarchingTets")
  endif()
  add_subdirectory("715_MeshImplicitFunction")
  add_subdirectory("716_HeatGeodesics")
endif()
