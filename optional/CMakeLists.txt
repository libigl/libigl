cmake_minimum_required(VERSION 2.6)
project(libigl)

SET(CMAKE_SKIP_RULE_DEPENDENCY 1)

# SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../lib)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../tutorial/cmake)
find_package(EIGEN REQUIRED)

add_definitions(-DIGL_STATIC_LIBRARY)
option(LIBIGL_VIEWER_WITH_NANOGUI "Enable nanogui " OFF)

#option(LIBIGL_NO_OPENGL "Disable OpenGL2 legacy features in main library" ON)
#option(LIBIGL_OPENGL_4 "Enable OpenGL4 features in optional libraries" ON)
#
#if(LIBIGL_NO_OPENGL)
#  add_definitions(-DIGL_NO_OPENGL)
#endif(LIBIGL_NO_OPENGL)
#if(LIBIGL_OPENGL_4)
#  add_definitions(-DIGL_OPENGL_4)
#endif(LIBIGL_OPENGL_4)

## Check for GLFW
find_package(GLFWH QUIET)
if (GLFW_FOUND)
  include_directories( ${GLFW_INCLUDE_DIR})
else (GLFW_FOUND)
  message(STATUS "GLFW not found.")
endif(GLFW_FOUND)

## Check for GLEW
find_package(GLEW QUIET)
if (GLEW_FOUND)
  include_directories( ${GLEW_INCLUDE_DIR})
else (GLEW_FOUND)
  message(STATUS "GLEW not found.")
endif(GLEW_FOUND)

find_package(OpenGL QUIET)

## Check for Anttweakbar
find_package(ANTTWEAKBAR QUIET)
if (ANTTWEAKBAR_FOUND)
  include_directories( ${ANT_TWEAK_BAR_INCLUDE_DIR})
else (ANTTWEAKBAR_FOUND)
  message(STATUS "AntTweakBar not found.")
endif(ANTTWEAKBAR_FOUND)

## Check for NANOGUI
if(LIBIGL_VIEWER_WITH_NANOGUI)
  find_package(NANOGUIH QUIET)
  #message(FATAL_ERROR "nanogui not found. ${NANOGUI_INCLUDE_DIRS}")
  if (NANOGUI_FOUND)
    include_directories( ${NANOGUI_INCLUDE_DIRS})
  else (NANOGUI_FOUND)
    message(STATUS "nanogui not found.")
  endif(NANOGUI_FOUND)
endif()

## Check for Boost
find_package(BOOST QUIET)
if (Boost_FOUND)
  include_directories( ${BOOST_INCLUDE_DIR})
endif(Boost_FOUND)


## Check for CoMiSo, if not available skip the examples that depends on it
find_package(LIBCOMISOH QUIET)
if (LIBCOMISO_FOUND)
  include_directories( ${LIBCOMISO_INCLUDE_DIRS})
endif(LIBCOMISO_FOUND)

## Check for MATLAB, if not available skip the examples that depends on it
find_package(MATLAB QUIET)
if (MATLAB_FOUND)
  include_directories( ${MATLAB_INCLUDE_DIR})
else (MATLAB_FOUND)
  message(STATUS "MATLAB not found.")
endif(MATLAB_FOUND)

## Check for EMBREE, if not available skip the examples that depends on it
find_package(EMBREEH QUIET)
if (EMBREE_FOUND)
  include_directories( ${EMBREE_INCLUDE_DIR})
else (EMBREE_FOUND)
  message(STATUS "Embree not found.")
endif(EMBREE_FOUND)

## Check for CGAL, if not available skip the examples that depends on it
find_package(CGAL QUIET)

## Check for mosek
find_package(MOSEK QUIET)
if(MOSEK_FOUND)
  include_directories(${MOSEK_INCLUDE_DIR})
else(MOSEK_FOUND)
  add_definitions(-DIGL_NO_MOSEK)
endif(MOSEK_FOUND)

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

## Check for PNG
find_package(PNG QUIET)
if(PNG_FOUND)
  include_directories( ${PNG_INCLUDE_DIR})
else(PNG_FOUND)
  message(STATUS "PNG not found")
endif(PNG_FOUND)

## Check for TETGEN
find_package(TETGEN QUIET)
if(TETGEN_FOUND)
  include_directories( ${TETGEN_INCLUDE_DIR})
else(TETGEN_FOUND)
  message(STATUS "Tetgen not found.")
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

## Check for YIMG
find_package(YIMG QUIET)
if(YIMG_FOUND)
  include_directories( ${YIMG_INCLUDE_DIR})
else(YIMG_FOUND)
  message(STATUS "YIMG not found")
endif(YIMG_FOUND)

IF(MSVC)
### Enable parallel compilation for Visual Studio
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
ELSE(MSVC)
#### Libigl requires a modern C++ compiler that supports c++11
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF(MSVC)

#### Compile the core library that depends only on EIGEN ####
include_directories( ${EIGEN_INCLUDE_DIR})
include_directories( ${PROJECT_SOURCE_DIR}/../include/)

file(GLOB SOURCES
  "${PROJECT_SOURCE_DIR}/../include/igl/*.cpp"
  "${PROJECT_SOURCE_DIR}/../include/igl/copyleft/*.cpp"
)

add_library(igl STATIC ${SOURCES})
set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "igl")

#### Compile the AntTweakBar part
if (ANTTWEAKBAR_FOUND)
  file(GLOB SOURCES_ANTTWEAKBAR
    "${PROJECT_SOURCE_DIR}/../include/igl/anttweakbar/*.cpp"
  )
add_library(iglanttweakbar STATIC ${SOURCES_ANTTWEAKBAR})
set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglanttweakbar")
endif (ANTTWEAKBAR_FOUND)

#### Compile the BBW part
file(GLOB SOURCES_BBW
  "${PROJECT_SOURCE_DIR}/../include/igl/bbw/*.cpp"
)
add_library(iglbbw STATIC ${SOURCES_BBW})
set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglbbw")

#### Compile the boolean and cgal parts
if (CGAL_FOUND)
#### Compile the boolean part
if (NOT CORK_FOUND)
  add_definitions(-DIGL_NO_CORK)
endif(NOT CORK_FOUND)
file(GLOB SOURCES_BOOLEAN
  "${PROJECT_SOURCE_DIR}/../include/igl/copyleft/boolean/*.cpp"
)
add_library(iglboolean STATIC ${SOURCES_BOOLEAN})
set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglboolean")
#### Compile the cgal part
  file(GLOB SOURCES_CGAL
    "${PROJECT_SOURCE_DIR}/../include/igl/copyleft/cgal/*.cpp"
  )
add_library(iglcgal STATIC ${SOURCES_CGAL})
set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglcgal")
endif (CGAL_FOUND)

#### Compile the comiso part
if (LIBCOMISO_FOUND)
  file(GLOB SOURCES_COMISO
    "${PROJECT_SOURCE_DIR}/../include/igl/comiso/*.cpp"
  )
  add_library(iglcomiso STATIC ${SOURCES_COMISO})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglcomiso")
endif (LIBCOMISO_FOUND)

#### Compile the embree part
if (EMBREE_FOUND)
  file(GLOB SOURCES_EMBREE
    "${PROJECT_SOURCE_DIR}/../include/igl/embree/*.cpp"
  )
  add_library(iglembree STATIC ${SOURCES_EMBREE})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglembree")
endif (EMBREE_FOUND)

#### Compile the lim part
if (LIM_FOUND)
  file(GLOB SOURCES_LIM
    "${PROJECT_SOURCE_DIR}/../include/igl/lim/*.cpp"
  )
  add_library(igllim STATIC ${SOURCES_LIM})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "igllim")
endif (LIM_FOUND)

#### Compile the matlab part
if (MATLAB_FOUND)
  file(GLOB SOURCES_MATLAB
    "${PROJECT_SOURCE_DIR}/../include/igl/matlab/*.cpp"
  )
  add_library(iglmatlab STATIC ${SOURCES_MATLAB})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglmatlab")
endif (MATLAB_FOUND)

#### Compile the mosek part
if (MOSEK_FOUND)
  file(GLOB SOURCES_MOSEK
    "${PROJECT_SOURCE_DIR}/../include/igl/mosek/*.cpp"
  )
add_library(iglmosek STATIC ${SOURCES_MOSEK})
set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglmosek")
endif (MOSEK_FOUND)

#### Compile the opengl parts
if (OPENGL_FOUND)
  file(GLOB SOURCES_OPENGL2
    "${PROJECT_SOURCE_DIR}/../include/igl/opengl2/*.cpp"
  )
  add_library(iglopengl2 STATIC ${SOURCES_OPENGL2})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglopengl2")
  file(GLOB SOURCES_OPENGL
    "${PROJECT_SOURCE_DIR}/../include/igl/opengl/*.cpp"
  )
  add_library(iglopengl STATIC ${SOURCES_OPENGL})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglopengl")
endif (OPENGL_FOUND)

#### Compile the png part
if (PNG_FOUND AND YIMG_FOUND AND OPENGL_FOUND)
  file(GLOB SOURCES_PNG
    "${PROJECT_SOURCE_DIR}/../include/igl/png/*.cpp"
  )
  add_library(iglpng STATIC ${SOURCES_PNG})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglpng")

endif (PNG_FOUND AND YIMG_FOUND AND OPENGL_FOUND)

#### Compile the tetgen part
if (TETGEN_FOUND)
  file(GLOB SOURCES_TETGEN
    "${PROJECT_SOURCE_DIR}/../include/igl/copyleft/tetgen/*.cpp"
  )

  add_library(igltetgen STATIC ${SOURCES_TETGEN})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "igltetgen")

endif (TETGEN_FOUND)

#### Compile the triangle part
if (TRIANGLE_FOUND)
  file(GLOB SOURCES_TRIANGLE
    "${PROJECT_SOURCE_DIR}/../include/igl/triangle/*.cpp"
  )

  add_library(igltriangle STATIC ${SOURCES_TRIANGLE})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "igltriangle")

endif (TRIANGLE_FOUND)

#### Compile the viewer
if(LIBIGL_VIEWER_WITH_NANOGUI)
  if (NANOGUI_FOUND)
    add_definitions(-DIGL_VIEWER_WITH_NANOGUI)
  endif (NANOGUI_FOUND)
endif()
file(GLOB SOURCES_VIEWER
  "${PROJECT_SOURCE_DIR}/../include/igl/viewer/*.cpp"
)
add_library(iglviewer STATIC ${SOURCES_VIEWER})
set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglviewer")

#### Compile the xml part
if (TINYXML2_FOUND)
  file(GLOB SOURCES_XML
    "${PROJECT_SOURCE_DIR}/../include/igl/xml/*.cpp"
  )
  add_library(iglxml STATIC ${SOURCES_XML})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} "iglxml")
endif (TINYXML2_FOUND)

# Pass the list of compiled libraries to the parent if there is one
if(NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} PARENT_SCOPE)
endif (NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
