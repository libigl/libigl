cmake_minimum_required(VERSION 2.8.12)
project(levenberg-marquadt)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(LIBIGL QUIET)
find_package(LIBHEDRA QUIET)

if (NOT LIBIGL_FOUND)
   message(FATAL_ERROR "libigl not found --- You can download it using: \n git clone --recursive https://github.com/libigl/libigl.git ${PROJECT_SOURCE_DIR}/../libigl")
endif()

if (NOT LIBHEDRA_FOUND)
   message(FATAL_ERROR "libhedra not found --- You can download it in https://github.com/avaxman/libhedra.git")
endif()

# Compilation flags: adapt to your needs 
if(MSVC)
  # Enable parallel compilation
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") 
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  # Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()

# libigl options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)

# add a customizable menu bar
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   OFF)

# libigl options: choose your dependencies (by default everything is OFF except opengl) 
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           ON)
option(LIBIGL_WITH_BBW         "Use BBW"            OFF)
option(LIBIGL_WITH_EMBREE      "Use Embree"         OFF)
option(LIBIGL_WITH_PNG         "Use PNG"            OFF)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       OFF)
option(LIBIGL_WITH_XML         "Use XML"            OFF)
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF) # This option is not supported yet
option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF) # This option is not supported yet
option(LIBIGL_WITH_CGAL        "Use CGAL"           OFF)
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, the cgal build system fails without it
  find_package(CGAL REQUIRED)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
  include(${CGAL_USE_FILE})
endif()

# Adding libigl: choose the path to your local copy libigl 
# This is going to compile everything you requested 
#message(FATAL_ERROR "${PROJECT_SOURCE_DIR}/../libigl/cmake")
add_subdirectory("${LIBIGL_INCLUDE_DIR}/../shared/cmake" "libigl")

# libigl information 
message("libigl includes: ${LIBIGL_INCLUDE_DIRS}")
message("libigl libraries: ${LIBIGL_LIBRARIES}")
message("libigl extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libigl extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libigl definitions: ${LIBIGL_DEFINITIONS}")

message("libhedra includes: ${LIBHEDRA_INCLUDE_DIRS}")

# Prepare the build environment
include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})

include_directories(${LIBHEDRA_INCLUDE_DIRS})

# Store location of data directory
set(DATA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../data CACHE PATH "location of mesh data")
add_definitions("-DDATA_PATH=\"${DATA_PATH}\"")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add your project files
FILE(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})
