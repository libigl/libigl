cmake_minimum_required(VERSION 2.6)
project(libigl)

### Compilation flags: adapt to your needs ###
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") ### Enable parallel compilation
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../" )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)
if(LIBIGL_USE_STATIC_LIBRARY)
  add_definitions(-DIGL_STATIC_LIBRARY)
endif()

### libIGL options: decide if you want the nanogui bar in the viewer
option(LIBIGL_WITH_VIEWER      "Use Nanogui Viewer" OFF)
if(LIBIGL_WITH_VIEWER)
  add_definitions(-DIGL_VIEWER_WITH_NANOGUI)
endif()

### libIGL options: choose your dependencies (by default everything is OFF, in this example we need the viewer) ###
option(LIBIGL_WITH_BBW         "Use BBW"            ON)
option(LIBIGL_WITH_EMBREE      "Use Embree"         ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_PNG         "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       ON)
option(LIBIGL_WITH_XML         "Use XML"            ON)

### Still need to be tested ----
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_BOOLEAN     "Use Cork boolean"   OFF)
option(LIBIGL_WITH_CGAL        "Use CGAL"           OFF)
#option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF)
#option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF)

### END TESTING

### Adding libIGL: choose the path to your local copy libIGL ###
### This is going to compile everything you requested ###
add_subdirectory("${PROJECT_SOURCE_DIR}/../cmake" "libigl")

### ligIGL information ###
message("libIGL includes: ${LIBIGL_INCLUDE_DIRS}")
message("libIGL libraries: ${LIBIGL_LIBRARIES}")
message("libIGL extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libIGL extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libIGL definitions: ${LIBIGL_DEFINITIONS}")
