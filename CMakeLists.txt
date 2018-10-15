cmake_minimum_required(VERSION 3.1)
project(libigl)


option(LIBIGL_BUILD_TESTS      "Build libigl unit test"        OFF)
option(LIBIGL_BUILD_TUTORIALS  "Build libigl tutorial"         OFF)
option(LIBIGL_BUILD_PYTHON     "Build libigl python bindings"  OFF)


### Available options ###
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ON)
option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         ON)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ${LIBIGL_BUILD_TUTORIALS})
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ${LIBIGL_BUILD_TUTORIALS})
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ${LIBIGL_BUILD_TUTORIALS})
option(LIBIGL_WITH_PNG               "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ${LIBIGL_BUILD_TUTORIALS})
option(LIBIGL_WITH_XML               "Use XML"            ON)
option(LIBIGL_WITH_PYTHON            "Use Python"         ${LIBIGL_BUILD_PYTHON})

# These ensure that lib*.a are placed in the directory where `cmake
# ../optional/` was issued.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

### Adding libIGL: choose the path to your local copy libIGL
include(libigl)


if(LIBIGL_BUILD_TUTORIALS)
	add_subdirectory(tutorial)
endif()

if(LIBIGL_BUILD_TESTS)
	add_subdirectory(tests)
endif()
