cmake_minimum_required(VERSION 3.1)
project(libigl)

# Detects whether this is a top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(LIBIGL_TOPLEVEL_PROJECT ON)
else()
	set(LIBIGL_TOPLEVEL_PROJECT OFF)
endif()

# Build tests, tutorials and python bindings
option(LIBIGL_BUILD_TESTS      "Build libigl unit test"        ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_BUILD_TUTORIALS  "Build libigl tutorial"         ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_BUILD_PYTHON     "Build libigl python bindings"  ${LIBIGL_TOPLEVEL_PROJECT})

# USE_STATIC_LIBRARY speeds up the generation of multiple binaries,
# at the cost of a longer initial compilation time
# (by default, static build is off since libigl is a header-only library)
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)

# All dependencies that are downloaded as cmake projects and tested on the auto-builds are ON
# (by default, all build options are off)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   ON)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   ON)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    ON)
option(LIBIGL_WITH_PNG               "Use PNG"                      ON)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
option(LIBIGL_WITH_XML               "Use XML"                      ON)
option(LIBIGL_WITH_PYTHON            "Use Python"                   ${LIBIGL_BUILD_PYTHON})
### End

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
	enable_testing()
	add_subdirectory(tests)
endif()

if(LIBIGL_WITH_PYTHON)
	add_subdirectory(python)
endif()
