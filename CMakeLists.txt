# Detects whether this is a top-level project
get_directory_property(LIBIGL_PARENT_DIR PARENT_DIRECTORY)
if(NOT LIBIGL_PARENT_DIR)
    set(LIBIGL_TOPLEVEL_PROJECT ON)
else()
    set(LIBIGL_TOPLEVEL_PROJECT OFF)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.16.0")
if(LIBIGL_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build Libigl is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

################################################################################
project(libigl VERSION 3.0.0)

# Build tests and tutorials
option(LIBIGL_BUILD_TESTS      "Build libigl unit test"      ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_BUILD_TUTORIALS  "Build libigl tutorial"       ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_EXPORT_TARGETS   "Export libigl CMake targets" ${LIBIGL_TOPLEVEL_PROJECT})

# USE_STATIC_LIBRARY speeds up the generation of multiple binaries,
# at the cost of a longer initial compilation time
# (by default, static build is off since libigl is a header-only library)
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" ${LIBIGL_TOPLEVEL_PROJECT})

# Available options
option(LIBIGL_WITH_CGAL              "Use CGAL"             ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_COMISO            "Use CoMiso"           ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_CORK              "Use Cork"             ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_EMBREE            "Use Embree"           ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_MATLAB            "Use Matlab"           ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_MOSEK             "Use MOSEK"            ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_OPENGL            "Use OpenGL"           ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"             ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"            ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_PNG               "Use PNG"              ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_TETGEN            "Use Tetgen"           ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"         ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_PREDICATES        "Use exact predicates" ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_WITH_XML               "Use XML"              ${LIBIGL_TOPLEVEL_PROJECT})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

list(PREPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_LIST_DIR}/cmake/
    ${CMAKE_CURRENT_LIST_DIR}/cmake/igl
    ${CMAKE_CURRENT_LIST_DIR}/cmake/recipes
)

# Include CMake helper functions
include(igl_add_library)
include(igl_add_tutorial)
include(igl_target_sources)

# Defines CMake targets for selected libigl modules
include(libigl)

if(LIBIGL_BUILD_TUTORIALS)
    add_subdirectory(tutorial)
endif()

if(LIBIGL_BUILD_TESTS)
    # include(CTest)
    # add_subdirectory(tests)
endif()
