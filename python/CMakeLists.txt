cmake_minimum_required(VERSION 2.8.12)
project(pybind)

### Compilation flags: adapt to your needs ###
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj /w") ### Enable parallel compilation
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations -Wno-unused-parameter -Wno-deprecated-register -Wno-return-type-c-linkage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated-declarations -Wno-unused-parameter -Wno-deprecated-register -Wno-return-type-c-linkage")
endif()

# Force a specific python version
# SET(PYTHON_LIBRARIES "D:/Python34/libs/python34.lib")
# SET(PYTHON_INCLUDE_DIR "D:/Python34/include")

# Force a specific python version
# SET(PYTHON_LIBRARIES "/usr/local/Cellar/python3/3.5.0/Frameworks/Python.framework/Versions/3.5/lib/libpython3.5m.dylib")
# SET(PYTHON_INCLUDE_DIR "/usr/local/Cellar/python3/3.5.0/Frameworks/Python.framework/Versions/3.5/include/python3.5m")

set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6)
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
  endif()
endif()

# Compile with compiler warnings turned on
# if(MSVC)
#   if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
#     string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
#   else()
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
#   endif()
# else()
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
# endif()

include_directories(${PYTHON_INCLUDE_DIR} include)

## include pybing
include_directories(${PROJECT_SOURCE_DIR}/../external/pybind11/include)

## include libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_NANOGUI          "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_CGAL             "Use CGAL"           OFF)
option(LIBIGL_WITH_BOOLEAN          "Use Cork boolean"   OFF)
option(LIBIGL_WITH_COMISO           "Use CoMiso"         ON)
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_LIM              "Use LIM"            ON)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          OFF)
option(LIBIGL_WITH_BBW              "Use BBW"            ON)
option(LIBIGL_WITH_OPENGL_AND_PNG   "Use OpenGL"         ON)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       ON)
option(LIBIGL_WITH_XML              "Use XML"            ON)

if(LIBIGL_WITH_CGAL) # Do not remove or move this block, cgal strange build system fails without it
  find_package(CGAL REQUIRED)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
  include(${CGAL_USE_FILE})
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/../shared/cmake" "libigl")


### Prepare the build environment

include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})

## Optional modules
if (LIBIGL_WITH_VIEWER)
  add_definitions(-DPY_VIEWER)
  list(APPEND SHARED_SOURCES "py_igl_viewer.cpp")
endif ()

if (LIBIGL_WITH_COMISO)
  add_definitions(-DPY_COMISO)
  list(APPEND SHARED_SOURCES "copyleft/py_igl_comiso.cpp")
endif ()


## Prepare the python library
add_library(pyigl SHARED
  python_shared.cpp
  py_vector.cpp
  py_igl.cpp
  py_doc.cpp
  ${SHARED_SOURCES}
)

set_target_properties(pyigl PROPERTIES PREFIX "")
set_target_properties(pyigl PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

target_link_libraries(pyigl ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})

if (WIN32)
  if (MSVC)
    # Enforce size-based optimization and link time code generation on MSVC (~30% smaller binaries in experiments)
    set_target_properties(pyigl PROPERTIES COMPILE_FLAGS "/Os /GL")
    set_target_properties(pyigl PROPERTIES LINK_FLAGS "/LTCG")
  endif()

  # .PYD file extension on Windows
  set_target_properties(pyigl PROPERTIES SUFFIX ".pyd")

  # Link against the Python shared library
  # message(FATAL_ERROR ${PYTHON_LIBRARY})
  # target_link_libraries(igl ${PYTHON_LIBRARY})
  target_link_libraries(pyigl ${PYTHON_LIBRARIES})


elseif (UNIX)
  # It's quite common to have multiple copies of the same Python version
  # installed on one's system. E.g.: one copy from the OS and another copy
  # that's statically linked into an application like Blender or Maya.
  # If we link our plugin library against the OS Python here and import it
  # into Blender or Maya later on, this will cause segfaults when multiple
  # conflicting Python instances are active at the same time.

  # Windows does not seem to be affected by this issue. The solution for Linux
  # and Mac OS is simple: we just don't link against the Python library. The
  # resulting shared library will have missing symbols, but that's perfectly
  # fine -- they will be resolved at import time.

  # .SO file extension on Linux/Mac OS
  set_target_properties(pyigl PROPERTIES SUFFIX ".so")

  #Enable flag if undefined symbols appear on pyigl module import to get notified about the missing symbols at link time
  option(CHECK_UNDEFINED        "Check for undefined symbols"    OFF)

  # Strip unnecessary sections of the binary on Linux/Mac OS 
  if(APPLE)
    set_target_properties(pyigl PROPERTIES MACOSX_RPATH ".")
    
    if (NOT CHECK_UNDEFINED)
      set_target_properties(pyigl PROPERTIES LINK_FLAGS "-undefined dynamic_lookup -dead_strip")
    endif()

    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET pyigl POST_BUILD COMMAND strip -u -r ${CMAKE_CURRENT_BINARY_DIR}/../pyigl.so)
    endif()
  else()

    if (CHECK_UNDEFINED)
      target_link_libraries(pyigl ${PYTHON_LIBRARIES})
      set_target_properties(pyigl PROPERTIES LINK_FLAGS "-Wl,--no-undefined")
    endif()

    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET pyigl POST_BUILD COMMAND strip ${CMAKE_CURRENT_BINARY_DIR}/../pyigl.so)
    endif()
  endif()
endif()
