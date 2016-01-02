cmake_minimum_required(VERSION 2.6)
project(609_Boolean)

find_package(CGAL REQUIRED)
set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
include(${CGAL_USE_FILE})
# CGAL's monkeying with all of the flags. Rather than change the CGAL_USE_FILE
# just get ride of this flag.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  remove_cxx_flag("-stdlib=libc++")
endif()
IF(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") ### Enable parallel compilation for Visual Studio
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
ELSE(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../" )
ENDIF(MSVC)

# This is absolutely necessary for Exact Construction

# for some reason must come after cgal include. I think that it's overwriting
# come flags like CXX_FLAGS
#set(CMAKELISTS_SHARED_INCLUDED FALSE)
include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES} ${CGAL_LIBRARIES})
