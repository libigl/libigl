cmake_minimum_required(VERSION 2.6)
project(106_Matlab)

include("../CMakeLists.shared")

find_package(MATLAB REQUIRED)
include_directories( ${MATLAB_INCLUDE_DIR} )

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

if(APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "-rpath ${MATLAB_INCLUDE_DIR}/../../bin/maci64"})
endif (APPLE) #APPLE

#message(FATAL_ERROR ${CMAKE_EXE_LINKER_FLAGS})

add_executable(106_Matlab ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(106_Matlab ${SHARED_LIBRARIES} ${MATLAB_LIBRARIES})
