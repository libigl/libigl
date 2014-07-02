cmake_minimum_required(VERSION 2.6)
project(602_Matlab)

include("../CMakeLists.shared")

find_package(MATLAB REQUIRED)
include_directories( ${MATLAB_INCLUDE_DIR} )

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

# Alec: This seems like a clang specific thing and also seems to have an extra
# bracket.
# if(APPLE)
#   set(CMAKE_EXE_LINKER_FLAGS "-rpath ${MATLAB_INCLUDE_DIR}/../../bin/maci64"})
# endif (APPLE) #APPLE

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES} ${MATLAB_LIBRARIES})
