cmake_minimum_required(VERSION 2.6)
project(100_FileIO)

include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(100_FileIO ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(100_FileIO ${SHARED_LIBRARIES})
