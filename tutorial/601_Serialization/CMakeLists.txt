cmake_minimum_required(VERSION 2.6)
project(601_Serialization)

include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(101_FileIO ${SOURCES} ${SHARED_SOURCES} ${TINYXML2_SOURCES})
target_link_libraries(101_FileIO ${SHARED_LIBRARIES})
