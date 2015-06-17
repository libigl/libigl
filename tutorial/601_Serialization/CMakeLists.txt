cmake_minimum_required(VERSION 2.6)
project(601_Serialization)

include("../CMakeLists.shared")

find_package(TINYXML2 REQUIRED)
include_directories( ${TINYXML2_INCLUDE_DIR})

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES} ${TINYXML2_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES})
