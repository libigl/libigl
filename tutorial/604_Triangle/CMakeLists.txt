cmake_minimum_required(VERSION 2.6)
project(604_Triangle)

include("../CMakeLists.shared")
find_package(TRIANGLE REQUIRED)

include_directories( ${TRIANGLE_INCLUDE_DIR})

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES} ${SHARED_SOURCES} ${TRIANGLE_SOURCES})
target_link_libraries(${CMAKE_PROJECT_NAME} ${SHARED_LIBRARIES})
