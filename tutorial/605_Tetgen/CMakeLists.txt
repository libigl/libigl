cmake_minimum_required(VERSION 2.6)
project(605_Tetgen)

include("../CMakeLists.shared")
find_package(TETGEN REQUIRED)

include_directories( ${TETGEN_INCLUDE_DIR})

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES} ${SHARED_SOURCES} ${TETGEN_SOURCES})
target_link_libraries(${PROJECT_NAME} ${SHARED_LIBRARIES})
