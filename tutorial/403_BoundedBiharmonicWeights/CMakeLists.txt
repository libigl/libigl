cmake_minimum_required(VERSION 2.6)
project(403_BoundedBiharmonicWeights)

include("../CMakeLists.shared")

add_definitions(-DIGL_NO_MOSEK)

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES})
