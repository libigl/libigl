cmake_minimum_required(VERSION 2.6)
project(607_Picking)

include("../CMakeLists.shared")

find_package(EMBREE REQUIRED)

include_directories(${EMBREE_INCLUDE_DIRS})

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES} ${EMBREE_LIBRARIES})
