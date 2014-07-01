cmake_minimum_required(VERSION 2.6)
project(608_LIM)

include("../CMakeLists.shared")
find_package(LIM REQUIRED)

include_directories(${LIM_INCLUDE_DIR})

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES} ${LIM_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES})
