cmake_minimum_required(VERSION 2.6)
project(104_Events)

include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(104_Events ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(104_Events ${SHARED_LIBRARIES})
