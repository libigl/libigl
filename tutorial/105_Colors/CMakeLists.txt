cmake_minimum_required(VERSION 2.6)
project(105_Colors)

include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(105_Colors ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(105_Colors ${SHARED_LIBRARIES})
