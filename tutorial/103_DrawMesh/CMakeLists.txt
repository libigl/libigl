cmake_minimum_required(VERSION 2.6)
project(103_DrawMesh)

include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(103_DrawMesh ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(103_DrawMesh ${SHARED_LIBRARIES})
