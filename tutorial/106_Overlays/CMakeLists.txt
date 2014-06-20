cmake_minimum_required(VERSION 2.6)
project(106_Overlays)

include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(106_Overlays ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(106_Overlays ${SHARED_LIBRARIES})
