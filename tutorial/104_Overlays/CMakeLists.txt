cmake_minimum_required(VERSION 2.6)
project(104_Overlays)

include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(104_Overlays ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(104_Overlays ${SHARED_LIBRARIES})
