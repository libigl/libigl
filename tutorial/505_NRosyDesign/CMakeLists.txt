cmake_minimum_required(VERSION 2.6)
project(505_NRoSyDesign)

include("../CMakeLists.shared")

find_package(LIBCOMISO REQUIRED)

include_directories( ${LIBCOMISO_INCLUDE_DIR} )

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${CMAKE_PROJECT_NAME} ${SHARED_LIBRARIES} ${LIBCOMISO_LIBRARY})
