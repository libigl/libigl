cmake_minimum_required(VERSION 2.6)
project(505_MIQ)

include("../CMakeLists.shared")

find_package(LIBCOMISO REQUIRED)

include_directories( ${LIBCOMISO_INCLUDE_DIR} )

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME} ${SHARED_LIBRARIES} ${LIBCOMISO_LIBRARY})
