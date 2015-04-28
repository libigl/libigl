cmake_minimum_required(VERSION 2.6)
project(504_NRosyDesign)

include("../CMakeLists.shared")

find_package(LIBCOMISO REQUIRED)

include_directories( ${LIBCOMISO_INCLUDE_DIR}
                     ${LIBCOMISO_INCLUDE_DIRS} )

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES} ${LIBCOMISO_LIBRARY})
