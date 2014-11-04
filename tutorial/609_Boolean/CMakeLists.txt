cmake_minimum_required(VERSION 2.6)
project(609_Boolean)

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

# for some reason must come after cgal include
include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES} ${CGAL_LIBRARIES})
