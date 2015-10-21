cmake_minimum_required(VERSION 2.6)
project(609_Boolean)

#find_package(CGAL REQUIRED)
#include(${CGAL_USE_FILE})

include_directories(${LIBIGL_INCLUDE_DIRS})
add_executable(${PROJECT_NAME}_bin
  main.cpp
${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})
