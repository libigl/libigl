cmake_minimum_required(VERSION 2.6)
project(610_CSGTree)

#find_package(CGAL REQUIRED)
#include(${CGAL_USE_FILE})
# This is absolutely necessary for Exact Construction

# for some reason must come after cgal include. I think that it's overwriting
# come flags like CXX_FLAGS
#set(CMAKELISTS_SHARED_INCLUDED FALSE)


include_directories(${LIBIGL_INCLUDE_DIRS})
add_executable(${PROJECT_NAME}_bin
  main.cpp
${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})
