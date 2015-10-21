cmake_minimum_required(VERSION 2.6)
project(608_LIM)

### OpenMP ###
find_package(OpenMP)
if(OPENMP_FOUND AND NOT WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(${LIBIGL_INCLUDE_DIRS})
add_executable(${PROJECT_NAME}_bin
  main.cpp
${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})
