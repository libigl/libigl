cmake_minimum_required(VERSION 2.6)
project(403_BoundedBiharmonicWeights)

include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})
add_executable(${PROJECT_NAME}_bin
  main.cpp
${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})
