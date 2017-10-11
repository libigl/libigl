cmake_minimum_required(VERSION 2.8.12)
project(403_BoundedBiharmonicWeights)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::viewer tutorials)
if(TARGET igl::mosek)
  target_link_libraries(${PROJECT_NAME}_bin igl::mosek)
endif()
