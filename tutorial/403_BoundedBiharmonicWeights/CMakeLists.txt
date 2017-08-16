cmake_minimum_required(VERSION 2.8.12)
project(403_BoundedBiharmonicWeights)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::viewer tutorials)
if(LIBIGL_WITH_MOSEK)
  target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_MOSEK_EXTRA_LIBRARIES})
endif()
