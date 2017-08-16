cmake_minimum_required(VERSION 2.8.12)
project(609_Boolean)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::cgal igl::viewer tutorials)
