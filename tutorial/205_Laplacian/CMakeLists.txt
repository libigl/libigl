cmake_minimum_required(VERSION 2.6)
project(205_Laplacian)

add_executable(${PROJECT_NAME}_bin
  main.cpp)
target_include_directories(${PROJECT_NAME}_bin PRIVATE ${LIBIGL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})
