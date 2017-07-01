cmake_minimum_required(VERSION 2.8.12)
project(402_PolyharmonicDeformation)

add_executable(${PROJECT_NAME}_bin
  main.cpp)
target_include_directories(${PROJECT_NAME}_bin PRIVATE ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(${PROJECT_NAME}_bin PRIVATE ${LIBIGL_DEFINITIONS})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_VIEWER_EXTRA_LIBRARIES} ${LIBIGL_OPENGL_EXTRA_LIBRARIES} ${LIBIGL_OPENGL_GLFW_EXTRA_LIBRARIES})
